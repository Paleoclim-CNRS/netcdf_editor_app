from datetime import datetime
import os
import shutil
import subprocess
import tempfile

import mosaix

from climate_simulation_platform import create_app

from climate_simulation_platform.db import get_file_path, save_file_to_db

def calculate_weights_mosaix(body):
    print(f"{datetime.now()} Calculating weights mosaix", flush=True)

    app = create_app()

    _id = body["_id"]

    # Load file
    with app.app_context():
        bathy_file = get_file_path(_id, "bathy", full=True)
        coords_file = get_file_path(_id, "weight_coords_mosaix", full=True)

    runner = mosaix.MosaixRunner(bathy_file, coords_file)
    runner.run()

    # Tar files directly into directory
    temp_name = runner.output_dir + ".tar.gz"
    temp_path = os.path.join(app.config["UPLOAD_FOLDER"], temp_name)

    print(f"{datetime.now()} Compressing files to tar", flush=True)
    subprocess.Popen(
        ["tar", "-cvzf", temp_path, os.path.join(runner.mosaix_dir, runner.temp_dir)],
        cwd=runner.mosaix_dir,
    ).wait()

    # Add file to db
    print(f"{datetime.now()} Saving new db file {temp_name} to database", flush=True)
    with app.app_context():
        save_file_to_db(_id, temp_name, "weights_mosaix")

    # Delete Folder
    print(f"{datetime.now()} Cleaning up")
    runner.cleanup()
