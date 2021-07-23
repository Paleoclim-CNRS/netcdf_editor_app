from datetime import datetime
import mosaic

import os
import uuid
import tempfile
import shutil
import subprocess

from climate_simulation_platform import create_app
from climate_simulation_platform.db import get_file_path, save_file_to_db


def calculate_weights(body):
    print(f"{datetime.now()} Calculating weights mosaic", flush=True)
    app = create_app()

    _id = body["id"]

    # Load file
    with app.app_context():
        bathy_file = get_file_path(_id, "bathy", full=True)
        coords_file = get_file_path(_id, "weight_coords", full=True)
        subbasins_file = get_file_path(_id, "sub_basins", full=True)

    print(f"{datetime.now()} Running Mosaic Runner", flush=True)
    runner = mosaic.MosaicRunner(
        root="/tmp", cpl_dir="/usr/src", user_name=str(uuid.uuid4())
    )
    runner.run(
        bathy_file=bathy_file, coords_file=coords_file, subbasins_file=subbasins_file
    )

    # Tar files directly into directory
    temp_name = next(tempfile._get_candidate_names()) + ".tar.gz"
    temp_path = os.path.join(app.config["UPLOAD_FOLDER"], temp_name)

    print(f"{datetime.now()} Compressing files to tar", flush=True)
    subprocess.Popen(
        ["tar", "-cvzf", temp_path, "IGCM", "DOMSK"],
        cwd=os.path.join("/", "home", runner.user_name),
    ).wait()
    # Add file to db
    print(f"{datetime.now()} Saving new db file {temp_name} to database", flush=True)
    with app.app_context():
        save_file_to_db(_id, temp_name, "weights")

    # Delete Folder
    print(f"{datetime.now()} Cleaning up")
    shutil.rmtree(os.path.join("/", "home", runner.user_name, "IGCM"))
    shutil.rmtree(os.path.join("/", "home", runner.user_name, "DOMSK"))
    runner.cleanup()
