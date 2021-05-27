import sqlite3
import os
import json
import tempfile

import click
from climate_simulation_platform.constants import dependant_files
from flask import current_app, g
from flask.cli import with_appcontext
from flask import flash

from werkzeug.utils import secure_filename

import xarray as xr

from werkzeug.security import generate_password_hash


def get_db():
    if "db" not in g:
        g.db = sqlite3.connect(
            current_app.config["DATABASE"], detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row

    return g.db


def close_db(e=None):
    db = g.pop("db", None)

    if db is not None:
        db.close()


def init_db():
    db = get_db()

    with current_app.open_resource("db_schema.sql") as f:
        db.executescript(f.read().decode("utf8"))


def update_user_password(user_id, password):
    print(f"Updating password for userid {user_id}", flush=True)
    db = get_db()
    db.execute(
        "UPDATE user SET password = ? WHERE id = ?",
        (generate_password_hash(password), user_id),
    )
    db.commit()


def add_user(username, password, init=False):
    db = get_db()

    try:
        row = db.execute(
            "SELECT id FROM user WHERE username = ?", (username,)
        ).fetchone()
        _id = row["id"]
    except (IndexError, TypeError):
        _id = None
    # The user is already in the database
    if _id is not None:
        # Check Updating own name
        if init:
            update_user_password(_id, password)
        elif "id" in g.user and g.user["id"] == _id:
            update_user_password(_id, password)
        else:
            flash(
                "User already in database and you don't have permission to change password"
            )

    else:
        db.execute(
            "INSERT INTO user (username, password) VALUES (?, ?)",
            (username, generate_password_hash(password)),
        )
    db.commit()


def load_file(_id, file_type=None, revision=-1):
    # Get filename
    filepath = get_file_path(_id, file_type=file_type, revision=revision)
    if filepath is None:
        return None
    # Load file
    try:
        ds = xr.open_dataset(filepath)
    except ValueError:
        ds = xr.open_dataset(filepath, decode_times=False)
    return ds


def overwrite_file_nc(_id, ds, file_type=None, revision=-1):
    # Get filename
    filepath = get_file_path(_id, file_type=file_type, revision=revision)
    if filepath is None:
        return None
    ds.to_netcdf(
        filepath,
        format="NETCDF3_64BIT",
    )


def get_coord_names(_id):
    ds = load_file(_id)
    return [name for name in ds.coords]


def set_data_file_coords(_id, longitude, latitude):
    db = get_db()
    db.execute(
        "UPDATE data_files SET longitude = ?, latitude = ? WHERE id = ?",
        (longitude, latitude, str(_id)),
    )
    db.commit()


def is_data_file(_id, file_type):
    name, extension = os.path.splitext(get_file_path(_id, file_type, full=False))
    if extension in [".nc"]:
        return True
    return False


def upload_file(file, data_file_id=None, file_type="raw", info={}):
    if file_type == "raw" and data_file_id is not None:
        raise AttributeError("Cannot upload a raw file to a different datafile")
    filename = secure_filename(file.filename)
    name, extension = os.path.splitext(filename)
    name = "_".join(name.split("."))
    filename = name + extension
    temp_name = next(tempfile._get_candidate_names()) + ".nc"
    # Save the file to the file system
    file.save(os.path.join(current_app.config["UPLOAD_FOLDER"], temp_name))
    # Add file to known data_files
    db = get_db()
    if data_file_id is None:
        data_file = db.execute(
            "INSERT INTO data_files (owner_id, filename, info)" " VALUES (?, ?, ?)",
            (g.user["id"], filename, json.dumps(info)),
        )
        db.commit()
        data_file_id = data_file.lastrowid
        revision = 0
    else:
        # Get revision
        revision = (
            db.execute(
                "SELECT MAX(revision) FROM revisions WHERE data_file_id = ? ORDER BY revision ASC",
                (data_file_id,),
            ).fetchone()["MAX(revision)"]
            + 1
        )
    # ADD the file to the revisions table
    db.execute(
        "INSERT INTO revisions (data_file_id, filepath, revision, file_type)"
        " VALUES (?, ?, ?, ?)",
        (data_file_id, temp_name, revision, file_type),
    )
    db.commit()

    return data_file_id


def add_info(_id, file_type, info={}, revision=-1):
    # If the file isnt a data file but a tar for example then return
    if not is_data_file(_id, file_type):
        return
    ds = load_file(_id, file_type, revision)
    for key, value in info.items():
        key = "climate_sim_" + key
        if key in ds.attrs:
            if str(value) not in ds.attrs[key]:
                ds.attrs[key] += "\n" + str(value)
        else:
            ds.attrs[key] = str(value)
    overwrite_file_nc(_id, ds, file_type, revision)


def get_info(_id, file_type):
    db = get_db()

    # Get initial datafile info
    query = "SELECT info FROM data_files WHERE id = ?"
    info = json.loads(db.execute(query, (str(_id),)).fetchone()["info"])

    # Get specific file modification info
    file_types_to_analyse = dependant_files(file_type)
    for file_type in file_types_to_analyse:
        query = "SELECT revision, info FROM revisions WHERE data_file_id = ? AND file_type = ? ORDER BY revision"
        results = db.execute(query, (str(_id), file_type)).fetchall()
        new_info = []
        for i in range(len(results)):
            res = results[i]
            if res["info"] is not None and len(res["info"]):
                new_info.append([i + 1, json.loads(res["info"])])
        for item in new_info:
            rev, d = item
            for key, item in d.items():
                key = key + "_" + file_type
                item = f"Revision {rev}: " + item
                if key in info.keys():
                    if item not in info[key]:
                        info[key] += "\n" + item
                else:
                    info[key] = item
    return info


def save_revision(_id, ds, file_type, info={}):
    temp_name = next(tempfile._get_candidate_names()) + ".nc"
    # Save the file to the file system
    # We use NETCDF3_64Bit files because on some of the calculators the later libraries have not been installed
    # to handle netcdf4(?) or 5(this is sure)
    # However Anta Sarr ran a simulation using netcdf5 files so the simulator appears to accept the files
    ds.to_netcdf(
        os.path.join(current_app.config["UPLOAD_FOLDER"], temp_name),
        format="NETCDF3_64BIT",
    )

    save_file_to_db(_id, temp_name, file_type, info)


def save_file_to_db(_id, filename, file_type, info={}):
    # ADD the file to the revisions table
    db = get_db()
    # Get latest revision
    query = "SELECT revision FROM revisions WHERE data_file_id = ? ORDER BY revision DESC LIMIT 0, 1"
    revision_nb = db.execute(query, (str(_id),)).fetchone()["revision"] + 1
    db.execute(
        "INSERT INTO revisions (data_file_id, filepath, revision, file_type, info)"
        " VALUES (?, ?, ?, ?, ?)",
        (str(_id), filename, revision_nb, file_type, json.dumps(info)),
    )
    db.commit()
    try:
        flash(f"Saved new version of {file_type}")
    except RuntimeError:
        pass


def step_seen(_id, step):
    db = get_db()

    query = "SELECT COUNT(id) from steps WHERE data_file_id = ? AND step = ?"
    results = db.execute(query, (str(_id), step)).fetchone()["COUNT(id)"]

    return results > 0


def step_up_to_date(_id, step):
    db = get_db()

    query = "SELECT up_to_date from steps WHERE data_file_id = ? AND step = ?"
    results = db.execute(query, (str(_id), step)).fetchone()["up_to_date"]

    # Up to date if equal to 1
    return bool(results)


def steps_seen(_id):
    db = get_db()

    query = "SELECT step from steps WHERE data_file_id = ?"
    results = db.execute(query, (str(_id),)).fetchall()
    steps = [res["step"] for res in results]

    return steps


def step_parameters(_id, step):
    db = get_db()

    query = "SELECT parameters from steps WHERE data_file_id = ? AND step = ?"
    parameters = db.execute(query, (str(_id), step)).fetchone()
    if len(parameters):
        parameters = json.loads(parameters["parameters"])

    return parameters


def save_step(_id, step, parameters, up_to_date=True):
    if type(parameters) is dict:
        parameters = json.dumps(parameters)

    db = get_db()
    # Get already seen steps
    query = "SELECT step FROM steps WHERE data_file_id = ?"
    steps = [s["step"] for s in db.execute(query, (str(_id),)).fetchall()]
    if step not in steps:
        db.execute(
            "INSERT INTO steps (data_file_id, step, parameters, up_to_date)"
            " VALUES (?, ?, ?, ?)",
            (str(_id), step, parameters, int(up_to_date)),
        )
        print(f" [*] Saving : {(str(_id), step, parameters, int(up_to_date))}")
    else:  # update parameters and set as up to date
        db.execute(
            "UPDATE steps SET parameters = ?, up_to_date = ? WHERE data_file_id = ? AND step = ?",
            (parameters, int(up_to_date), str(_id), step),
        )
        print(f" [*] Saving : {(parameters, int(up_to_date), str(_id), step)}")
    db.commit()
    try:
        flash(f"Updated Step {step} for {_id}")
    except RuntimeError:
        pass


def invalidate_step(_id, step):
    db = get_db()

    db.execute(
        "UPDATE steps SET up_to_date = ? WHERE data_file_id = ? AND step = ?",
        (0, str(_id), step),
    )
    print(f" [*] Invalidating : {(0, str(_id), step)}")

    db.commit()

    try:
        flash(f"Invalidated Step {step} for {_id}")
    except RuntimeError:
        pass


def get_owner_id(data_file_id):
    db = get_db()

    query = "SELECT owner_id FROM data_files WHERE id = ?"
    owner_id = db.execute(query, (data_file_id,)).fetchone()["owner_id"]
    return owner_id


def get_latest_file_versions():
    query = (
        "SELECT created, filename, df.id FROM"
        + " (SELECT MAX(revision) as revision, created, data_file_id FROM revisions GROUP BY data_file_id) as r"
        + " JOIN (SELECT * FROM data_files WHERE owner_id = ?) AS df ON r.data_file_id = df.id"
        + " JOIN user u ON df.owner_id = u.id"
        + " ORDER BY created DESC"
    )
    db = get_db()
    data_files = db.execute(query, (g.user["id"],)).fetchall()
    return data_files


def get_file_type_counts(_id):
    db = get_db()
    revisions = db.execute(
        "SELECT file_type, count(file_type) \
        FROM revisions WHERE data_file_id = ? GROUP BY file_type;",
        (str(_id),),
    ).fetchall()
    file_type_counts = {}
    for r in revisions:
        file_type_counts[r[0]] = r[1]
    return file_type_counts


def get_file_types(_id):
    db = get_db()
    file_types = db.execute(
        "SELECT DISTINCT file_type FROM revisions WHERE data_file_id = ?", (str(_id),)
    ).fetchall()

    file_types = [ft["file_type"] for ft in file_types]
    return file_types


def get_filename(_id):
    query = "SELECT filename FROM data_files WHERE id = ?"
    db = get_db()
    return db.execute(query, (str(_id),)).fetchone()["filename"]


def get_file_path(_id, file_type=None, full=True, revision=-1):
    db = get_db()
    if file_type is None:
        revisions = db.execute(
            "SELECT revision FROM revisions WHERE data_file_id = ? ORDER BY revision ASC",
            (str(_id),),
        ).fetchall()
    else:
        revisions = db.execute(
            "SELECT revision FROM revisions WHERE data_file_id = ? AND file_type = ? ORDER BY revision ASC",
            (str(_id), file_type),
        ).fetchall()
    if not len(revisions):
        return None
    revisions = [rev["revision"] for rev in revisions]
    revision_nb = revisions[revision]
    filepath = db.execute(
        "SELECT filepath FROM revisions WHERE data_file_id = ? AND revision = ?",
        (
            str(_id),
            str(revision_nb),
        ),
    ).fetchone()["filepath"]
    if not full:
        return filepath
    return os.path.join(current_app.instance_path, filepath)


def get_all_file_paths(_id, full=True):
    db = get_db()
    filepath = db.execute(
        "SELECT filepath FROM revisions WHERE data_file_id = ?", (str(_id),)
    ).fetchall()
    filepath = [fp["filepath"] for fp in filepath]
    if not full:
        return filepath
    return [os.path.join(current_app.instance_path, fp) for fp in filepath]


def remove_data_file(_id):
    # Get filepath
    file_paths = get_all_file_paths(_id)
    for filepath in file_paths:
        os.remove(filepath)
    # Remove entry from DB
    db = get_db()
    db.execute("DELETE FROM data_files WHERE id = ?", (str(_id),))
    db.execute("DELETE FROM revisions WHERE data_file_id = ?", (str(_id),))
    db.commit()


def get_lon_lat_names(_id):
    db = get_db()
    return db.execute(
        "SELECT longitude, latitude FROM data_files WHERE id = ?", (str(_id),)
    ).fetchone()


@click.command("init-db")
@with_appcontext
def init_db_command():
    """Clear the existing data and create new tables."""
    init_db()
    click.echo("Initialized the database.")


@click.command("create-user")
@click.argument("username")
@click.argument("password")
@with_appcontext
def create_user_command(username, password):
    add_user(username, password)


def init_app(app):
    app.teardown_appcontext(close_db)
    app.cli.add_command(init_db_command)
    app.cli.add_command(create_user_command)
