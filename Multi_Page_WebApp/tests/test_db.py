import sqlite3

import pytest
from netcdf_editor_app.db import get_db, get_file_types, get_file_type_counts


def test_get_close_db(app):
    with app.app_context():
        db = get_db()
        assert db is get_db()

    with pytest.raises(sqlite3.ProgrammingError) as e:
        db.execute("SELECT 1")

    assert "closed" in str(e.value)


def test_init_db_command(runner, monkeypatch):
    class Recorder(object):
        called = False

    def fake_init_db():
        Recorder.called = True

    monkeypatch.setattr("netcdf_editor_app.db.init_db", fake_init_db)
    result = runner.invoke(args=["init-db"])
    assert "Initialized" in result.output
    assert Recorder.called


def test_get_file_type_counts(app):
    with app.app_context():
        file_type_counts = get_file_type_counts(1)
    assert type(file_type_counts) == dict
    for key, val in file_type_counts.items():
        assert type(key) is str
        assert type(val) is int
    assert file_type_counts["raw"] == 2
    assert file_type_counts["routing"] == 1
    assert "fake_file_type" not in file_type_counts.keys()


def test_get_file_types(app):
    with app.app_context():
        file_types = get_file_types(1)
    assert type(file_types) == list
    assert "raw" in file_types
    assert "routing" in file_types
