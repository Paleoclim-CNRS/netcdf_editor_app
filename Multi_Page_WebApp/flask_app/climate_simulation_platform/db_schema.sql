DROP TABLE IF EXISTS user;
DROP TABLE IF EXISTS data_file;

CREATE TABLE user (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  username TEXT UNIQUE NOT NULL,
  password TEXT NOT NULL
);

CREATE TABLE data_files (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  owner_id INTEGER NOT NULL,
  filename TEXT NOT NULL,
  longitude TEXT,
  latitude TEXT,
  info TEXT,
  state TEXT,
  FOREIGN KEY (owner_id) REFERENCES user (id)
);

CREATE TABLE revisions (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  data_file_id INTEGER NOT NULL,
  created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  filepath TEXT NOT NULL,
  revision INTEGER NOT NULL,
  file_type TEXT NOT NULL,
  info TEXT,
  FOREIGN KEY (data_file_id) REFERENCES data_files (id)
);

CREATE TABLE steps (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  data_file_id INTEGER NOT NULL,
  step TEXT NOT NULL,
  parameters TEXT,
  up_to_date INTEGER NOT NULL DEFAULT 1,
  FOREIGN KEY (data_file_id) REFERENCES data_files (id)
);