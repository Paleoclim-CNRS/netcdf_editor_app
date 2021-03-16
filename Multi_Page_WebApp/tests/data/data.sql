INSERT INTO user (username, password)
VALUES
  ('test', 'pbkdf2:sha256:50000$TCI4GzcX$0de171a4f4dac32e3364c7ddc7c14f3e2fa61f2d17574483f7ffbb431b4acb2f'),
  ('other', 'pbkdf2:sha256:50000$kJPKsz6N$d2d4784f1b030a9761f5ccaeeaca413f27f2ecb76d6168407af962ddce849f79');

INSERT INTO data_files (id, owner_id, filename, longitude, latitude)
VALUES
  (1, 1, 'filename', 'lon', 'lat');

INSERT INTO revisions (id, data_file_id, filepath, revision, file_type)
VALUES
  (1, 1, 'tmp1.tmp', 0, 'raw'),
  (2, 1, 'tmp2.tmp', 1, 'routing');