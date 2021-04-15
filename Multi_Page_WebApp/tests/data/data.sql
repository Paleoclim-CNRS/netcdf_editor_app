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
  (2, 1, 'tmp2.tmp', 1, 'routing'),
  (3, 1, 'tmp1_1.tmp', 2, 'raw');

INSERT INTO steps (id, data_file_id, step, parameters, up_to_date)
VALUES
  (1,	9,	'regrid',	'{"id": 9, "limits": "default", "Longitude Step": "1", "Latitude Step": "1", "interpolator": "nearest"}',	1),
  (2,	9,	'routing',	'{"id": 9, "topo_var": "Z"}',	1),
  (3,	2,	'regrid',	'{"id": 2, "limits": "default", "Longitude Step": "1", "Latitude Step": "1", "interpolator": "nearest"}', 1),
  (4,	2,	'routing',	'{"id": 2, "topo_var": "relief"}', 1),
  (5,	2,	'pft',	'{"id": 2, "data": "{\"dataArray\":[[15,0,75,25,0,0,0,0,0,0,0,0,0,0],[35,0,15,55,0,0,0,0,0,0,30,0,0,0],[50,0,0,0,0,70,30,0,0,0,0,0,0,0],[80,0,0,0,0,40,30,0,30,0,0,0,0,0],[90,0,0,0,0,0,30,0,0,40,30,0,0,0]]}"}',	1),
  (6,	2,	'heatflow',	'{"id": 2, "invalidated": "yes", "has_params": "no"}', 1),
  (7,	2,	'ahmcoef',	'{"id": 2, "invalidated": "yes", "has_params": "no"}', 1);