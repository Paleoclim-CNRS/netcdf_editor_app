from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, current_app
)
from werkzeug.exceptions import abort
from werkzeug.utils import secure_filename

import os
import tempfile

from netcdf_editor_app.auth import login_required
from netcdf_editor_app.db import get_db

bp = Blueprint('app', __name__)


@bp.route('/')
def index():
    db = get_db()
    data_files = db.execute(
        'SELECT created, filename, username'
        ' FROM data_files df JOIN user u ON df.owner_id = u.id'
        ' ORDER BY created DESC'
    ).fetchall()

    return render_template('app/index.html', data_files=data_files)

def allowed_file(filename):
    ALLOWED_EXTENSIONS = {'nc'}
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@bp.route('/upload', methods=('GET', 'POST'))
@login_required
def upload():
    if request.method == 'POST':
        # check if the post request has the file part
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        # if user does not select file, browser also
        # submit an empty part without filename
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            temp_name = next(tempfile._get_candidate_names()) + ".nc"
            file.save(os.path.join(current_app.config['UPLOAD_FOLDER'], temp_name))
            db = get_db()
            db.execute(
                'INSERT INTO data_files (owner_id, filename, filepath)'
                ' VALUES (?, ?, ?)',
                (g.user['id'], filename, temp_name)
            )
            db.commit()
            return redirect(url_for('index'))

    return '''
    <!doctype html>
    <title>Upload new File</title>
    <h1>Upload new File</h1>
    <form method=post enctype=multipart/form-data>
      <input type=file name=file>
      <input type=submit value=Upload>
    </form>
    '''