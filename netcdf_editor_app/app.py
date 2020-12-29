from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for
)
from werkzeug.exceptions import abort

from netcdf_editor_app.auth import login_required

bp = Blueprint('app', __name__)


@bp.route('/')
def index():
    return render_template('app/index.html')
