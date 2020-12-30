import os

from flask import Flask, session

from werkzeug.exceptions import abort


def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'netcdf_editor.sqlite'),
        UPLOAD_FOLDER=app.instance_path,
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # a simple page that says hello
    @app.route('/hello')
    def hello():
        return 'Hello, World! from NetCDF Editor App'

    @app.route('/session')
    def debug_session():
        if app.env == 'development':
            return str(dict(session))
        return abort(403, "Session debug only available in debug mode")

    from . import db
    db.init_app(app)

    from . import auth
    app.register_blueprint(auth.bp)

    from . import app as editor_app
    app.register_blueprint(editor_app.bp)
    app.add_url_rule('/', endpoint='index')

    return app
