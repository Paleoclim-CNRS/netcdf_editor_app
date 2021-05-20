import os

from flask import Flask, session
from flask.templating import render_template

from werkzeug.exceptions import abort

from ._version import __version__


def create_app(test_config=None):
    # create and configure the app
    app = Flask(
        __name__, instance_relative_config=True, instance_path="/usr/src/app/instance"
    )
    app.config.from_mapping(
        SECRET_KEY=b'\xd8\xb7\xc5 \xdc\xac\x92\xa6\xfd"\xc2a\xe4k*\x17',
        DATABASE=os.path.join(app.instance_path, "netcdf_editor.sqlite"),
        UPLOAD_FOLDER=app.instance_path,
        AUTH = os.environ.get("AUTH", 'basic').lower(),
        THANKS = os.environ.get("THANKS", "")
    )

    if app.config['AUTH'] == 'logged_in':
        with app.app_context():
            from .db import add_user
            import random
            import string
            password = ''.join(random.choice(string.ascii_lowercase) for i in range(10))
            username = os.environ.get("CSP_USERNAME", 'admin')
            password = os.environ.get("CSP_PASSWORD", password)
            print("Username: ", username, flush=True)
            print("Password: ", password, flush=True)
            add_user(username, password)
        

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile("config.py", silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # a simple page that says hello
    @app.route("/hello")
    def hello():
        return "Hello, World! from The Climate Simulation Platform"

    @app.route("/about")
    def about():
        return render_template("app/about.html")

    @app.route("/session")
    def debug_session():
        if app.env == "development":
            return str(dict(session))
        return abort(403, "Session debug only available in debug mode")

    from . import db

    db.init_app(app)

    from . import auth

    app.register_blueprint(auth.bp)

    from . import message_broker

    message_broker.init_app(app)

    from . import app as editor_app

    app.register_blueprint(editor_app.bp)
    app.add_url_rule("/", endpoint="index")

    return app
