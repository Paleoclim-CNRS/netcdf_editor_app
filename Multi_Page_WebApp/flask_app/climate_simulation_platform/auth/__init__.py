import functools
import os

from flask import Blueprint, redirect, url_for, g

# Module wide variables
bp = Blueprint("auth", __name__, url_prefix="/auth")

def login_required(view):
    @functools.wraps(view)
    def wrapped_view(**kwargs):
        if g.user is None:
            return redirect(url_for("auth.login"))

        return view(**kwargs)

    return wrapped_view

# Load correct auth based of from enviroment varaible

#TODO in csp.__init__ this is already added to the flask config
# we should probably use the same everywhere
auth_config = os.environ.get("AUTH", "basic").lower()
if auth_config == "basic":
    from .auth_default import DefaultAuth as obj
elif auth_config == 'logged_in':
    from .auth_logged_in import LoggedInAuth as obj

# Set global
@bp.route("/register", methods=("GET", "POST"))
def register():
    return obj.register()

@bp.route("/login", methods=("GET", "POST"))
def login():
    return obj.login()

load_logged_in_user = obj.load_logged_in_user

logout = obj.logout