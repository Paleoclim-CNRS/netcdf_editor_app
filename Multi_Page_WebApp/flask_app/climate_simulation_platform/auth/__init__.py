import functools
import os
from climate_simulation_platform.db import get_owner_id

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


def user_required(view):
    @functools.wraps(view)
    def wrapped_view(**kwargs):
        print("kw: ", kwargs, flush=True)
        print(g.user["id"])
        owner_id = get_owner_id(kwargs["_id"])
        print(owner_id)
        if g.user["id"] != owner_id:
            return redirect(url_for("index"))

        return view(**kwargs)

    return wrapped_view
