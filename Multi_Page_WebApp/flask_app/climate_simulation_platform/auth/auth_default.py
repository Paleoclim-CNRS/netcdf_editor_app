from climate_simulation_platform.auth import bp

from flask import (
    flash,
    g,
    redirect,
    render_template,
    request,
    session,
    url_for,
)
from werkzeug.security import check_password_hash

from climate_simulation_platform.db import add_user, get_db

class DefaultAuth(object):

    @staticmethod
    def register():
        if request.method == "POST":
            username = request.form["username"]
            password = request.form["password"]
            error = ''

            if not username:
                error += "Username is required."
            elif not password:
                error += "Password is required."

            if len(error) == 0:
                add_user(username, password)
                return redirect(url_for("auth.login"))

        return render_template("auth/register.html")

    @staticmethod
    def login():
        if request.method == "POST":
            username = request.form["username"]
            password = request.form["password"]
            db = get_db()
            error = None
            user = db.execute(
                "SELECT * FROM user WHERE username = ?", (username,)
            ).fetchone()

            if user is None:
                error = "Incorrect username."
            elif not check_password_hash(user["password"], password):
                error = "Incorrect password."

            if error is None:
                session.clear()
                session["user_id"] = user["id"]
                return redirect(url_for("index"))

            flash(error)
        return render_template("auth/login.html")

    @staticmethod
    @bp.before_app_request
    def load_logged_in_user():
        user_id = session.get("user_id")

        if user_id is None:
            g.user = None
        else:
            g.user = (
                get_db().execute("SELECT * FROM user WHERE id = ?", (user_id,)).fetchone()
            )

    @staticmethod
    @bp.route("/logout")
    def logout():
        session.clear()
        return redirect(url_for("index"))