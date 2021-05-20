from .auth_default import DefaultAuth

from climate_simulation_platform.auth import login_required

class LoggedInAuth(DefaultAuth):

    register = login_required(DefaultAuth.register)