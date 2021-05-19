import setuptools

# Get long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Get version string embedded inside package
import re

VERSIONFILE = "climate_simulation_platform/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setuptools.setup(
    name="csp",
    version=verstr,
    license="MIT",
    author="Wesley Banfield",
    author_email="banfield@cerege.fr",
    description="Flask app for Climate simulations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CEREGE-CL/netcdf_editor_app/",
    project_urls={
        "Bug Tracker": "https://github.com/CEREGE-CL/netcdf_editor_app/issues",
        "Documentation": "https://cerege-cl.github.io/netcdf_editor_app/",
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    install_requires=[
        # TODO
    ],
    extras_require={
        "test": ["pytest", "coverage"],
    },
    python_requires=">=3.6",
    zip_safe=False,
    include_package_data=True,
    package_data={
        "": ["templates/*/*.html", "templates/*.html", "static/*", "db_schema.sql"],
    },
)
