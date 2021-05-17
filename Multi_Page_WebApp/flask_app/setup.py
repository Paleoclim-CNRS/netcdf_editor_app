import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="csp",
    version="1.0.0",
    license='MIT',
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
        'Development Status :: 3 - Alpha',
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    install_requires=[
        #TODO
        ],
    extras_require={
        "test":  ["pytest", "coverage"],
    },
    python_requires=">=3.6",
    zip_safe=False,
    include_package_data=True,
    package_data={
      '': ['templates/*/*.html', 'templates/*.html', 'static/*', 'db_schema.sql'],
   },
)