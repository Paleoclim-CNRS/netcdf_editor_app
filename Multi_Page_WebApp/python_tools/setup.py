import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="climpy",
    version="0.0.1",
    license="MIT",
    author="Wesley Banfield",
    author_email="banfield@cerege.fr",
    description="Python tools for Climate Simulation Sciences",
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
    include_package_data=True,
    package_data={
        "": ["bc/ipsl/*.npy", "bc/ipsl/*.nc"],
    },
    zip_safe=False,
    install_requires=[
        "setuptools-git"
        # TODO
    ],
    extras_require={
        "test": ["pytest", "coverage"],
    },
    python_requires=">=3.6",
)
