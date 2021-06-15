from setuptools import setup, find_packages

# Load the README file
with open(file="README.md", mode="r") as readme_handle:
    long_description = readme_handle.read()

setup(
    # Define the library name, this is what is used along with 'pip install'
    name="radarscatter",

    # Define the author of the repository
    author="Dylan Shiltz",
    author_email="djs2112@rit.edu",

    # Define the version of this library: <Major Version>.<Minor Version>.<Maintenance Version>
    version="0.0.1",

    # Here is a small description of the library.  This appears
    # when someone searches for the library on https://pypi.org/search
    description="Implements a variety of microwave radar scattering and soil dielectric models for remote sensing applications",

    # Long description just references the README file from above
    long_description=long_description,
    long_description_content_type="text/markdown",

    # Here is the URL where you can find the code, in this case on GitHub
    url="https://github.com/djshiltz/radarscatter",

    # These are the dependencies I had in place when building this package
    # Newer versions may work, but are not guaranteed
    install_requires=[
        "numpy==1.20.2",
        "matplotlib==3.4.1",
    ],

    # Here are the keywords of my library
    keywords="IEM, SAR, Microwave Remote Sensing",

    # Here are the packages I want to "build"
    packages=find_packages(
        include=["radarscatter", "radarscatter.*"]
    ),

    # # Here we specify any package data, if required
    # package_data={}
    # include_package_data=True,

    # Here I can specify the Python version needed to run this library
    python_requires=">=3.7",

    # # Include any classifiers that will offer more information to
    # # potential users.  For a complete list of classifiers, visit
    # # https://pypi.org/classifiers/
    # classifiers=[],
)