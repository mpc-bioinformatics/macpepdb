import setuptools

with open("Readme.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MaCPepDB",
    version="1.2.3",
    author="Dirk Winkelhardt",
    author_email="dirk.winkelhardt@rub.de",
    description="Package for building large peptide databases.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mpc-bioinformatics/macpepdb",
    packages=setuptools.find_packages(),
    classifiers=[
    ],
    python_requires='>=3.9, <3.10',
    install_requires = [
        "alembic >=1, <2",
        "bjoern >= 3, <4",
        "Flask >=1, <2",
        "markupsafe ==2.0.1",       # fix for `ImportError: cannot import name 'soft_unicode'
        "matplotlib >=3, <4",
        "piwikapi ==0.3",
        "PyYAML >=5, <6",
        "psycopg2 >=2, <3"
    ]
)