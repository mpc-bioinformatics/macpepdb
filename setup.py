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
    python_requires='>=3.7',
    install_requires = [
        "psycopg2==2.8.6",
        "SQLAlchemy==1.3.19",
        "alembic==1.4.2",
    ]
)