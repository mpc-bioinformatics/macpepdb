# MaCPepDB - Mass Centric Peptide Database

## Description
Creates a peptide databases by digesting proteins stored in FASTA-/Uniprot-Text-files. 

## Development
### Prepare development environment
```bash
# Install necessary dependencies for your distro
sudo pacman -S python python-pip postgresql postgresql-libs

# It is recommended to use a pyenv to make sure the python version is matching
# Follow the instruction to install pyenv on https://github.com/pyenv/pyenv#installation

# Install the correct python version. You can find the needed python version in .python-version at the beginning of the string (.python-version contains the actual name of the python environment).
# The following command will extract the python version from .python-version for you and install it
pyenv install $(cat .python-version | awk 'BEGIN { FS = "/" } ; { print $1 }')

# Create an environment. The following command will extract the python version and environment name from .python-version for you and install it
pyenv virtualenv $(cat .python-version | awk 'BEGIN { FS = "/" } ; { print $1 }') $(cat .python-version | awk 'BEGIN { FS = "/" } ; { print $3 }')

# And activate the environment (later pyenv will do it for you if you visit a folder with a .python-version file)
pyenv activate

# Update pip (to make sure its the newest version)
pip install --upgrade pip

# Install needed python modules
pip install -r ./requirements.txt

# Create necessary folders for your database files in the project folder. This folder are already ignored by GIT.
mkdir -p .db/pgsql/data .db/pgsql/run .db/pgsql/config .db/pgsql/wal
# Initialize Postgresql-database. Use password `developer` to match the existing Procfile for Foreman.
initdb -D .db/pgsql/data -X $(pwd)/.db/pgsql/wal -U postgres -W
# Start the database
postgres -D $(pwd)/.db/pgsql/data -h 127.0.0.1 -p 5433 -k $(pwd)/.db/pgsql/run

# In another terminal(tab) create the databases
psql -h 127.0.0.1 -p 5433 -U postgres -c "create database macpepdb_dev;"
psql -h 127.0.0.1 -p 5433 -U postgres -c "create database macpepdb_test;"
# Run migrations
MACPEPDB_DB_URL=postgresql://postgres:developer@127.0.0.1:5433/macpepdb_test alembic upgrade head
```

If you add some new python modules, make sure you update `requirements.txt` by running `pip freeze > requirements.txt`

### Running tests
```bash
TEST_MACPEPDB_URL=postgresql://postgres:developer@127.0.0.1:5433/macpepdb_test python -m unittest tests/*_test_case.py
```
### Run the modules CLI
Run `python -m macpepdb --help` in the root-folder of the repository.

## Usage

### Native installation
For a native installation you need Python 3.7 or higher and the PostgreSQL headers.
Than update pip with `pip install --upgrade pip` and run `pip install -e git+https://github.com/mpc-bioinformatics/macpepdb.git@<MACPEPDB_GIT_TAG>#egg=MaCPepDB` to install MaCPepDB.
Then you can use MacPepDB by running `python -m macpepdb`. 
Appending `--help` shows the available command line parmeter.

### Docker installation
To create a Docker image use: `docker build --tag macpepdb-py .` . You can use the image to start a container with
`docker run -it --rm macpepdb-py --help`.
To access your files in the container mount your files to `/usr/src/macpepdb/data` with `-v YOUR_DATA_FOLDER:/usr/src/macpepdb/data` (add it before the `macpepdb-py`). Keep in mind your working in a container, so all files pathes are within the container.   
If you intend to create a protein/peptide database and your Postgresql server is running in a Docke container too, make sure both, the  Postgresql server and the MacPepDB container have access to the same Docker network by adding `--network=YOUR_DOCKER_NETWORK` (before the ´macpepdb-py´).

### Building a databse
#### Prepare the database
Run `MACPEPDB_DB_URL=postgresql://<USER>:<PASSWORD>@<HOST>:<PORT>/<DATABASE> test alembic upgrade head`    
If you use the docker container, run the command in a temporary container: `docker run --rm -it macpepdb-py sh`

#### Fill the database
First create a work folder with the following structure:
```
|_ work_dir
   |__ protein_data
   |__ taxonomy_data
   |__ logs
```
Place your protein data files as `.dat`- or `.txt`-files, containing the proteins in UniProt-text-format, in the `protein_data`-folder.
If you like to use the web interface as well, download the `taxdump.zip` from [NCBI](https://ftp.ncbi.nih.gov/pub/taxonomy/) and put the contained `.dmp`-files in the `taxonomy_data`-folder.

Than start the database maintenance job with `python -m macpepdb database ...`. Run `python -m macpepdb database --help` to see the required arguments. Remember to use the container internal paths when using a docker container.