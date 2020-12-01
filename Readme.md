# TrypperDB

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
psql -h 127.0.0.1 -p 5433 -U postgres -c "create database trypperdb_dev;"
psql -h 127.0.0.1 -p 5433 -U postgres -c "create database trypperdb_test;"
# Run migrations
TRYPPERDB_DB_URL=postgresql://postgres:developer@127.0.0.1:5434/trypperdb_test alembic upgrade head
```

If you add some new python modules, make sure you update `requirements.txt` by running `pip freeze > requirements.txt`

### Running tests
```bash
TRYPPERDB_TEST_DB_URL=postgresql://postgres:developer@127.0.0.1:5433/trypperdb_test python -m unittest tests/*_test_case.py
```

## Usage

### Native installation
You can follow the development instruction to `pip install -r ./requirements.txt`. Then you can use TrypperDB by running `python -m trypperdb`. 
Appending `--help` shows the existing command line parmeter.

### Docker installation
To create a Docker image use: `docker build --tag trypperdb-py .` . You can use the image to start a container with
`docker run -it --rm trypperdb-py --help`.
To access your files in the container mount your files to `/usr/src/trypperdb/data` with `-v YOUR_DATA_FOLDER:/usr/src/trypperdb/data` (add it before the `trypperdb-py`). Keep in mind your working in a container, so all files pathes are within the container.   
If you intend to create a protein/peptide database and your Postgresql server is running in a Docke container too, make sure both, the  Postgresql server and the TrypperDB container have access to the same Docker network by adding `--network=YOUR_DOCKER_NETWORK` (before the ´trypperdb-py´).

### Building a databse
#### Prepare the database
Run `TRYPPERDB_DB_URL=postgresql://<USER>:<PASSWORD>@<HOST>:<PORT>/<DATABASE> test alembic upgrade head`    
If you use the docker container, run the command in a temporary container: `docker run --rm -it trypperdb-py sh`

#### Digest
To build the database you need a running postgresql server with an empty database. As input file for the digestion you can use FASTA-files or the [UniProt-text-files](https://www.uniprot.org/docs/userman.htm#linetypes). Then run the digestion command: `python -m trypperdb digestion ...`, check out `python -m trypperdb digestion --help` for more information.

#### Include taxonomy trees
This step is only necessary if you want to run the webinterface. Download the `taxdump.zip` from [NCBI](https://ftp.ncbi.nih.gov/pub/taxonomy/). Then run the include command: `python -m trypperdb taxonomy-maintenance ...`, check out `python -m trypperdb taxonomy-maintenance --help` for more information.
