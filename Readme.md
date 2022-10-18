# MaCPepDB 2.0 - Mass Centric Peptide Database

## Description

Creates a peptide databases by digesting proteins stored in FASTA-/Uniprot-Text-files.

## Dependencies

Only necessary for development and non-Docker installation

- GIT
- Build tools (Ubuntu: `build-essential`, Arch Linux: `base-devel`)
- C/C++-header for PostgreSQL (Ubuntu: `libpq-dev`, Arch Linux: `postgresql-libs`)
- C/C++-header for libev (Ubuntu: `libev-dev`, Arch Linux: `libev`)
- Rust Compiler
- Docker & Docker Compose
- Python 3.x
- [pyenv](https://github.com/pyenv/pyenv)
- [pipenv](https://pipenv.pypa.io/en/latest/)

## Development

Make sure `pipenv` finds `pyenv`

### Prepare development environment

```bash
# Install correct python version and create environment
pipenv install -d

# Change to environment
pipenv shell

# Start the database
docker-compose up

# Run migrations
MACPEPDB_URL=postgresql://postgres:developer@127.0.0.1:5433/macpepdb_dev pipenv run db-migrate
```

Use `pipenv` to install or uninstall Python modules

### Running tests

```bash
TEST_MACPEPDB_URL=postgresql://postgres:developer@127.0.0.1:5433/macpepdb_dev pipenv run tests
```

### Run the modules CLI

Run `python -m macpepdb --help` in the root-folder of the repository.

## Usage

### Native installation

Than update pip with `pip install --upgrade pip` and run `pip install -e git+https://github.com/mpc-bioinformatics/macpepdb.git@<MACPEPDB_GIT_TAG>#egg=MaCPepDB` to install MaCPepDB.
Then you can use MacPepDB by running `python -m macpepdb`.
Appending `--help` shows the available command line parameter.

### Docker installation

To create a Docker image use: `docker build --tag macpepdb-py .` . You can use the image to start a container with
`docker run -it --rm macpepdb-py --help`.
To access your files in the container mount your files to `/usr/src/macpepdb/data` with `-v YOUR_DATA_FOLDER:/usr/src/macpepdb/data` (add it before the `macpepdb-py`). Keep in mind your working in a container, so all files pathes are within the container.  
If you intend to create a protein/peptide database and your Postgresql server is running in a Docke container too, make sure both, the Postgresql server and the MacPepDB container have access to the same Docker network by adding `--network=YOUR_DOCKER_NETWORK` (before the ´macpepdb-py´).

### Building a database

#### Prepare the database

1. Follow the [Citus documentation](http://docs.citusdata.com/en/v10.0/installation/multi_node.html) to setup a Citus cluster.
2. Run `psql -h <CITUS_CONTROLLER> -U <DB_USER> -c "ALTER DATABASE <DB_NAME> SET citus.multi_shard_modify_mode = 'sequential';"` and `psql -h <CITUS_CONTROLLER> -U <DB_USER> -c "ALTER DATABASE <DB_NAME> SET citus.shard_count = 100;"` to configure the database
3. Run `MACPEPDB_DB_URL=postgresql://<USER>:<PASSWORD>@<HOST>:<PORT>/<DATABASE> alembic upgrade head`, if you use the docker container, run the command in a temporary container: `docker run --rm -it macpepdb sh`

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

For example:

```
python -m macpepdb database --work-dir workdir -t 2 -e trypsin -d postgresql://postgres:developer@127.0.0.1:5433/macpepdb_dev
```

### WebAPI

Create a new config file with the default config

```
python -m macpepdb web write-config-file <PATH_TO_CONFIG_YAML>
```

Adjust the YAML file to your needs. Than start the WebAPI with

```
python -m macpepdb web serve -e production -c <PATH_TO_CONFIG_YAML>
```

For high availability in production use start multiple WebAPI and combine them with NginX (have a look in `nginx.example.conf`)

## Upgrading

### 1.x to 2.x

Due to changes of the database schema and the database engine, version 2.x is not compatible with version 1.x. You have to recreate the database.

## Citation and Publication

- **MaCPepDB: A Database to Quickly Access All Tryptic Peptides of the UniProtKB**  
   Julian Uszkoreit, Dirk Winkelhardt, Katalin Barkovits, Maximilian Wulf, Sascha Roocke, Katrin Marcus, and Martin Eisenacher  
   Journal of Proteome Research 2021 20 (4), 2145-2150  
   [DOI: 10.1021/acs.jproteome.0c00967](https://doi.org/10.1021/acs.jproteome.0c00967)
