![FunCoup Logo](./FunCoup_logo.png)

[FunCoup 6](https://funcoup6.scilifelab.se) is a functional association network tool developed using Python (v3.10.8) and the Django (v4.1.1) framework. It processes biological data to train networks and visualize results in a web application. The data is stored in a PostgreSQL database, and the frontend leverages Bootstrap and D3.js to provide an interactive user experience. FunCoup 6 is containerized and deployed using Docker and Docker Compose for seamless development and production environments.

## Authors

- Davide Buzzao ([GitHub](https://github.com/davidebuzzao))
- Emma Persson ([GitHub](https://github.com/emmape))

## Directory Structure

- **FunCoup/data/**: Contains the backend components, including `dataCollection` and `dataTraining`. These modules are responsible for collecting and scoring the data, and for training the networks using bin-free redundancy weighting and naive Bayesian statistics. `dataQuery` contains scripts to query the database. `networkAnalysis` contains scripts to get pathway/tissue gene annotations, as well as performing orthology-transfering of whole networks for 618 [InparanoiDB 9](https://inparanoidb.sbc.su.se) species.
- **FunCoup/website/**: Contains the frontend components, built using Python, JavaScript, HTML, and CSS. It utilizes Bootstrap components for UI and D3.js for drawing the network for each query.
- **FunCoup/config/**: Stores configuration files for species, gold standards, evidences, and parameters. These configuration files are specific to each database and website instance.
- **results/**: Stores the results of the analysis performed for the paper, which is currently available as a preprint on bioRxiv at [https://doi.org/10.1101/2024.09.13.612391](https://doi.org/10.1101/2024.09.13.612391). To be able to download this folder, install Git LFS by following the installation instructions for your operating system: [Installing Git LFS](https://git-lfs.github.com/). Descriptive statistics and benchmarking are performed using R (v4.3.2). To reproduce this analysis, a Conda environment `r4_env` needs to be set up by running `conda env create -f results/environment.yml` and `conda activate r4_env`. ***Please, refer to this [bitbucket repo](https://bitbucket.org/sonnhammergroup/) to download all data and results if you want to reproduce the results***

For security reasons, database tokens are not included in this repository.

## Installation

### Step 1: Install Dependencies

FunCoup is implemented in Python and tested within a version-controlled Conda environment. To set up the environment, clone the FunCoup repository, and install [Miniconda](https://docs.conda.io/en/latest/miniconda.html). Then, generate the required environment with:

```bash
conda env create -f environment.yml
```

### Step 2: Program Execution

Activate the environment before running FunCoup:

```bash
conda activate pyfuncoup
```

To update the environment after modifying `environment.yml`:

```bash
conda env update -f environment.yml --prune
```

If you want to export the Conda environment:

```bash
conda env export --name pyfuncoup --no-builds | grep -v "prefix" > environment.yml
```

### Step 3: Database Setup

FunCoup 6 runs on a PostgreSQL database, which must be started within a Docker container before launching the application. The container is defined in the `docker-compose.yml` file.

To install Docker on Ubuntu, follow the instructions [here](https://docs.docker.com/engine/install/ubuntu/).

Once Docker is installed, start the database with:

```bash
sudo docker compose up -d
```

The database stores its files in the `database/` directory. Ensure this directory is not checked into Git as it may grow large over time. To transfer the database between machines, copy the database directory and start the container on the new machine.

To access the database manually using `psql`, install the PostgreSQL client and connect as follows (credentials are stored in `docker-compose.yml`):

```bash
psql -U FunCoup -h localhost FunCoup
```

### Step 4: Running the Application

The main file for running FunCoup is `funCoup.py`. This script provides several functionalities, including generating instances, training networks, constructing networks, generating statistics, and running network analysis.

You can run the following commands from `funCoup.py`:

- To generate an instance by fetching data or using existing data from the database:
  
  ```bash
  python funCoup.py -g --configDir <config_directory>
  ```

- To train and construct a network from an instance:
  
  ```bash
  python funCoup.py -t --configDir <config_directory>
  ```

- To perform network analysis from an instance:
  
  ```bash
  python funCoup.py -a --configDir <config_directory>
  ```

- To generate statistics for the database:
  
  ```bash
  python funCoup.py -s --configDir <config_directory>
  ```

## Cytoscape App

FunCoup also offers a Cytoscape app connected to its API. The app is available as a JAR file at [Bitbucket](https://bitbucket.org/sonnhammergroup/funcoup_cytoscape/) and can be downloaded from the [Cytoscape App Store](https://apps.cytoscape.org/apps/funcoup). The preprint for this application is available on bioRxiv at [https://doi.org/10.1101/2024.10.04.616627](https://doi.org/10.1101/2024.10.04.616627).

## FunCoup 5

FunCoup 5 can be reached at [https://funcoup6.scilifelab.se](https://funcoup6.scilifelab.se).

### Contacts 

* Davide Buzzao (davide.buzzao@scilifelab.se)
* Erik L.L. Sonnhammer (erik.sonnhammer@scilifelab.se)