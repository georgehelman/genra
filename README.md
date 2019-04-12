# GenRA

This repository contains the code for [https://doi.org/10.1016/j.yrtph.2016.05.008](GenRA) analysis tools.  


# Installation

### Operating system 
This system has only been tested under linux: [Red Hat Enterprise Linux v6/v7](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) and [Ubuntu Server 16.04 LTS](http://www.ubuntu.com). It should be possible to run it under macOS and possibly Windows (if you can install the following requirements). 

### MongoDB
Unless already available, install the latest version of [mongodb](http://www.mongodb.com). MongoDB is a document-oriented database that is quite handy for storing complex data structures [here's](https://docs.mongodb.com/getting-started/shell/introduction/) an introduction to MongoDB. Read the [tutorial on mongodb security](https://docs.mongodb.com/manual/tutorial) and confirm that authentication switched on (edit the /etc/mongod.conf add "authorization:enabled" under security, then restart mongodb).

## Create the mongo database 

Login to mongodb as the administrative user and create a database by the name "organtox_v1". [Create the user](https://docs.mongodb.com/manual/reference/method/db.createUser/) devel (password: devel) that will be used by the code to access the organtox_v1 database. 


Test out the mongodb installation by connecting to the db using the mongo command line client.

# Python

This code has been tested on [Python 2.7](http://python.org). The easiest way to install Python 2.7 is via [Anaconda](https://www.continuum.io/downloads). Install the python packages given in lib/requirements.txt as follows:

```
unix> pip install -r lib/requirements.txt
```

## Jupyter
[Jupyter notebook](http://jupyter.org/) is an interactive computing environment based and is freely available. It can be installed using Anaconda. After jupyter is installed read the [quickstart instructions](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/) to create a configuration file. Set the following variables in your jupyter_notebook_config.py:

```python
c.NotebookApp.port = 7777
```

Enter the notebooks directort and start the notebook server from the command line using (I always run this in a [screen](https://www.gnu.org/software/screen/manual/screen.html)):

```
unix> jupyter notebook
```


# Testing the system

After you have completed the above steps open [this page](http://localhost:7777) in your browser to run different steps of the analysis. The machine learning analysis given in "organ-tox-ml.ipynb", the impact of different factors on F1 performance is given in "organ-tox-stats.ipynb", and the code for reproducing all figures / tables is given in "organ-tox-figs.ipynb".

#Directory structure

notebooks - This folder contains all of the notebooks related to the project. It is loosely organized into subdirectories based on subject matter
	db - pertains to MongoDB ETL operations
	rax - pertains to GenRA analyses
lib - This folder contains all of the user-defined modules used by the notebooks. The modules are loosely organized into subdirectories based on the subject matter
	db - pertains to MongoDB ETL operations
	rax - pertains to GenRA hazard prediction 
	chm - pertains to the creation of chemical fingerprints
	bio - pertains to the creation of biological fingerprints
	utl - includes shell utilities and MySQL queries
scripts - Scripts for updating the MongoDB databases. These can be executed by cronjobs.
services - Related web services

Note: Many notebooks refer to data and fig folders. These are not included in the project as they are too large in size. The notebooks will create these folders at the root directory of the project if they do not already exist.

# Help

Contact shah.imran@epa.gov
