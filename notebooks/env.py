%load_ext autoreload
%autoreload 2
%pylab inline
%matplotlib inline

import matplotlib.text as text
import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import sys
import rpy2 
import os 

#Set environment variables


# Set up the local source files
#TOP = '/'.join(os.getcwd().split('/')[:-2])+'/'
TOP = "/share/home/ishah/ipynb/chiron/genra-analysis/"

LIB = TOP+'lib'
if not LIB in sys.path: 
    sys.path.insert(0,LIB)

os.environ['PYTHONPATH']=LIB


DAT_DIR = TOP + '/data/'
FIG_DIR = TOP + '/figs/'

if not os.path.exists(DAT_DIR): os.mkdir(DAT_DIR)
if not os.path.exists(FIG_DIR): os.mkdir(FIG_DIR)


from db.mongo import *

DB = openMongo(db='genra_dev_v4')