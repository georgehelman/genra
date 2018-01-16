import matplotlib.text as text
import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import sys
import os 
import pymongo
import importlib

#Set environment variables


# Set up the local source files
TOP = '/'.join(os.getcwd().split('/')[:-1])+'/'

LIB = TOP+'lib'
if not LIB in sys.path: 
    sys.path.insert(0,LIB)
print(LIB)
os.environ['PYTHONPATH']=LIB


DAT_DIR = TOP + '/data/'
FIG_DIR = TOP + '/figs/'

if not os.path.exists(DAT_DIR): os.mkdir(DAT_DIR)
if not os.path.exists(FIG_DIR): os.mkdir(FIG_DIR)


from db.mongo import *
print(importlib.import_module('db.mongo',package='mongo'))

DB = openMongo(db='genra_dev_v4')
print(DB.chm_fp.find_one())
from db.fpsim import *
from rax.genrapred import *
