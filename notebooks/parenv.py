import matplotlib.text as text
import pandas as pd
import numpy as np
import pylab as pl
import scipy as sp
import sys
import os 
import pymongo
import importlib
from rdkit import Chem as chm
from rdkit.Chem import Lipinski as lip
from __future__ import division

#Set environment variables
print('hello')

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
from db.fpsim import *
from rax.genrapred import *

def get_phys_fp(compound):
    c=[]
    c.append(compound['mol_weight']/500)
    logp=get_logp(compound['dsstox_sid'])
    logp= logp/10 if logp else logp
    c.append(logp)
    m=chm.MolFromSmiles(compound['smiles'])
    c.append(lip.NumHDonors(m)/5)
    c.append(lip.NumHAcceptors(m)/10)
    return c
