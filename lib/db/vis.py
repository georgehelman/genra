from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import re
from rdkit import Chem
import StringIO
from rdkit.Chem.Draw import cairoCanvas

import numpy as np
import pymongo

import cairo

def smiles2svg(smiles,size=(100,100)):
    if smiles=='FAIL' or smiles is None:
        return ""
    imageData = StringIO.StringIO()
    surf = cairo.SVGSurface(imageData,size[0],size[1])
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)
    Chem.SanitizeMol(mol,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    ctx = cairo.Context(surf)
    canv = cairoCanvas.Canvas(ctx=ctx, size=size, imageType='svg')
    Chem.Draw.MolToImage(mol, size=size, canvas=canv)
    canv.flush()
    surf.finish()
    return imageData.getvalue()
