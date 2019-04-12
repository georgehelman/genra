from flask import Flask, jsonify, request, json, make_response
import pymongo
import numpy as np
import pandas as pd
import math
import re
from lib import *


con = pymongo.MongoClient("mongodb://devr:devr@pb.epa.gov/genra_dev_v4")
DB = con['genra_dev_v4']
dsstox = DB['compound']

NO_SVG = """<?xml version="1.0" encoding="iso-8859-1"?>
<svg:svg version="1.1" baseProfile="full"
        xmlns:svg="http://www.w3.org/2000/svg"
        xmlns:xlink="http://www.w3.org/1999/xlink"
        xml:space="preserve" width="100px" height="100px" >
<svg:polygon fill="rgb(255,255,255)" stroke="none" stroke-width="0" points="0.00,0.00 100.00,0.00 100.00,100.00 0.00,100.00 0.00,0.00"></svg:polygon>
</svg:svg>"""

app = Flask(__name__)

@app.route('/hello/')
def hello():
    return "hello"

@app.route('/draw')
def draw():
    return app.send_static_file('index.html')

@app.route('/substance/')
def substance():

    sid=request.args.get('sid')

    C = dsstox.find_one({'dsstox_sid':sid})

    if C and 'viz' in C:
        svg = re.sub('>(\S+)</svg:text>',
                     '>&#8226;</svg:text>',
                     C['viz'])
        svg = re.sub('encoding="iso-8859-1"',
                     'encoding="UTF-8"',
                     svg)
        svg = re.sub('font-size="7.00"',
                     'font-size="20.00"',
                     svg)

    else:
        svg = NO_SVG

    response = make_response(svg)
    response.content_type = 'image/svg+xml'
    return response

@app.route('/graph/')
def graph():
    target_sid=request.args.get('target')
    neighbor_sids=request.args.get('neighbors').split(',')
    sims=radial_sims(target_sid,neighbor_sids)
    names=[]
    for sid in neighbor_sids:
        names.append(dsstox.find_one({'dsstox_sid':sid})['name'])

    k0=len(neighbor_sids)
    s0=.1
    W=700
    H=700
    rs=1.0
    img_w=60
    img_h=60
    fp='chm_mrgn'
    rdst='equal'

    r_max=240
    C=[0,0]
    O=[-1*W*.5,-1*H*.5]
    th0=1.32*math.pi
    th_tot=1.9*math.pi

    NN = pd.DataFrame({'dsstox_sid': neighbor_sids, 'jaccard': sims, 'name':names})
    target_df=pd.DataFrame([[target_sid,1,dsstox.find_one({'dsstox_sid':target_sid})['name']]],columns=['dsstox_sid','jaccard','name'])
    NN=pd.concat([target_df,NN],ignore_index=True)
    def shorten(s):
        if len(s)>23:
            s=s[0:21]+'..'
        return s
    NN['name']=NN['name'].apply(shorten)

    NN['d'] = 1 - NN.jaccard
    k0 = NN.shape[0]
    dth = th_tot / k0

    NN.sort_values(by='d', inplace=True)
    # NN.d[NN.d<=0.7]=0.7
    # NN['r'] = r_max*NN.d*rs
    # For now ...
    NN['r'] = r_max
    NN['th'] = th0 + dth * np.arange(0, k0)
    NN['x'] = NN.r * np.cos(NN.th)
    NN['y'] = NN.r * np.sin(NN.th)

    # Add a shorter loc for ending the edge
    NN['xb'] = 0.75 * NN.r * np.cos(NN.th)
    NN['yb'] = 0.75 * NN.r * np.sin(NN.th)

    # Root
    NN.ix[0, 'r'] = 0
    NN.ix[0, 'x'] = C[0]
    NN.ix[0, 'y'] = C[1]
    NN.ix[0, 'th'] = 0

    # Add image coordinates
    NN['v_img_w'] = img_w
    NN['v_img_h'] = img_h
    NN['v_img_x'] = NN.x - img_w * 0.5
    NN['v_img_y'] = NN.y - img_h * 0.5

    # Change to screen coordinates
    NN['v_x'] = NN.x - O[0]
    NN['v_y'] = NN.y - O[1]
    NN['xb'] = NN.xb - O[0]
    NN['yb'] = NN.yb - O[1]

    NN['v_img_x'] = NN.v_img_x - O[0]
    NN['v_img_y'] = NN.v_img_y - O[1]

    Vi = NN.ix[0]

    # Edges
    E2 = NN[['dsstox_sid', 'v_x', 'v_y', 'jaccard', 'xb', 'yb']].ix[1:].copy()
    E2.rename(columns=dict(dsstox_sid='vj', v_x='xj', v_y='yj'), inplace=True)
    E2['vi'] = target_sid
    E2['xi'] = Vi.v_x
    E2['yi'] = Vi.v_y
    E2['lx'] = 0.5 * (E2.xi + E2.xj)
    E2['ly'] = 0.5 * (E2.yi + E2.yj)
    E2['label'] = E2.jaccard.apply(lambda i: '%3.2f%s' % (i, fp[0]))

    return jsonify(dict(edges=E2.to_dict('records'),
                        nodes=NN.to_dict('records'),
                        k0=k0, s0=s0,
                        fp=fp,
                        W=W, H=H,
                        root=target_sid,
                        n=NN.shape[0]))

@app.route('/template')
def template():
    target=request.args.get('dsstox_sid')
    s=request.args.get('s')
    k=request.args.get('k')

    neighbors=searhCollByFP(sid,s0=s,k0=k,)

if __name__ == '__main__':
    app.debug = True
    app.run(host='localhost',port=8010)