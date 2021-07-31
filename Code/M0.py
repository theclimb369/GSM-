#!/usr/bin/env python
# coding: utf-8

############
#Calculate misfit between data (after stacking) and synthetics for all the models in specific library (after cc time shift)
# Created by YW, May 14, 2019
#############
import os
from glob import glob
import numpy as np
from numpy import linalg as LA
import math
from obspy import read
from scipy.signal import correlate
import time
start_time = time.time()

startnm = 1
nm = 10000
w = 2

dname = '/Tdata/yuweili/Adaptivestacking/Boxtest/Ba60_64_o224_232_ab'
sname = '/Tdata/yuweili/GSM/modellib_460_D0.2/'
refernm = dname.split('/')[-1]

os.chdir(sname)

datalist = glob('%s/ITD_St*respl' % dname)

datalib = dict()     # Create dict excluding gcarc missing from data
gcarc_r = np.arange(70, 85.5, 0.5)  # Initial list of gcarc for creating keys
gnm = len(gcarc_r)   # Length of gcarc list
sampin = 5

# Read data and create dictionary to store gcarc, data
for data in datalist:
    tr = read(data)[0]
    trd = tr.data[0:500]
    if np.max(np.abs(trd)) > 0: trd /= np.max(np.abs(trd))
    gcarc = tr.stats.sac.gcarc
    if gcarc in gcarc_r:
        try:
            t6 = tr.stats.sac.t6
            t7 = tr.stats.sac.t7
        except:
            datalib[gcarc] = [trd, 0, 0]
        else:
            sstart = np.floor(t6) * sampin + 1
            send = np.ceil(t7) * sampin + 1
            datalib[gcarc] = [trd, int(sstart), int(send)]

       # print gcarc, datalib[gcarc][1]

# Loop through all models to CC syn and data and shift syn
# Then calculate misfits

misfit = np.empty((gnm,3))  #####3D array to store gcarc, 2 types of misfits

for i in range(startnm,nm+startnm):
    os.chdir('model'+str(i))
    synlist = glob('seis*.cut')
    for j, syn in enumerate(synlist):
        sr = read(syn)[0]
        srd = sr.data[0:500]
        if np.max(np.abs(srd)) > 0: srd /= np.max(np.abs(srd))
        gcarc = round(sr.stats.sac.gcarc,1)
        misfit[j, 0] = gcarc
        if gcarc in datalib.keys():
            data = datalib[gcarc][0]
            sst = datalib[gcarc][1]
            sed = datalib[gcarc][2]
            lag = np.argmax(correlate(data, srd))
            c_sig = np.roll(srd, shift=int(np.ceil(lag)))
            ## part 1
            if sst:
                seg1 = np.sum((c_sig[0:sst]-data[0:sst]) ** 2)
                seg2 = np.sum((c_sig[sst:sed]-data[sst:sed]) ** 2)
                seg3 = np.sum((c_sig[sed:]-data[sed:]) ** 2)

                misfit[j, 1] = math.sqrt((seg1 + seg2 * w + seg3) / (len(data)+(w-1)*(sed-sst)))
                misfit[j, 2] = misfit[j, 1] / (LA.norm(data) / math.sqrt(len(data)))
#                misfit[j, 2] = (LA.norm(c_sig[0:sst]-data[0:sst]) + w * LA.norm(c_sig[sst:sed]-data[sst:sed]) + LA.norm(c_sig[sed:]-data[sed:])) / (len(data)+(w-1)*(sed-sst))
            else:
                misfit[j, 1] = LA.norm(c_sig-data) / math.sqrt(len(data))
                misfit[j, 2] = LA.norm(c_sig-data)/LA.norm(data)
        else:
            misfit[j, 1] = np.nan
            misfit[j, 2] = np.nan
    misfit =  misfit[np.argsort(misfit[:,0])]
#    average = np.nanmean(misfit[:,1]), np.nanmean(misfit[:,2])
    np.savetxt(refernm+'misfitw2_'+str(i), misfit, fmt=['%0.1f', '%0.4f', '%0.4f'])
    os.chdir(sname)
    print ("The %d th model done") % i
    print ("--- %s seconds ---" % (time.time() - start_time))
