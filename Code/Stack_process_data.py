#!/usr/bin/env python
# coding: utf-8
################################################
# Linear stacking every 0.5 degrees and write into St*date.sac files
# Created by YW on May 12, 2019
################################################



### Read deconvolution files *.sac ###
import numpy as np
from obspy import read
from glob import glob
from obspy.io.sac import SACTrace
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.signal import correlate
from obspy.signal.cross_correlation import xcorr_max


dname = '/Tdata/yuweili/Adaptivestacking/CentralAmerica/20080903_112514'
rdep = 560
#dname = '/Tdata/yuweili/GSM/20110902refer'
os.chdir(dname)
#Decon_list = glob('/Tdata/yuweili/Adaptivestacking/20110902_134709/2011245.*.sac') + glob('/Tdata/yuweili/Adaptivestacking/20110620_163601/2011171.*.mul')
Decon_list = glob('%s/20?????.*.sac' % dname)
#Decon_list = glob('%s/ITD.*.sac' % dname)
Decon_list.sort()
print len(Decon_list)
evedate = Decon_list[0].split('/')[-1].split('.')[0]
print "event date: ", evedate

##### stack from 70 to 85 degrees every 0.5 degrees
gcarc_min = 69.75
gcarc_max = 85.5
gcarc_r = np.arange(gcarc_min, gcarc_max, 0.5)


for n in range(0, len(gcarc_r)-1):
    stack = 0
    snm = 0   # stacking traces for each bin for caluating average
    #.tolist()
    print "In this distance range:", gcarc_r[n], gcarc_r[n+1]

    for Deconf in Decon_list:
        decon = read(Deconf, format='SAC')[0]
#        decon.detrend('linear')
#        decon.detrend('demean')
#        decon.filter('bandpass', corners=4, zerophase=True, freqmin=0.02, freqmax=0.2)
        dep = decon.stats.sac.evdp
        # Depth correction
        gcarc = 0.002 * (dep - rdep) + decon.stats.sac.gcarc
        t3 = decon.stats.sac.t3
#     for index, files in dictl:
        if gcarc > gcarc_r[n] and gcarc < gcarc_r[n+1]:
            #begtime = decon.stats.starttime - decon.stats.sac.b + t3 - 20
            #endtime = begtime + 120
            #decon.trim(begtime, endtime, pad = 1, fill_value=0.0)
            #print begtime, endtime
            if type(stack) == int:
                stack = stack + decon.data
                print stack
                print "add file:", Deconf.split('/')[-1]
            if type(stack) != int:
                shift, value = xcorr_max(correlate(decon.data, stack))
                stack = stack + decon.data
                #lag = np.argmax(correlate(stack, decon.data))
                #c_sig = np.roll(decon.data, shift=int(np.ceil(lag)))
                #stack = stack + c_sig
                #print "add file:", Deconf.split('/')[-1], "corr_coe: ", value
            snm = snm + 1
            print "There are %d traces" % snm

    if type(stack) != int:   ###Avoid no data in certain distance range
        stack = stack / float(snm)
        stack_gcarc = 0.5 * (gcarc_r[n] + gcarc_r[n+1])
    # #     print "The stacking result is:", stack

        #### Write the stacking results in SAC format
        header = {'evla': decon.stats.sac.evla, 'evlo': decon.stats.sac.evlo, 'evdp': decon.stats.sac.evdp, 'nzyear': decon.stats.sac.nzyear,
          'nzjday': decon.stats.sac.nzjday, 'nzhour': decon.stats.sac.nzhour, 'nzmin': decon.stats.sac.nzmin, 'nzsec': decon.stats.sac.nzsec,
          'nzmsec': decon.stats.sac.nzmsec, 'delta': decon.stats.sac.delta, 'gcarc': stack_gcarc, 'dist': decon.stats.sac.dist, 'az': 30.0, 'kcmpnm': decon.stats.sac.kcmpnm, 'knetwk': decon.stats.sac.knetwk, 'o': decon.stats.sac.o}

        ## write data files after deconvolution
        fname = 'ITD_St_%s_%.2f.sac' % (evedate, stack_gcarc)
        sac = SACTrace(data = np.asarray(stack), **header)
        sac.write(fname)
        print "Writing stacking file:", fname




###### Seismograms plot #######################
plt.figure(figsize=(10, 5))
gs1 = gridspec.GridSpec(1, 2) ###gs1: two seismograms
gs1.update(left=0.10, right=0.95, bottom=0.08, top=0.70)
gs2 = gridspec.GridSpec(1, 2) ###gs2: linear and squared stackings before and afterwards
gs2.update(left=0.10, right=0.95, bottom=0.75, top=0.95)

Stack_list = glob('%s/ITD_St*.sac' % dname)

### Read deconvolution results
dictl = list()
for Deconf in Decon_list:
    decon = read(Deconf, format='SAC')[0]
#     decon.detrend('linear')
#     decon.detrend('demean')
#     decon.filter('bandpass', corners=4, zerophase=True, freqmin=0.02, freqmax=0.2)
    gcarc = decon.stats.sac.gcarc
    dicl = gcarc, decon.data
    dictl.append(dicl)
delta =  decon.stats.sac.delta
### Read stacking results
stackl = list()
for Stackf in Stack_list:
    stackf = read(Stackf, format='SAC')[0]
#     stackf.detrend('linear')
#     stackf.detrend('demean')
#     stackf.filter('bandpass', corners=4, zerophase=True, freqmin=0.033, freqmax=0.2)
    gcarc = stackf.stats.sac.gcarc
    stackinfo = gcarc, stackf.data
    stackl.append(stackinfo)


### plot raw data
ax1 = plt.subplot(gs1[0, 0])
for tr in dictl:
    data_norm = tr[1] * 0.25
    if np.max(np.abs(data_norm)) > 0: data_norm /= np.max(np.abs(data_norm))
    time1 = np.arange(len(data_norm)) * delta
#     ax1.fill_between(time1, tr[0], tr[0]+data_norm, where=(data_norm<0), color='grey', alpha=0.5)
#     ax1.fill_between(time1, tr[0], tr[0]+data_norm, where=(data_norm>=0), color='blue', alpha=0.5)
    plt.plot(time1, tr[0]+data_norm, c="k", linewidth=0.5)
# plt.ylim(dis_max+1,dis_min-1)
# plt.ylim(90, 80)



### plot deconvolved results
ax2 = plt.subplot(gs1[0, 1], sharex=ax1, sharey=ax1)
for dr in stackl:
    tmp = dr[1] * 0.25
    if np.max(np.abs(tmp)) > 0: tmp /= np.max(np.abs(tmp))
    time = np.arange(len(dr[1])) * delta
#     ax2.fill_between(time, dr[0], dr[0]+tmp, where=(tmp<0), color='grey', alpha=0.5)
#     ax2.fill_between(time, dr[0], dr[0]+tmp, where=(tmp>=0), color='blue', alpha=0.5)
    plt.plot(time, dr[0]+tmp, c="r", linewidth=0.5, ls='--')

plt.xlim(0, 120)
# plt.ylim(dis_max+1, dis_min-1)
plt.ylim(gcarc_min, gcarc_max)

plt.savefig('Stacking_%.1f_%.1f.pdf' % (gcarc_min, gcarc_max))
#plt.show()




## Sheng, removing the max and min value of each time point within each bin
## Results not as good as linear stacking, some Scd features got diminished
#      stack = np.linspace(0, 0, len(dictl[0][1]) )
    #lst_idx_to_stack = [idxfile for (idxfile,files) in enumerate(dictl) if files[0] > gcarc_r[n] and files[0] < gcarc_r[n+1] ]
    #print(type(stack), stack.size )
    #for idx_t, junk in enumerate(stack):
    #    tmp = [ dictl[idx_to_stack][1][idx_t] for idx_to_stack in lst_idx_to_stack]
    #    if len(tmp) > 2:
    #        value = (np.sum(tmp) - np.min(tmp) - np.max(tmp) ) / (len(tmp) - 2)
    #    elif len(tmp) > 0:
    #        value = np.sum(tmp) / len(tmp)
    #    else:
    #        value = 0.0
    #    stack[idx_t] = value
    #print "add file:", files[0]


# In[ ]:




