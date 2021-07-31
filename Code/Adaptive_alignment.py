#!/usr/bin/env python
# coding: utf-8
##########Align all the seismograms using Adaptive alignment stacking##########
##########                New arrivals are stored in t3              ###########

#! /usr/bin/env python

from obspy import read
from obspy.core import UTCDateTime, Stream
from glob import glob
from obspy.taup import TauPyModel
import matplotlib.pyplot as plt
import numpy as np

model = TauPyModel('ak135')

####### readin arguments from gen_aq_YW (contain the path to data) for writing AQ file
arg_fname = 'Path/gen_aq_YW.cmd'
fid = open(arg_fname, 'r')
print 'Arguments read from:', arg_fname
# -- data waveform --
dname = fid.readline().split()[0]
print dname
if dname[-1] == '/': dname = dname[:-1]
origtime = UTCDateTime().strptime(dname.split('/')[-1], '%Y%m%d_%H%M%S')
print origtime
#origtime = UTCDateTime().strptime(dname.split('/')[-1], '%Y%m%d_%H%M%S.a')
#print ' - dname:', dname
#print ' - origin time:', origtime
# -- phase of interest --
channel = fid.readline().split()[0]
print ' - channel:', channel
# -- phase of interest --
akph = fid.readline().split()[0]
print ' - phase of interest:', akph
# -- data window --
offset, wlen = fid.readline().split()[0:2]
offset, wlen = float(offset), float(wlen)
print ' - window offset and length:', offset, wlen
# -- frequency bandpass --
freq_min, freq_max = fid.readline().split()[0:2]
freq_min, freq_max = float(freq_min), float(freq_max)
print ' - frequency bandpass:', freq_min, freq_max
# -- distance range --
dist_min, dist_max = fid.readline().split()[0:2]
dist_min, dist_max = float(dist_min), float(dist_max)
print ' - distance range:', dist_min, dist_max

#fname_list = glob('%s/*%s*snr' % (dname, channel))
#### for event 20110902, with HHT
#fname_list = glob('%s/*%s*snr' % (dname, channel))
fname_list = glob('%s/*HT*snr' % dname)
fname_list.sort()

tshift= list () #theoretical phase arrival time (ak135) for all stations
gcarca = list() #epicentral distance for all staions
ttime = list()
data_stream = Stream()

######################     Readin waveforms from SAC files  ######################
for fname in fname_list:
    try:
        ## computing window and cutting
        tr = read(fname, format='SAC')[0]
        gcarc = tr.stats.sac.gcarc
        evdp = tr.stats.sac.evdp
        if gcarc < dist_min or gcarc > dist_max: continue
        tr.resample(20)  # the sampling rate, the delta is 0.05 here
        arvs = model.get_travel_times(evdp, gcarc, [akph])
        if len(arvs) == 0: continue
        tt = min([ar.time for ar in arvs]) ##tt is the theoretical S arrival time by ak135
        ttime.append(tt)
#         ttto = 0.5 * (min(ttime) + max(ttime))
        arvs0 = model.get_travel_times(evdp, 0.5*(dist_min+dist_max), [akph])
        ttto = min([ar.time for ar in arvs0])
        begtime = origtime + ttto + offset  #Trace begin time
        endtime = begtime + wlen    #Trace end time
        tr.trim(begtime, endtime, pad = 1, fill_value=0.0)
        ts = ttto - tt

        if tr.stats.npts < int(wlen / tr.stats.delta):
            continue #What is this for??
        tshift.append(ts)
        gcarca.append(gcarc)
        data_stream.append(tr)
        tr.stats.sac.t2 = tt
#         tr.write(fname, format='sac')
#         print 'theoretical S arrival time (ak135) is:', tr.stats.sac.t2, 'the moveout correction is:', ts
        print gcarc, tt
    except Exception as ex:
        print fname, ex

## filtering
for tr in data_stream:
    try:
        tr.detrend('linear')
        tr.detrend('demean')
        tr.filter('bandpass', corners=4, zerophase=True, freqmin=freq_min, freqmax=freq_max)
    except Exception as ex:
        print ex
        data_stream.remove(tr)
print 'Number of traces being read:', len(data_stream)

########################   write out waveform in AQ format      ######################
aq_fname = '%s/%04d-%02d-%02d_%s.aq' % (dname, origtime.year, origtime.month, origtime.day, akph)
print "The AQ file name is:", aq_fname
fid = open(aq_fname, 'w')
# -- AQ header
fid.write('%3d\n' % len(data_stream))
fid.write('%8.4f %8.4f %8.2f\n' % (tr.stats.sac.evla, tr.stats.sac.evlo, tr.stats.sac.evdp))
fid.write('%6d %6d %6d\n' % (begtime.year, begtime.month, begtime.day))
fid.write('%6d %5d %7.4f %10.4f\n' % (begtime.hour, begtime.minute, begtime.second+begtime.microsecond / 1e6, ttto+offset))
fid.write('%.4f %3s\n' % (tr.stats.delta, akph))

# - metadata file format
meta_fname = aq_fname.replace('.aq', '.meta')
fid_meta = open(meta_fname, 'w')
# - meta header
fid_meta.write('%s # initial AQ file\n' % aq_fname.replace('.aq', '.asi'))
fid_meta.write('%s # final AQ file\n' % aq_fname.replace('.aq', '.asf'))
fid_meta.write('%.2f  %.2f # event latitude and longitude\n' % (tr.stats.sac.evla, tr.stats.sac.evlo))
fid_meta.write('%.2f  %.2f # window offset from the reference phase and length\n' % (offset, wlen))
fid_meta.write('%.2f  %.2f # frequency band\n' % (freq_min, freq_max))
fid_meta.write('%.2f  %.2f # epicentral distance range\n' % (dist_min, dist_max))
fid_meta.write('%s # phase of interest\n' % akph)
fid_meta.write('%-10s%-10s%-10s%-10s%-10s%-10s\n' % ('stat', 'gcarc', 'evdp', 'lat', 'lon', 'netw'))

# - for each trace
for tr, ts in zip(data_stream, tshift):
    # -- AQ data
    fid.write('%3d %8d %8.4f %6s\n' % (1, len(tr.data), ts, tr.stats.station))
    for d in tr.data:
        fid.write('%12.4e' % d)
    fid.write('\n')
    # -- meta data
    hdr = tr.stats
    fid_meta.write('%-10s%-10.2f%-10.2f%-10.2f%-10.2f%-10s\n' %         (hdr.station, hdr.sac.gcarc, hdr.sac.evdp, hdr.sac.stla, hdr.sac.stlo, hdr.network))

fid.close()
fid_meta.close()

print 'AQ data file:', aq_fname
print 'Metadata file:', meta_fname

##########################  plot seismograms    ######################

fig = plt.figure(figsize=(7.5, 7))
print 'ttto:', ttto
for tr, ts in zip(data_stream, tshift):
    sampin = tr.stats.delta
    time = np.arange(len(tr.data)) * sampin + ttto + offset + ts
    gcarc = tr.stats.sac.gcarc
    tmp = tr.data
    if np.max(np.abs(tmp)) > 0: tmp /= np.max(np.abs(tmp)) #normalisation

    plt.fill_between(time, gcarc, gcarc+0.3*tmp, where=(tmp<0), color='grey', alpha=0.5)
    plt.fill_between(time, gcarc, gcarc+0.3*tmp, where=(tmp>=0), color='blue', alpha=0.5)

plt.savefig('../testfig.pdf')

#plt.show()

###############################  Change directory to dname and and write input tcas.cmd file ######################
import os

os.chdir(dname)
tcmd_fname = './tcas.cmd'

fidcmd = open(tcmd_fname, 'w')
i_num = 5 #number of iterations
pst_index = 3.0 #index for phase stack
coe_perr = 1.15 #coefficient for pick error
Err_min = 25.0 #min error limits(ms)
Err_max = 150.0 #max error limits(ms)
swin_sta = 90.0 #stack window start
swin_len = 40.0 #stack window LENGTH
AQ = '%04d-%02d-%02d_%s.aq' % (origtime.year, origtime.month, origtime.day, akph) #INPUT event file
difft_min = -5.0 #diff time
difft_max = 5.0  #diff time

fidcmd.write('%d number of iterations\n' % i_num)
fidcmd.write('%.1f index for phase stack\n' % pst_index)
fidcmd.write('%.2f coefficient for pick error\n' % coe_perr)
fidcmd.write('%.1f %.1f min, max error limits(ms)\n' % (Err_min, Err_max))
fidcmd.write("%.1f %.1f stack window(start,length)\n" % (swin_sta, swin_len))
fidcmd.write("%s\n" % AQ)
fidcmd.write("%.1f diff time\n" % difft_min)
fidcmd.write("%.1f diff time" % difft_max)
fidcmd.close()

############ Run tcas ############
import os
import subprocess

os.chdir(dname)
# os.getcwd()
cmd = "Path/tcas"
subprocess.call(cmd, shell=True)  # returns the exit code in unix
######### generating asi, asf under dname (data directory)#########

########################  Function: Read from Output AQfiles (asi, asf) for plotting #####################
def read_aq(fname):
    fid = open(fname, 'r')
    ## read aq file header
    line = fid.readline()
    nsta = int(line)
    line = fid.readline()
    evlat, evlon, evdep = np.array(line.split()).astype(float)
    line = fid.readline()
    iyr, imon, iday = np.array(line.split()).astype(int)
    line = fid.readline()
    ihr, imin = np.array(line.split()[0:2]).astype(int)
#    sec  = float(line.split()[2])
    line = fid.readline()
    sampin = float(line.split()[0])

    tshift = list()
    csta = list()
    w = list()
    npts = list()
    data = list()
    for i in range(nsta):
        line = fid.readline()
        w.append(float(line.split()[0])) # weight
        npts.append(int(line.split()[1])) # number of points
        tshift.append(float(line.split()[2])) # time shift
        csta.append(line.split()[3]) # station name

        line = fid.readline()
        tmp = np.array(line.split()).astype(float) * w[i]
        if np.max(np.abs(tmp)) > 0: tmp /= np.max(np.abs(tmp))
        data.append(tmp)
    fid.close()
    return csta, tshift, sampin, data

########################  Plotting Output AQ files (asi, asf) #####################
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import os

print os.getcwd()

if __name__ == '__main__':

    # -- init and final aq files --
    init_fname = '%04d-%02d-%02d.asi' % (origtime.year, origtime.month, origtime.day)
    print ' - initial aq:', init_fname
    final_fname = '%04d-%02d-%02d.asf' % (origtime.year, origtime.month, origtime.day)
    print ' - final aq:', final_fname
#
    STEP = 1

    csta, tshift, sampin, data = read_aq(init_fname)
    csta, tshiftf, sampin, dataf = read_aq(final_fname)

#     print tshift
#     print tshiftf

    ## plot traces
    fig = plt.figure(figsize=(7.5, 7))
    gs1 = gridspec.GridSpec(1, 2)
    gs1.update(left=0.10, right=0.95, bottom=0.08, top=0.70)
    gs2 = gridspec.GridSpec(2, 2)
    gs2.update(left=0.10, right=0.95, bottom=0.75, top=0.95)

    ## Initial - individual traces
    ax1 = plt.subplot(gs1[0, 0])

    print len(gcarca), len(data)

    for n in range(0, len(gcarca), STEP):
        norm_data = data[n] * 0.5
        time = np.arange(len(norm_data)) * sampin + ttto + tshift[n] + offset
#         time = np.arange(len(norm_data)) * sampin + offset
        ax1.fill_between(time, gcarca[n], gcarca[n]+norm_data, where=(norm_data<0), color='grey', alpha=0.5)
        ax1.fill_between(time, gcarca[n], gcarca[n]+norm_data, where=(norm_data>=0), color='blue', alpha=0.5)
    ## Initial - stacked traces
    # Squared stack
    ax1_ss = plt.subplot(gs2[1, 0], sharex=ax1)
    squared_stack = data[-1] * 0.5
    time = np.arange(len(squared_stack)) * sampin + ttto + tshift[-1]  + offset
    ax1_ss.plot(time, squared_stack, lw=1.0, c='k')
    # Linear stack
    ax1_ls = plt.subplot(gs2[0, 0], sharex=ax1)
    linear_stack = data[-2] * 0.5
    time = np.arange(len(linear_stack)) * sampin + ttto + tshift[-2]  + offset
    ax1_ls.plot(time, linear_stack, lw=1.0, c='k')
    ## Labelling
    ax1.set_ylabel('Epi. distance (deg)')
    ax1.set_xlabel('Time (s)')
    ax1.set_title('Initial', fontsize=10)
    ax1.plot([ttto, ttto  ], [dist_min, dist_max], c='b', lw=1.0)
    ax1.annotate(akph, xy=(0, dist_max), color='b', va='bottom', ha='center', fontsize=10)
    ## Final - invidual traces
    ax2 = plt.subplot(gs1[0, 1], sharex=ax1, sharey=ax1)
    for n in range(0, len(gcarca), STEP):
        norm_data = dataf[n] * 0.5
        time = np.arange(len(norm_data)) * sampin + ttto + tshiftf[n]  + offset
#         time = np.arange(len(norm_data)) * sampin - tshiftf[n] + tshift[n] + offset
        ax2.fill_between(time, gcarca[n], gcarca[n]+norm_data, where=(norm_data<0), color='grey', alpha=0.5)
        ax2.fill_between(time, gcarca[n], gcarca[n]+norm_data, where=(norm_data>=0), color='blue', alpha=0.5)
    ## Final - stacked traces
    # Squared stack
    ax2_ss = plt.subplot(gs2[1, 1], sharex=ax1)
    squared_stack = dataf[-1] * 0.5
    time = np.arange(len(squared_stack)) * sampin + ttto + tshiftf[-1]  + offset
    ax2_ss.plot(time, squared_stack, lw=1.0, c='k')
    # Linear stack
    ax2_ls = plt.subplot(gs2[0, 1], sharex=ax1)
    linear_stack = dataf[-2] * 0.5
    time = np.arange(len(linear_stack)) * sampin + ttto + tshiftf[-2]  + offset
    ax2_ls.plot(time, linear_stack, lw=1.0, c='k')
    ## Labelling
    ax2.set_xlabel('Time (s)')
    ax2.set_title('Final', fontsize=10)
    ax2.set_ylim(dist_min, max(gcarca)+1.0)
    ax2.plot([ttto, ttto], [dist_min, dist_max], c='b', lw=1.0)
    ax2.annotate(akph, xy=(0, dist_max), color='b', va='bottom', ha='center', fontsize=10)
    ax2.set_xlim(ttto+offset, ttto+offset+wlen)
#     bname = meta_fname.split('/')[-1].replace('.meta', '')
#     ####plot seismograms
#     ####Initial traces
#     ax1 = plt.subplot(1,2,1)
#     for n in range(0, len(gcarca), 1):
#         norm_data = data[n] * 0.3
#         time = np.arange(len(norm_data)) * sampin - tshift[n] + offset
# #         time = np.arange(len(norm_data)) * sampin + offset
#         ax1.fill_between(time, gcarca[n], gcarca[n]+norm_data, where=(norm_data<0), color='green', alpha=0.5)
#         ax1.fill_between(time, gcarca[n], gcarca[n]+norm_data, where=(norm_data>=0), color='red', alpha=0.5)

#     ####Final traces

#     ax2 = plt.subplot(1,2,2)
#     for n in range(0, len(gcarca), 1):
#         norm_data = data[n] * 0.3
#         time = np.arange(len(norm_data)) * sampin - tshiftf[n] + offset
# #         time = np.arange(len(norm_data)) * sampin - tshiftf[n] + tshift[n] + offset
#         ax2.fill_between(time, gcarca[n], gcarca[n]+norm_data, where=(norm_data<0), color='green', alpha=0.5)
#         ax2.fill_between(time, gcarca[n], gcarca[n]+norm_data, where=(norm_data>=0), color='red', alpha=0.5)

    img_fname = '../%04d-%02d-%02d-g%d_%d.pdf' % (origtime.year, origtime.month, origtime.day, dist_min, dist_max)
    plt.savefig(img_fname)
#    plt.show()


########################  Write new theoretical S arrival times (t3 = ttto - tshiftf) to the header of original files #####################
for fname in fname_list:
    tr = read(fname, format='SAC')[0]
    i = 0
    kstnm = str(tr.stats.sac.kstnm).strip()
    print type(kstnm), kstnm

    for stnm, tti, ttf in zip(csta, tshift, tshiftf):
#         print type(stnm), stnm
        stnm = stnm.strip()
#         print kstnm, id(kstnm), stnm, id(stnm)
        i+=1
        if kstnm == stnm:
            tr.stats.sac.t3 = ttto - ttf
            tr.stats.sac.t4 = ttto + ttf
            print "found the correponding station:", i, kstnm, tr.stats.sac.gcarc
            print kstnm, tr.stats.sac.t2, tr.stats.sac.t3, tr.stats.sac.t4
    tr.write(fname, format='sac')

# csta, tshift, sampin, data = read_aq(init_fname)
# csta, tshiftf, sampin, dataf = read_aq(final_fname)
# print zip(csta, ttto-tshift, ttto-tshiftf)
