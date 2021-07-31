### Find the minimum log likelihood for each value to certain parameter, e.g. for D1, when D1=D10, search through othe parameters andfind the minimum L.
### Calculate ration = Lc - Lco (-2log(L/Lo))
### Input: Likelihood files (Model# in order, log likelihood, likelihood)
# Output: Array of paremeter and mimimum Lc(-2log(L)) , log Likelihood ratio Lc- Lco (-2log(L/Lo)) e.g. D10, Lc, Lc-Lco
## Created by YW on April 13, 2020
## Last modified on April 15, 2020

import numpy as np
#import time
#start_time=time.time()

###### new  Function ######
def likelihood_ratio_1D_cut(para, para_list, msfl):
    """
    Calculate 1D LR for para
    e.g. D1_L = likelihood_ratio_1D(D1, D1_list, msfl)
    """
    lglike = msfl
    Lo = min(lglike)  #Likelihood for optimal model

    likelihood = np.empty((len(para_list),3))
    for j, pari in enumerate(para_list):
        mask = (para == pari)
        Li = min(lglike[mask])
        #Ri = Lo / Li
        Ri = Li - Lo
        likelihood[j,:] = pari, Li, Ri
    return likelihood

def likelihood_ratio_2D(para1, para2, para_list1, para_list2, msfl):
    """
    Calculate 2D LR for para1, para2
    e.g. D1V1_L = likelihood_ratio_2D(D1, V1, D1_list, V1_list)
    """
    lglike = msfl
    Lo = min(lglike)  #Likelihood for optimal model

    likelihood = np.empty((len(para_list1)*len(para_list2), 4))
    j = 0
    for par1 in para_list1:
        for par2 in para_list2:
            mask = (para1 == par1) * (para2 == par2)
            if len(lglike[mask]) != 0:
                Li = min(lglike[mask])
                Ri = Lo/Li
                #Ri = Li - Lo
                likelihood[j:] = par1, par2, Li, Ri
                j += 1
    return likelihood


###### Old  Function ######
def likelihood_ratio_1D(para, para_list, msfl):
    """
    Calculate 1D LR for para
    e.g. D1_L = likelihood_ratio_1D(D1, D1_list, msfl)
    """
    lglike = np.loadtxt(msfl, usecols=[1])

    Lo = min(lglike)  #Likelihood for optimal model

    likelihood = np.empty((len(para_list),3))
    for j, pari in enumerate(para_list):
        mask = (para == pari)
        Li = min(lglike[mask])
        #Ri = Lo / Li
        Ri = Li - Lo
        likelihood[j,:] = pari, Li, Ri
    return likelihood

###### Main Section ######
dir = ""
os.chdir(dir)
np.set_printoptions(suppress=True)

mdfl = "modelinfoall"  #model information file

# Load para list from model info file
D1 = np.loadtxt(mdfl, usecols=[1])
D2 = np.loadtxt(mdfl, usecols=[2])
V1 = np.loadtxt(mdfl, usecols=[3])
V2 = np.loadtxt(mdfl, usecols=[4])
V3 = np.loadtxt(mdfl, usecols=[5])
VI = V2-V1
#VI = (V2-V1)*100/V1
G = (V2-V1)/(D2-D1)
Th = D2-D1

# Load Likelihood from likelihood file
lglike = np.loadtxt(msfl, usecols=[1])
Lo = min(lglike)  #Likelihood for optimal model

assert(len(D1) == len(lglike))
# parameter ranges
V1_list = np.unique(V1)
V2_list = np.unique(V2)
V3_list = np.unique(V3)
D1_list = np.unique(D1)
D2_list = np.unique(D2)
VI_list = np.unique(VI)
G_list = np.unique(G)
Th_list = np.unique(Th)


# Calculate 1D L and ratios for each value in each parameter
#D1_L = likelihood_ratio_1D(D1, D1_list)
#print D1_L
#D2_L = likelihood_ratio_1D(D2, D2_list)
#V1_L = likelihood_ratio_1D(V1, V1_list)
#V2_L = likelihood_ratio_1D(V2, V2_list)
#V3_L = likelihood_ratio_1D(V3, V3_list)
VI_L = likelihood_ratio_1D(VI, VI_list)
#print VI_L
G_L = likelihood_ratio_1D(G, G_list)
Th_L = likelihood_ratio_1D(Th, Th_list)
#print G_L
"""
"""
# Calculate 2D L and ratios
D2V2_L = likelihood_ratio_2D(D2, V2, D2_list, V2_list)
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
fig = plt.figure()
ax = plt.axes(projection='3d')
points = ax.scatter(D2V2_L[:,0], D2V2_L[:,1], D2V2_L[:,3], c=D2V2_L[:,3], cmap='viridis', linewidth=0.5)
X, Y = np.meshgrid(D2_list, V2_list)
z = 3.84*np.ones(len(np.ravel(X)))
Z = z.reshape(X.shape)
ax.plot_surface(X, Y, Z, cmap='binary')
plt.colorbar(points)
plt.show()
#print D2V1_L
#D1V1_L = likelihood_ratio_2D(D1, V1, D1_list, V1_list)
#D1V2_L = likelihood_ratio_2D(D1, V2, D1_list, V2_list)
#D1V3_L = likelihood_ratio_2D(D1, V3, D1_list, V3_list)
#D1D2_L = likelihood_ratio_2D(D1, D2, D1_list, D2_list)
#V1D2_L = likelihood_ratio_2D(V1, D2, V1_list, D2_list)
#V1V2_L = likelihood_ratio_2D(V1, V2, V1_list, V2_list)
#V1V3_L = likelihood_ratio_2D(V1, V3, V1_list, V3_list)
#V2D2_L = likelihood_ratio_2D(V2, D2, V2_list, D2_list)
#V2V3_L = likelihood_ratio_2D(V2, V3, V2_list, V3_list)
#V3D2_L = likelihood_ratio_2D(V3, D2, V3_list, D2_list)
#
#print ('===%s seconds===' % (time.time() - start_time))



#Plot notes (2D and 1D):
#2D Plot
plt.figure(figsize=(20,10))
st = plt.suptitle("2D Likelihood", fontsize="x-large")
st.set_y(0.95)

plt.subplot(4,4,1)
points = plt.scatter(D1V1_L[:,0], D1V1_L[:,1], c=D1V1_L[:,3], s=50)
plt.colorbar(points)
plt.show()

#1D Plot
import matplotlib.pyplot as plt
fig, axes = plt.subplots()
fig.subplots_adjust(hspace=0.6,wspace=0.3)
ax1 = plt.subplot(311)
ax1.scatter(VI_L[:,0], VI_L[:,2], marker='o', s=10)
ax1.set(xlabel='VI (V2-V1)', ylabel='-2log(Lo/Lmin)')

ax2 = plt.subplot(312)
ax2.scatter(G_L[:,0], G_L[:,2], marker='o', s=10)
ax2.set(xlabel='G (s^-1)', ylabel='-2log(Lo/Lmin)')

ax3 = plt.subplot(313)
ax3.scatter(Th_L[:,0], Th_L[:,2], marker='o', s=10)
ax3.set(xlabel='Th (D2-D1) (km)', ylabel='-2log(Lo/Lmin)')

plt.savefig("Para_Ll1D.pdf")
plt.show()


fig, axes = plt.subplots()
fig.subplots_adjust(hspace=0.6,wspace=0.3)
ax1 = plt.subplot(321)
ax1.plot(D1_L[:,0], D1_L[:,2], marker='o', markersize=3)
ax1.set(xlabel='D1 (km)', ylabel='Lo/Lmin')
ax2 = plt.subplot(323, sharex = ax1, sharey = ax1)
ax2.plot(D2_L[:,0], D2_L[:,2], marker='o', markersize=3)
ax2.set(xlabel='D2 (km)', ylabel='Lo/Lmin')
ax3 = plt.subplot(322, sharey = ax1)
ax3.plot(V1_L[:,0], V1_L[:,2], marker='o', markersize=3)
ax3.set(xlabel='V1 (km)', ylabel='Lo/Lmin')
ax4 = plt.subplot(324, sharex = ax3, sharey = ax1)
ax4.plot(V2_L[:,0], V2_L[:,2], marker='o', markersize=3)
ax4.set(xlabel='V2 (km)', ylabel='Lo/Lmin')
ax5 = plt.subplot(326, sharex = ax3, sharey = ax1)
ax5.plot(V3_L[:,0], V3_L[:,2], marker='o', markersize=3)
ax5.set(xlabel='V3 (km)', ylabel='Lo/Lmin')

