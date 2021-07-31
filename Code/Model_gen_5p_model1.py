
# coding: utf-8
# Build inputcards of grid search models for running wkbj
# Remember to change inputcard info before running this. The source depth layer is different. Card 5d (two places) and first 5e (The max ray par is set to depth layer), also 5b, 5c if the source is shallow or deeper

#import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
#import os
import itertools


inputcard = '/Tdata/yuweili/GSM/Voon_Model1_210r_0.2'
modeldir = '/Tdata/yuweili/GSM/Voon_model_p5_210_D0.2/'
H = range(2291,2900,20)
D0 = 2291
D3 = 2891
V0 = 7.05

#deifine Vp and rho
Vp = np.linspace(13.1, 13.6, num=len(H)+1, endpoint=True)
rho = np.linspace(5.2, 5.42, num=len(H)+1, endpoint=True)

#print D_len
#print H[0]
i = 0
N = []

for D1 in range(H[0]+20,H[-1],20):
    for D2 in range(D1+20,2891,20):
        for V1 in np.arange(6.8, 7.332, 0.08):
            for V2 in np.arange(V1, 7.5, 0.08):
                for V3 in np.arange(6.9, 7.4, 0.08):

                    i = i + 1
                    D = [D0, D1, D2, D3]
                    V = [V0, V1, V2, V3]
                    N.append(i) #number of models

                    #interpolation: generate models
                    f = interpolate.interp1d(D,V)

                    #Search D1 location and create as D" interface (7th)
                    H = range(2291,2900,20)
                    k = H.index(D1)
                    H.insert(k, D1)
                    Vs = f(H)

                    #Write models to cards (~)
                    #np.savetxt('Modelspl%s.txt' % i, zip(H, Vm), fmt="%i \t %5.4f")

                    fmodel = open(modeldir+"Model%s.txt" % i, "w+") #create and write output card
                    finput = open(inputcard, "r")             #Add info from input card
                    card = open(modeldir+"Model%s.txt" % i, "a+")   #append to output card

                    Minfo = open(modeldir+"modelinfo", "a+") #create and write model information
                    Minfo.write("%d %.1f %.1f %5.4f %5.4f %5.4f\n" % (i, D1, D2, V1, V2, V3))

                    #Add first part info (including Card1 and model above D" (ak135)) from input card (newcard here)
                    for lines in itertools.islice(finput, 0, 56):
                        fmodel.write(lines)
                    fmodel.close()

                    #Add interpolated models from 2291 to 2891km
                    for j in range(len(H)):
                        card.write("%5.3f %5.4f %5.4f %5.4f\n" % (H[j], Vp[j], Vs[j], rho[j]))

                    #Add last line in model (CMB interface) and the rest cards(Card3~Card8)
                    for linee in itertools.islice(finput, 32, None):
                         card.write(linee)

                    card.close()


#                    plt.plot(V, D, 'o', Vs, H, '-')
#                    plt.gca().invert_yaxis()
#                    plt.xlim(6.5, 7.8)
#                    plt.show()

Minfo.close
#print N[0]
print len(N)

