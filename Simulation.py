from math import *
import numpy as np
import scipy.special as sc
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def Meca (R,V) : #R en kPC et V en km/s
    G=6.67*10**(-11)
    #L = getL(16,1,2.2,R)
    MTotDyna = ((kPCToKm(rkPC[-1]) * (vit[-1]**2) * 10**9)/G)
    MDyna = ((kPCToKm(R) * (V**2) * 10**9)/G)
    rapport = MDyna/MTotDyna
    MTotBar = 3.6727e41
    MBar = rapport * MTotBar

    v=sqrt((G*MBar)/(kPCToKm(R)*10**3))
    v=v*10**(-3)
    calR.append(R)
    calVit.append(v)

def kPCToKm(kPCToKm):
    km = kPCToKm * 3.086e16
    return(km)


#on prend les valeurs absolues des différentes vitesses radiales dans le tableau suivant

#Voici un tableau avec les vitesses mesurées par effet doppler (les vitesses sont en km/s)

def datasetRead(datasetPATH):
    DataFileIn = open(datasetPATH, "r")
    DataList = DataFileIn.readlines()
    DataList.sort()
    DataFileIn.close()
    return DataList

def datasetSplit(DataList):
    for i in range(len(DataList)):
        traitement = DataList[i].split(" ")
        for i2 in DataList[i].split(" "):
            if(i2 == ''):
                traitement.remove(i2)
        if(len(traitement) < 2):
            break
        rkPC.append(float(traitement[0]))
        vit.append(float(traitement[1].split("\n")[0]))

    rkPC.pop(0)
    vit.pop(0)

    
def getL(Amp_eff,R_eff,n,R): # R in kPC
    s1 = Sersic1D(Amp_eff,R_eff,n)
    I = s1(R)
    Ie = 0
    Re = 0
    bn = ( (log(I) + log(Ie))/( ((R/Re)**(1/n)) - 1) )

    x = bn *((R/Re)**(1/n))
    gammaInc = sc.gammainc(2*n,x)*sc.gamma(2*n)
    L = Ie*(Re**2)*2*pi*n*( (exp(bn))/(bn**(2*n)) )*gammaInc

    return(L)

def showPlot():
    plt.title('Comparaison de l\'observation des vitesses des étoiles et des vitesses simulées')
    blue_patch = mpatches.Patch(color='blue', label='Simulation')
    red_patch = mpatches.Patch(color='red', label='Observation')
    plt.legend(handles=[blue_patch,red_patch])
    plt.xlabel('kpc')
    plt.ylabel('km/s')

    plt.plot(calR,calVit,'b.',markersize=1)
    plt.plot(rkPC,vit,'r.',markersize=1)

    plt.show()


L = 0
rkPC = []
vit = []
calR = []
calVit = []

data = datasetRead("/home/toor/PythonSim/SimuDarkMatter/0224.dat")
datasetSplit(data)

for i in range(len(rkPC)):
    Meca(rkPC[i],vit[i])

showPlot()