from math import *
import numpy as np
import scipy.special as sc
#from scipy.special import gammainc
from astropy.modeling.models import Sersic1D
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def Sim (Amp_eff,R_eff,n,R,V) : #R en kPC, V en km/s, Amp_eff R_eff et n du modèle sersic étudié
    #s1 = Sersic1D(Amp_eff,R_eff,n)
    R = R*1e3 #R en pc
    G=6.67*10**(-11)
    #a = 40
    #b = a*sqrt(1-(0.56**2))
    #Lsol = 3.839e26
    Msol = 1.98892e30
    Lbulge = Lr[int((R/1e3)/0.001)-2]
    #Lbulge = getLbulge(Amp_eff,R_eff,n,R)
    Lbulgelist.append(Lbulge)
    #Ldisk = s1(R)*2*pi*(b/a)*(R**2)
    #Ldisklist.append(Ldisk)

    #if(Lbulge > Ldisk):
        #Lmax.append(Lbulge)
    #else:
        #Lmax.append(Ldisk)

    #Mcentre = 3.3812e38
    Mbulge = 5*Lbulge*Msol
    #Mdisk = 5*Ldisk*Msol
    
    MBar = Mbulge*10
    MBarlist.append(MBar)
    Mdyna = (((kPCToKm(R))*(V*1e3)**2)/G)
    #print(Mdyna)
    Mdynalist.append(Mdyna)

    v=sqrt((G*MBar)/(kPCToKm(R))) #kPCToKm(R) donne des m puisque R est en pc
    v=v*10**(-3)
    calR.append(R/(1e3))
    calVit.append(v)

    #print(G*MBar)
    #print((kPCToKm(R)))

def kPCToKm(kPCToKm):
    km = kPCToKm * 3.086e16
    return(km)


#on prend les valeurs absolues des différentes vitesses radiales dans le tableau suivant

#Voici un tableau avec les vitesses mesurées par effet doppler (les vitesses sont en km/s)

def datasetRead(datasetPATH): #Le chemin relatif ne marche pas, METTRE DONC CHEMIN ABSOLUE
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

    
def getLbulge(Amp_eff,R_eff,n,R): # R en pc
    s2 = Sersic1D(Amp_eff,R_eff,n)
    a = 1
    b = a*sqrt(1-(0.56**2))
    #I = s2(R)
    Re = R_eff
    Ie = s2(Re)
    bn = 1.9992*n - 0.3271

    #x = bn *((R/Re)**(1/n))
    #gammaInc = gammainc(2*n,x)*sc.gamma(2*n)

    L = (2*pi*(b/a)*Ie*(R**2)*exp(bn)*n*sc.gamma(2*n))/(bn**(2*n))
    print(L)

    #L = Ie*(Re**2)*2*pi*n*( (exp(bn))/(bn**(2*n)) )*gammaInc

    return(L)

def showPlot():
    plt.title('Comparaison de l\'observation des vitesses des étoiles et des vitesses simulées')
    blue_patch = mpatches.Patch(color='blue', label='Simulation')
    red_patch = mpatches.Patch(color='red', label='Observation')
    plt.legend(handles=[blue_patch,red_patch])
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('kpc')
    plt.ylabel('km/s')

    plt.plot(calR,calVit,'b.',markersize=1)
    plt.plot(rkPC,vit,'r.',markersize=1)


    #blue_patch = mpatches.Patch(color='blue', label='Lbulge')
    #green_patch = mpatches.Patch(color='green', label='Ldisk')
    #red_patch = mpatches.Patch(color='red', label='Lmax')
    #plt.legend(handles=[blue_patch,green_patch,red_patch])

    #plt.plot(rkPC,Lbulgelist,'b.',markersize=1)
    #plt.plot(rkPC,Ldisklist,'g.',markersize=1)
    #plt.plot(rkPC,Lmax,'r.',markersize=1)




    #blue_patch = mpatches.Patch(color='blue', label='MBar')
    #red_patch = mpatches.Patch(color='red', label='Mdyna')
    #plt.legend(handles=[blue_patch,red_patch])

    #plt.plot(rkPC,MBarlist,'b.',markersize=1)
    #plt.plot(rkPC,Mdynalist,'r.',markersize=1)

    plt.show()


def iterativeLuminosity(Amp_eff,R_eff,n):
    s2 = Sersic1D(Amp_eff,R_eff,n)
    a = 1
    b = a*sqrt(1-(0.56**2))
    #Re = R_eff
    #Ie = s2(Re)
    #bn = 1.9992*n - 0.3271

    r = np.arange(0,(max(rkPC)),0.001)
    for i in r:
        Ilist.append(s2(i))
    for index in range(len(Ilist)-1):
        r1 = r[index]
        r2 = r[index+1]
        #I1 = Ilist[index]
        I2 = Ilist[index+1]
        L2 = 2*pi*(b/a)*I2*(((r2*1e3)**2)-((r1*1e3)**2))
        if (index == 0):
            Lr.append(L2)
        else:
            Lr.append(L2 + Lr[index-1])

L = 0
Ilist = []
Lr = []
MBarlist = []
Mdynalist = []
Lbulgelist = []
Ldisklist = []
Lmax = []
rkPC = []
vit = []
calR = []
calVit = []

data = datasetRead("/home/toor/PythonSim/SimuDarkMatter/0224.dat")
datasetSplit(data)

iterativeLuminosity(16,1,2.2)

for i in range(len(rkPC)):
    Sim(15.77,0.82,2.18,rkPC[i],vit[i])

showPlot()