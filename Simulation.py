from math import *
import numpy as np
import matplotlib.pyplot as plt
def Meca (R) :
    Ma=3.4*10**38
    G=6.67*10**(-11)
    v=sqrt((G*Ma)/R)
    v=v*10**(-3)
    return (v)

#on prend les valeurs absolues des différentes vitesses radiales dans le tableau suivant

#Voici un tableau avec les vitesses mesurées par effet doppler (les vitesses sont en km/s)

name = np.array(['Nu And', 'Mayall II','Zet And' ])
vit = np.array([10.30, 332.0, 24.43])







