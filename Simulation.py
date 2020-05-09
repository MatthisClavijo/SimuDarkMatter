from math import *

def Meca (R) :
    Ma=3.4*10**38
    G=6.67*10**(-11)
    v=sqrt((G*Ma)/R)
    v=v*10**(-3)
    return (v)

print(Meca(1.2*10**(21)))