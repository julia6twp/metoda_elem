import numpy as np
from dane import *
from macierz_H import *
from macierz_C import *
from macierz_HBC import *
from wektor_P import *
def punktowanie(g:Grid,gl:GlobalData): #obliczanie współrzędnych punktów
    nd = g.Nodes
    e_tab = g.Elements

    x = np.zeros([int(gl.data["Elements number"]),4]) #tablice dla wsp
    y = np.zeros([int(gl.data["Elements number"]),4])

    for i in range(int(gl.data["Elements number"])):
        x[i] = [nd[int(p)-1].x for p in e_tab[i].id]     #wsp x
        y[i] = [nd[int(p) - 1].y for p in e_tab[i].id]  #wsp y
    return x,y

def MacierzeHiC(x,y,elem,gd): #liczenie macierzy H dla każdego elementu
    Hm = np.zeros([int(gd.data["Elements number"]), 4, 4])
    Cm = np.zeros([int(gd.data["Elements number"]), 4, 4])
    for i in range(int(gd.data["Elements number"])):
        Hm[i] = H(x[i],y[i],elem,gd.data["Conductivity"]).H
        Cm[i] = C(x[i],y[i],elem,gd).C
    return Hm,Cm

def HBCiP(gd,grid,elem): #2,3,0,1 - kolejnosc scian
    HBCw = np.zeros([int(gd.data["Elements number"]), 4, 4])
    Pw = np.zeros([int(gd.data["Elements number"]), 4])
    for i in range(int(gd.data["Elements number"])):
        sciany  = np.array([grid.Elements[i].id[2],grid.Elements[i].id[3],grid.Elements[i].id[0],grid.Elements[i].id[1]])
        det = []
        pom = []
        for j in range(sciany.size):
            j1 = j+1
            if j1 == sciany.size: j1=0
            if np.any(grid.BC == sciany[j]) and np.any(grid.BC == sciany[j1]): pom.append(j)
            det.append( dlugosc_odc([grid.Nodes[int(sciany[j]) - 1],grid.Nodes[int(sciany[j1]) - 1]])/2)
        print(pom)
        hbcclass = HBC(gd,elem,det,grid)
        hbcclass.sciany(np.array(pom))
        hbcclass.oblicz()
        HBCw[i] = hbcclass.HBC
        pclass = P(gd,elem,det,grid)
        pclass.sciany(np.array(pom))
        pclass.oblicz()
        Pw[i] = pclass.wekP
    return np.flip(HBCw,0),np.flip(Pw,0)