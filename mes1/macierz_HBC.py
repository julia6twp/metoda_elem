import numpy as np
from macierz_H import *
from dane import *

#macierz HBC

class HBC:
    def __init__(self,gd:GlobalData,elem:ElementUniwersalny,detJ,g:Grid):
        self.alfa = gd.data["Alfa"]
        self.elem = elem
        self.detJ = detJ
        self.s = np.array([],dtype=int);
        self.bc=g.BC
        self.hbc = np.zeros([4,4,4]) #sciany
        self.HBC = np.zeros([4,4])  #suma scian
        self.pom = np.zeros([elem.n, 2])
        self.p = np.zeros([4, elem.n, 2])


    def sciany(self, n:np.array(int),): #na ktore sciany mam nakladac warunek brzegowy
        self.s = n
        for i in range(self.elem.n):
            self.pom[i] = [-1, self.elem.tab[i]]
        self.p[3] = self.pom
        for i in range(self.elem.n):
            self.pom[i] = [self.elem.tab[i], 1]
        self.p[2] = self.pom
        for i in range(self.elem.n):
            self.pom[i] = [1, self.elem.tab[i]]
        self.p[1] = self.pom
        for i in range(self.elem.n):
            self.pom[i] = [self.elem.tab[i], -1]
        self.p[0] = self.pom


    def oblicz(self):
        for i in self.s:
            for j in range(self.elem.n):
                N1 = 0.25 * (1 - self.p[i][j][0]) * (1 - self.p[i][j][1])
                N2 = 0.25 * (1 + self.p[i][j][0]) * (1 - self.p[i][j][1])
                N3 = 0.25 * (1 + self.p[i][j][0]) * (1 + self.p[i][j][1])
                N4 = 0.25 * (1 - self.p[i][j][0]) * (1 + self.p[i][j][1])
                pom = np.array([N1,N2,N3,N4])[np.newaxis]
                self.hbc[i] +=self.elem.waga[j]*(pom.T*pom)*self.alfa #hbc dla kazdego punktu
            self.hbc[i] *=self.detJ[int(i)] #mnożenie przez jakobian
            self.HBC +=self.hbc[i] #sumowanie dla każdej sciany dla jednego elementu

def dlugosc_odc(n:[]):
    return np.sqrt((n[0].x-n[1].x)**2+(n[0].y-n[1].y)**2)