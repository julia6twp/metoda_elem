import numpy as np
from elem_uniwersalny import ElementUniwersalny
from Jakobian import *
from dane import *

#macierz C
class C:
    def __init__(self,x,y, elem:ElementUniwersalny,g:GlobalData):
        self.ro = g.data["Density"]
        self.cw = g.data["SpecificHeat"]
        self.elem = elem
        self.jakobian = []
        self.Cmac=np.zeros([self.elem.n**2,4,4])
        self.C = np.zeros([4, 4])


        for i in range(self.elem.n**2):
            self.jakobian.append(Jakobian(x,y,elem,i))

        x = 0
        for i in range(self.elem.n):
            for j in range(self.elem.n):
                N1 = 0.25 * (1 - self.elem.tab[j]) * (1 - self.elem.tab[i])
                N2 = 0.25 * (1 + self.elem.tab[j]) * (1 - self.elem.tab[i])
                N3 = 0.25 * (1 + self.elem.tab[j]) * (1 + self.elem.tab[i])
                N4 = 0.25 * (1 - self.elem.tab[j]) * (1 + self.elem.tab[i])
                b = np.array([N1, N2, N3, N4])[np.newaxis]
                self.Cmac[i] += (self.ro * self.cw * self.jakobian[x].det *self.elem.waga[i]*self.elem.waga[j]*(b*b.T))
                x+=1

            self.C+=self.Cmac[i]







