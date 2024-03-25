import numpy as np
from elem_uniwersalny import ElementUniwersalny
from Jakobian import *

#macierz H
class H:
    def __init__(self,x,y, elem:ElementUniwersalny,k):
        self.elem = elem
        self.k = k
        self.jakobian = []
        for i in range(self.elem.n**2):
            self.jakobian.append(Jakobian(x,y,elem,i))
        self.Nx = np.zeros([self.elem.n ** 2, 4])
        self.Ny = np.zeros([self.elem.n ** 2, 4])
        self.Hmac=np.zeros([self.elem.n**2,4,4])
        self.H = np.zeros([4, 4])
        for i in range(self.elem.n**2):
            for j in range(4):
                self.Nx[i][j]=self.jakobian[i].inv[0][0] * self.elem.nKsi[i][j] + self.jakobian[i].inv[0][1] * elem.nEta[i][j]
                self.Ny[i][j]=self.jakobian[i].inv[1][0] * self.elem.nKsi[i][j] + self.jakobian[i].inv[1][1] * self.elem.nEta[i][j]
        for i in range(elem.n **2):
            a= np.array(self.Nx[i])[np.newaxis] #z wektora tworze dwuwymiarowa macierz, aby moc go transponowac
            b= np.array(self.Ny[i])[np.newaxis] # funkcja dodaje nowa oś (wymiar)

            px = a * a.T # T - funkcja, ktora dokonuje transpozycji
            py = b* b.T
            # macierz H dla każdego punktu całkowania:
            self.Hmac[i] =self.k* (px+py) * self.jakobian[i].det

        pom = 0
        for i in range(self.elem.n):
            for j in range(self.elem.n):
                self.H += self.Hmac[pom] * self.elem.waga[j] * self.elem.waga[i]
                pom += 1




