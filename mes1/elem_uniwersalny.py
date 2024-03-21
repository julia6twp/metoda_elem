import numpy as np
def fksi1(eta):
    return -1/4*(1-eta)
def fksi2(eta):
    return 1/4*(1-eta)
def fksi3(eta):
    return 1/4*(1+eta)
def fksi4(eta):
    return -1/4*(1+eta)

def feta1(ksi):
    return -1/4*(1-ksi)
def feta2(ksi):
    return -1/4*(1+ksi)
def feta3(ksi):
    return 1/4*(1+ksi)
def feta4(ksi):
    return 1/4*(1-ksi)

class ElementUniwersalny: #pochodne w pkt - Funkcje kształtu elementu czterowęzłowego
    def __init__(self, n:int):
        self.n = n
        self.tab, self.waga = np.polynomial.legendre.leggauss(n)
        self.nKsi = np.zeros([n**2, 4]) #nazwy
        self.nEta = np.zeros([n**2, 4]) #
        self.pom = 0                     #
        for i in range(n**2):
            if (i != 0 and i % n == 0): self.pom += 1
            self.nKsi[i][0] = fksi1(self.tab[self.pom])
            self.nKsi[i][1] = fksi2(self.tab[self.pom])
            self.nKsi[i][2] = fksi3(self.tab[self.pom])
            self.nKsi[i][3] = fksi4(self.tab[self.pom])

        self.pom = 0
        for i in range(n**2):
            self.pom = i%n
            self.nEta[i][0] = feta1(self.tab[self.pom])
            self.nEta[i][1] = feta2(self.tab[self.pom])
            self.nEta[i][2] = feta3(self.tab[self.pom])
            self.nEta[i][3] = feta4(self.tab[self.pom])
