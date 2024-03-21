import numpy as np
from elem_uniwersalny import ElementUniwersalny

class Jakobian:
    def __init__(self, x: [], y: [], e: ElementUniwersalny, punkt: int):
        self.x = np.array(x)
        self.y = np.array(y)
        # interpolacja
        self.X_Ksi = np.sum([self.x[i] * e.nKsi[punkt][i] for i in range(np.size(self.x))])
        self.X_Eta = np.sum([self.x[i] * e.nEta[punkt][i] for i in range(np.size(self.x))])
        self.Y_Ksi = np.sum([self.y[i] * e.nKsi[punkt][i] for i in range(np.size(self.x))])
        self.Y_Eta = np.sum([self.y[i] * e.nEta[punkt][i] for i in range(np.size(self.x))])

        self.jac = np.array([[self.X_Ksi,self.Y_Ksi],[self.X_Eta,self.Y_Eta]]) #macierz jacobiego
        self.det = np.linalg.det(self.jac)  #wyznacznik macierzy
        self.inv = np.linalg.inv(self.jac) #odwrócona macierz jacobiego
        print("inv")
        print(self.inv)



