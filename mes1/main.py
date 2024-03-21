import numpy as np
from dane import *
from elem_uniwersalny import *
from Jakobian import *
from macierz_HBC import *
from obliczanie_dla_elementow import *
from macierz_C import *
from globalne import *
import csv

# j=Wczytaj("Test1_4_4.txt")
# j = Wczytaj("Test2_4_4_MixGrid.txt")
j=Wczytaj("Test3_31_31_kwadrat.txt")
# j=Wczytaj("Test4_31_31_trapez.txt")
g, gd = j.wczytaj()

def f(x):
    return 5*x**2+3*x+6

def f2(x,y):
    return 5*x**2*y**2 + 3*x*y + 6

# print(g.BC)
# print(g.Nodes)
# print(g.Elements)
# print(gd.data)

#GAUSS
# c = gauss(3,1)
# print(c.gauss1d(f))

# d = gauss(2,2)
# print(d.gauss2d(f2))


# e = ElementUniwersalny(2)

# x = [0,0.025,0.025,0]
# y = [0,0,0.025,0.025]
# j = Jakobian(x,y,e,0)

# print(j.jac)
# print(j.inv)

# hbc = HBC(gd,e,[1,1,1,1],g)
# pom = [g.Elements[2].id[2],g.Elements[2].id[3],g.Elements[2].id[0],g.Elements[2].id[1]]
# C=C(x,y,e,gd)
# hbc.sciany(np.array(pom))
# hbc.oblicz()
# dlugosc_odc()

# x,y = punktowanie(g,gd)
# h,c = MacierzeHiC(x,y,e,gd)
# hbc,p = HBCiP(gd,g,e)
# print("------------H------------:")
# for i in h:
#     print(i)
#     print()
# print("------------C------------:")
# for i in c:
#     print(i)
#     print()
# print("-----------HBC------------:")
# for i in hbc:
#     print(i)
#     print()
# print("-----------P------------:")
# for i in p:
#     print(i)
#     print()
# print(globMatrix(g,gd,2)[2])
# print()
temp =[float(gd.data["InitialTemp"]) for i in range(len(g.Nodes))]
# s = solve(g,gd,2,temp)
tab = np.zeros([21, len(g.Nodes)])
tab[0] = temp
for i in range(1,21):
    temp = solve(g,gd,2,temp)
    tab[i]=temp
# plikvtk(g,tab[0 ], "test.vtk")
# wynikiminmax(tab)
# for i in range(11):
#     plikvtk(g,tab[i],f"./vtk/test{i}.vtk")
for i in tab:
    print(f"{max(i)} {min(i)}")


