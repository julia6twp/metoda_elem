import re

import numpy as np


class GlobalData:
    def __init__(self, pom : dict):
        self.data = pom

class Node:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return repr([self.x,self.y])

class Element:
    def __init__(self, element : np.array(int)):
        self.id = np.array(element)

    def __repr__(self):
        return repr(self.id)


class Grid:
    def __init__(self, x:np.array(Node), y:np.array(Element), bc:np.array(int)):
        self.Nodes = x
        self.Elements = y
        self.BC = bc  #konce odcinka na ktory nakladamy war. brzegowy

class Wczytaj:
    def __init__(self, s):
        self.path = s
        self.nodes = []
        self.elements = []
        self.bc = []
        self.current_section = None
        self.gd = dict()
    def wczytaj(self):
        with open(self.path,"r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("*"):
                    self.current_section = line[1:]
                elif self.current_section == "Node":
                    values = re.split(r"\s+",line.strip())
                    cv = [float(value.rstrip(",")) for value in values]
                    self.nodes.append(cv)
                elif self.current_section == "Element, type=DC2D4":
                    values = re.split(r"\s+",line.strip())
                    cv = [int(value.rstrip(",")) for value in values]
                    self.elements.append(cv)
                elif self.current_section == "BC":
                    values = re.split(r"\s+", line.strip())
                    cv = [float(value.rstrip(",")) for value in values]
                    self.bc = cv
                elif self.current_section == None:
                    values = re.split(r"\s",  line.strip())
                    if values[0] == "Elements" or values[0] == "Nodes":
                        self.gd[values[0] + " " + values[1]] = float(values[2])
                    else:
                        self.gd[values[0]] = float(values[1])

        file.close()

        pomnode = []
        pomelement = []

        for i in self.nodes:
            pomnode.append(Node(i[1], i[2]))
        for i in self.elements:
            pomelement.append(Element(i[1:]))

        grid = Grid(pomnode, pomelement, np.array([int(self.bc[i]) for i in range(len(self.bc))]))
        globalData  = GlobalData(self.gd)

        return grid, globalData


def plikvtk(g:Grid, temp, filename):
    tab =[
       ("# vtk DataFile Version 2.0"),
        ("MES Results"),
                ("ASCII"),
          ("DATASET UNSTRUCTURED_GRID"),
          (""),
   ]
    le = len(g.Nodes)
    dod="POINTS "+str(le)+" float"
    tab.append(dod)
    for i in range(len(g.Nodes)):
       t = str(g.Nodes[i].x)+" "+str(g.Nodes[i].y)+" 0.0"
       tab.append(t)
    tab.append("")
    le= len(g.Elements)
    dod= "CELLS "+str(le)+" "+str(le*5)
    tab.append(dod)
    tpom=""
    for i in range(len(g.Elements)):
        for j in g.Elements[i].id:
            tpom=tpom +" "+ str(j-1)
        t = "4"+ tpom
        tpom=""
        tab.append(t)
    tab.append("")
    dod="CELL_TYPES"+str(len(g.Elements))
    tab.append(dod)
    for i in range(len(g.Elements)):
        tab.append("9")
    le = len(g.Nodes)
    tab.append("")
    dod = "POINT_DATA "+str(le)
    tab.append(dod)
    tab.append("SCALARS Temp float 1")
    tab.append("LOOKUP_TABLE default")
    for i in range(le):
        tab.append(str(temp[i]))
    # print(t)
    with open(filename, 'w') as vtkfile:
        for i in tab:
           vtkfile.write(i+"\n")
    print(tab)

