import tkinter
import numpy as np
from gmesh_recmesh import makerec, reqrecMsh
from pydec import simplicial_complex

root = tkinter.Tk()
sz = 500

C = tkinter.Canvas(root, bg="blue", height=sz+10, width=sz+10)


V = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\vert.txt")*sz
E = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\tria.txt",dtype='int32')
V,E = makerec(sz,sz,sz/10)
sc = simplicial_complex(V,E)

edges = [[a[0],a[1]] for a in E]
edges.extend([[a[1],a[2]] for a in E])
edges.extend([[a[2],a[0]] for a in E])
edges = [[a[1],a[0]] if a[0]>a[1] else a for a in edges]
edges = np.unique(edges, axis=0)
ccs = 0.5*(V[edges[:,0]]+V[edges[:,1]])

sc[0].d.getcol(6).indices #enges of vertice 6?
sc[1].num_simplices #number of edges
[sc[0].simplex_to_index[i] for i in sc[1].index_to_simplex[34].boundary()] #boundary simplice indexes of edge num 34
# or just sc[1].simplices...
marks = np.zeros(sc[1].num_simplices)

def drawtria(tr):
    C.create_line(V[tr[0]][0], V[tr[0]][1], V[tr[1]][0], V[tr[1]][1])
    C.create_line(V[tr[1]][0], V[tr[1]][1], V[tr[2]][0], V[tr[2]][1])
    C.create_line(V[tr[2]][0], V[tr[2]][1], V[tr[0]][0], V[tr[0]][1])

def draw():
    # C.delete("all")
    # for tria in trias:
    #     drawtria(tria)
    for i,a in enumerate(sc[1].simplices):
        if marks[i]==1:
            C.create_line(V[a[0]][0],V[a[0]][1],V[a[1]][0],V[a[1]][1], fill="red")            
        else:
            C.create_line(V[a[0]][0],V[a[0]][1],V[a[1]][0],V[a[1]][1])

def flip(event):
    global E
    x = event.x
    y = event.y
    idx = np.sum(np.abs(ccs - [x,y])**2,axis=-1)
    b = idx.argmin()
    ts = [edges[b][0] in t and edges[b][1] in t for t in E]
    trias = E[ts]
    n1 = [e not in edges[b] for e in trias[0]]
    n2 = [e not in edges[b] for e in trias[1]]
    trias[0][trias[0]==edges[b][0]] = trias[1][n2]
    trias[1][trias[1]==edges[b][1]] = trias[0][n1]
    En = E[[not b for b in ts]]
    E = np.concatenate((En, trias))
    draw(E)

def detect(event):
    global marks
    edg = np.sum((sc[1].circumcenter-[event.x,event.y])**2, axis=-1).argmin()
    marks[edg]=1
    draw()

# C.bind("<Button-1>", flip)
C.bind("<Button-1>", detect)
draw()
C.pack()
root.mainloop()

