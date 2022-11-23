import tkinter
import numpy as np
from gmesh_recmesh import makerec, reqrecMsh
from pydec import simplicial_complex

root = tkinter.Tk()
sz = 500

C = tkinter.Canvas(root, bg="white", height=sz+10, width=sz+10)

# V = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\vert.txt")*sz
# E = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\tria.txt",dtype='int32')
stretch  = 2
V,E = reqrecMsh(sz/stretch,sz,sz/5)
V[:,[1,0]]=V[:,[0,1]]
V[:,1]*=stretch
sc = simplicial_complex(V,E)


marks = np.zeros(sc[1].num_simplices, dtype='int32')

def draw():
    o = 5
    for i,a in zip(marks, sc[1].simplices):
        colors = ["black","blue","orange","red"]
        C.create_line(V[a[0]][0]+o,V[a[0]][1]+o,V[a[1]][0]+o,V[a[1]][1]+o, fill=colors[i],  width=5)            

def updateAll(marks):
    for i in range(sc[1].num_simplices):
        updateComplete(i, marks)

def updateComplete(edg, marks):
    #nodes that are only missing one edge
    for vert in sc[1].simplices[edg]:
        vertEdges = sc[0].d.getcol(vert).indices
        if np.count_nonzero(marks[vertEdges]==3) == vertEdges.size-1:
            marks[vertEdges] = np.maximum(marks[vertEdges], 2*np.ones(vertEdges.size))
    #faces that are only missing one edge
    for vert in sc[1].d.getcol(edg).indices:
        faceEdges1 = sc[2].index_to_simplex[vert].boundary()
        edges1 = [sc[1].simplex_to_index[i] for i in faceEdges1]
        if np.count_nonzero(marks[edges1]==3) == 2:
            marks[edges1] = np.maximum(marks[edges1], np.ones(3))

def detect(event):
    global marks
    edg = np.sum((sc[1].circumcenter-[event.x,event.y])**2, axis=-1).argmin()
    marks[edg]=3
    updateComplete(edg, marks)
    draw()

def iterate(event):
    global marks
    marks = 3*(marks>0)
    updateAll(marks)
    draw()

# C.bind("<Button-1>", flip)
C.bind("<Button-3>", iterate)
C.bind("<Button-1>", detect)
draw()
C.pack()
root.mainloop()


# sc[0].d.getcol(6).indices #enges of vertice 6?
# sc[1].num_simplices #number of edges
# [sc[0].simplex_to_index[i] for i in sc[1].index_to_simplex[34].boundary()] #boundary simplice indexes of edge num 34

# def flip(event):
#     global E
#     x = event.x
#     y = event.y
#     idx = np.sum(np.abs(ccs - [x,y])**2,axis=-1)
#     b = idx.argmin()
#     ts = [edges[b][0] in t and edges[b][1] in t for t in E]
#     trias = E[ts]
#     n1 = [e not in edges[b] for e in trias[0]]
#     n2 = [e not in edges[b] for e in trias[1]]
#     trias[0][trias[0]==edges[b][0]] = trias[1][n2]
#     trias[1][trias[1]==edges[b][1]] = trias[0][n1]
#     En = E[[not b for b in ts]]
#     E = np.concatenate((En, trias))
#     draw(E)

# def drawtria(tr):
#     a = 15
#     C.create_line(V[tr[0]][0]+a, V[tr[0]][1]+a, V[tr[1]][0]+a, V[tr[1]][1]+a)
#     C.create_line(V[tr[1]][0]+a, V[tr[1]][1]+a, V[tr[2]][0]+a, V[tr[2]][1]+a)
#     C.create_line(V[tr[2]][0]+a, V[tr[2]][1]+a, V[tr[0]][0]+a, V[tr[0]][1]+a)

