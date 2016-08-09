# -*- coding: utf-8 -*-



"""

Created on Mon Nov 16 16:07:24 2015



@author: vrajagopal and maxmillian bode



"""
import numpy as np
polyfile = 'ellipse.poly'
minor_axis = 1.0
major_axis = 1.0
num_heminodes = 50 #number of nodes to make one half of the ellipse
#create array of x-coordinates from -minor->+minor
x_cords_top = np.linspace((-1*minor_axis), (minor_axis),num_heminodes)
x_cords_bot = x_cords_top[1:(num_heminodes-1)]
x_cords_bot = np.fliplr([x_cords_bot])[0]
x_cords = np.concatenate((x_cords_top,x_cords_bot),axis=0)

#calcualte corresponding y-coords based on equation of ellipse
y_cords_top = []
for i in range(len(x_cords_top)):
    y = ((1.0-((x_cords_top[i]**2.0))/(minor_axis**2.0))*(major_axis**2))**(1.0/2.0)
    y_cords_top.append(y)

#create vector of mirrored nodes about x-axis for bottom of ellipse
y_cords_top=np.array(y_cords_top)
y_cords_bot = (y_cords_top[1:(num_heminodes-1)]) *-1
y_cords_bot = np.fliplr([y_cords_bot])[0]
y_cords = np.concatenate((y_cords_top,y_cords_bot))

nodes = len(x_cords)
numbering = range(1 ,(nodes+1))
numbering = np.array(numbering)
zeros = [0]*nodes


#write to file 
f = open(polyfile, 'w')
nodesheader_line = [(nodes), 2, 0, 0]
nodesheader_line = '   '.join(map(str,nodesheader_line))
f.write(nodesheader_line)
f.write('\n')

list_nodes = np.matrix((numbering, x_cords, y_cords, zeros))
list_nodes = list_nodes.transpose()
templist = []
for row in range(0,len(list_nodes)):
    temp = [str(int(list_nodes[row, 0])), str(list_nodes[row, 1]), str(list_nodes[row, 2]), str(int(list_nodes[row, 3]))]
    templist.append(temp)
    f.write('   '.join(temp))
    f.write('\n')

facetheader_line = [nodes, 0]
facetheader_line = '   '.join(map(str, facetheader_line))
f.write(facetheader_line)
f.write('\n')
facet_numbers = np.array(range(1,(nodes+1)))
facetnodei = np.array(range(1,(nodes+1)))
facetnodej = np.array(range(2, (nodes+1)))
facetnodej = np.append(facetnodej,1)
facetbdmarker = [0]*(nodes)
facetlist = np.matrix((facet_numbers, facetnodei, facetnodej, facetbdmarker))
facetlist = facetlist.transpose()
templist = []
for row in range(0,len(facetlist)):
    temp = [str(int(facetlist[row, 0])), str(int(facetlist[row, 1])), str(int(facetlist[row, 2])), str(int(facetlist[row, 3]))]
    templist.append(temp)
    f.write('   '.join(temp))
    f.write('\n')

f.write('0')
f.close()
