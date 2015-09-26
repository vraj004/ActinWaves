#!/usr/bin/python
import numpy as np
import scipy as sp
import math as mth

def createMito(radii=[0.25,0.25,1],centre=[0.0,0.0,0.0],discretization=[2,3],node_offset=0,facet_marker_offset=0):
	"#A function to create a mito poly file of a user-defined size and at a user-defined location"
	#script to create an ellipsoidal poly file
	a = radii[0] #x-axis radius
	cx = centre[0] #x-center
	b = radii[1] #y-axis radius
	cy = centre[1] #y-center
	c = radii[2] #z-axis radius
	cz = centre[2] #z-center
	numz_steps = discretization[1]
	numx_steps=discretization[0]

	#ellipsoid eqn - x^2/a^2 + y^2/b^2 + z^2/c^2 = 1
	
	#parameterised eqn:
	#x = a cosbeta coslambda
	#y = b cosbeta sinlambda
	#z = c sinbeta
	#-pi/2<beta<pi/2; -pi<lambda<pi

	#loop over z from -1 to 1
	node=1+node_offset
	nodeind = 0
	mitoNodes=[]
	mitoNodes.append(np.array([node,(0+cx),(0+cy),(-1.0*c+cz)]))
	nodeind+=1

	z_step = 2.0*c/numz_steps
	z = -1.0*c+z_step

	while z<(c-1e-09):
		r = mth.sqrt((1-(z**2/c**2))*a**2)
		x = -1*r
		x_step = r*2/numx_steps
		while (x<=r):
			ysq = (((1-(z**2/c**2))*a**2)-x**2)
			if(ysq<0):
				print ysq
				if mth.fabs(ysq)<0.00001:
					y = 0
				else:
					die
			else:
				y = mth.sqrt(ysq)	
			node+=1
			mitoNode = np.array([node,(x+cx),(y+cy),(z+cz)])
			mitoNodes.append(mitoNode)
			nodeind+=1
			x+=x_step
		x = r-x_step
		while (x>-1*r):
			y = -1*mth.sqrt(((1-(z**2/c**2))*a**2)-x**2	)
			node+=1
			mitoNode = np.array([node,(x+cx),(y+cy),(z+cz)])
			mitoNodes.append(mitoNode)
			nodeind+=1
			x-=x_step
		z=z+z_step
	#adding the last point in the ellipsoid at the other end at z = 1
	node+=1
	mitoNode = np.array([node,0+cx,0+cy,(c+cz)])
	mitoNodes.append(mitoNode)
	nodeind+=1	

####create face segment list.#######
	numFacesperZsec = numx_steps*2 #This is the number of facets for each z-section. If you have 2 x steps, it is 4 dots on the circle. This means there will be 4 facets that make up a surface around these dots
	numZsec = numz_steps 
	totalFacets = numZsec*numFacesperZsec
	facetList = [0]*totalFacets
	facet_marker=facet_marker_offset+1
	zsec = 1
	collapsednode = int(mitoNodes[0][0])
	endnodeInd = 1
	facetnumber=0
	#set up the first set of facets.
	for facetCntr in range(1,numFacesperZsec):
		facetList[facetnumber] = np.array([str(3),str(collapsednode),str(int(mitoNodes[endnodeInd+facetCntr-1][0])),str(int(mitoNodes[endnodeInd+facetCntr][0])),str(facet_marker)]) #adding boundary marker label at the end of the list
		facetnumber+=1

	facetCntr = numFacesperZsec
	facetList[facetnumber] = np.array([str(3),str(collapsednode),str(int(mitoNodes[endnodeInd+facetCntr-1][0])),str(int(mitoNodes[endnodeInd][0])),str(facet_marker)])
	facetnumber+=1

	#fill out the facet list for the mid-z-sections.
	strtnodeInd = 1
	endnodeInd = strtnodeInd+numFacesperZsec
	for zsec in range(2,numz_steps):
		for facetCntr in range(1,(numFacesperZsec)): #just do the first three facets in the zsection - the last one is a bit special.
			facetList[facetnumber] = np.array([str(4),str(int(mitoNodes[strtnodeInd+facetCntr-1][0])),str(int(mitoNodes[endnodeInd+facetCntr-1][0])),str(int(mitoNodes[endnodeInd+facetCntr][0])),str(int(mitoNodes[strtnodeInd+facetCntr][0])),str(facet_marker)])
			facetnumber+=1

		facetCntr = numFacesperZsec
		facetList[facetnumber] = np.array([str(4),str(int(mitoNodes[strtnodeInd+facetCntr-1][0])),str(int(mitoNodes[endnodeInd+facetCntr-1][0])),str(int(mitoNodes[endnodeInd][0])),str(int(mitoNodes[strtnodeInd][0])),str(facet_marker)])
		facetnumber+=1
		strtnodeInd = endnodeInd
		endnodeInd = strtnodeInd+numFacesperZsec

	zsec = numz_steps #Let's look at the last z-section as this is special like the first z-section.
	collapsednode = int(mitoNodes[len(mitoNodes)-1][0]) #last node is the collapsing node on the other end.
	#set up the last set of facets.
	for facetCntr in range(1,numFacesperZsec):
		facetList[facetnumber] = ([str(3),str(int(mitoNodes[strtnodeInd+facetCntr-1][0])),str(collapsednode),str(int(mitoNodes[strtnodeInd+facetCntr][0])),str(facet_marker)])
		facetnumber+=1

	facetCntr = numFacesperZsec
	facetList[facetnumber] = ([str(3),str(int(mitoNodes[strtnodeInd+facetCntr-1][0])),str(collapsednode),str(int(mitoNodes[strtnodeInd][0])),str(facet_marker)])
	return {'mitoNodes':mitoNodes,'mitoFacets':facetList}

def writePoly(polyfile,nodes,facets,holes):
	"A function to write out a poly file containing the nodes and facts of the cell and its mitos"
	#write out the nodes into the poly file

	cellPoly = open(polyfile,'w+')
	cellPoly.write('#Part 1 - node list\n')
	cellPoly.write('#node count, 3 dim, no attribute, boundary markers on\n')
	nodeHeader = [str(len(nodes)),str(3),str(0),str(1)]
	cellPoly.write(' '.join(nodeHeader))
	cellPoly.write('\n')
	for node in range(len(nodes)):
		NodeArr = nodes[node]
		output=[str(int(NodeArr[0])),str(NodeArr[1]),str(NodeArr[2]),str(NodeArr[3]),str(0)] #last entry is for boundary markers - set to 0 so that tetgen can assign.
		cellPoly.write(' '.join(output))
		cellPoly.write('\n')
		print output
# write out facets into the poly file
	totalFacets = len(facets)
	cellPoly.write('#Part 2 - facet list\n')
	cellPoly.write('#facet count, boundary markers on\n')
	cellPoly.write(str(totalFacets)+' '+str(1)+'\n')
	cellPoly.write('#facets\n')
	for facet in range(len(facets)):
		facetArr = facets[facet]
		facetArrSize = len(facetArr)
		output = [facetArr[0]]
		for ind in range(1,facetArrSize-1):
			output.append(facetArr[ind])
		print output
		
		cellPoly.write(str(1)+' '+str(0)+' '+facetArr[facetArrSize-1]+'\n')
		
#		if int(facetArr[0]==str(3)):
#			output = [str((facetArr[0])),str((facetArr[1])),str((facetArr[2])),str((facetArr[3]))]
#		elif int(facetArr):
#			output = [str((facetArr[0])),str((facetArr[1])),str((facetArr[2])),str((facetArr[3])),str((facetArr[4]))]
		cellPoly.write(' '.join(output))
		cellPoly.write('\n')
	cellPoly.write('#part 3 - hole list\n')
	numHoles = len(holes)
	if(numHoles==0):
		cellPoly.write(str(0))
		cellPoly.write('\n')
	else:
		cellPoly.write(str(numHoles))
		cellPoly.write('\n')
		for hole in range(len(holes)):
			holecoords = holes[hole]
			output = [str(hole+1),str(holecoords[0]),str(holecoords[1]),str(holecoords[2])]
			cellPoly.write(' '.join(output))
			cellPoly.write('\n')
	
	cellPoly.write('#part 4 - region list\n')
	cellPoly.write(str(0)+' #no region\n')	
		
		
	cellPoly.close()
	return

def readPolyNodes(polyfile):
	"A function to read in a poly file. This is initially being written to write out boundary condition files (grouping of nodes for BCs)"
	
	cellPoly = open((polyfile+'.node'),'r')
	polyData = cellPoly.readlines()
	nodeHeader = (polyData[0]).rsplit()
	numNodes = int(nodeHeader[0])
	dim = int(nodeHeader[1])
	boundaryMarker = int(nodeHeader[3])
	nodes = [0]*numNodes
	for line in range(1,numNodes+1):
		node = polyData[line].rsplit()
		nodes[line-1] = list(node)
	return {'numNodes': numNodes,'boundaryMarker': boundaryMarker, 'nodes': nodes}
	
