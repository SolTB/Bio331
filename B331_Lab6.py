'''
Author: Sol Taylor-Brill
Date: 10/30/19
Project: Lab 6: Network Motif Finding
'''

from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session

import itertools
import random

def main():
	toyfile = "yeast-GRN-full.txt"
	nodelist, edgelist, nodestring = TextReader(toyfile)
	AR_List, AR_count = count_AR(nodelist,edgelist)
	MCL_List, MCL_count = count_MCL(nodelist,edgelist)
	FFL_Count = count_FFL(nodelist,edgelist)
	p_AR, p_MCL, p_FFL = p_function(nodelist, edgelist, 20, 20)
	Graph_Maker(nodelist, edgelist, MCL_List, p_AR, p_MCL, p_FFL, MCL_count)
	return

def TextReader(file):
	f=open(file,'r')
	print("Opening " + file)
	edgelist = []
	nodelist = []
	nodestring = ""
	for line in f:
		l = str(line)
		edgeuntrimmed = l.split(';') #splits the line by semicolon and makes a list

		nodeA = edgeuntrimmed[0] #the first node
		nodeB = edgeuntrimmed[1][0:4] #second node - trimmed
		edge = (nodeA,nodeB) #a tuple containing the edge without "\n"
		if edge not in edgelist:
			edgelist.append(edge) #appends the edge which is a tuple to the edgelist
		if nodeA not in nodelist:
			nodelist.append(nodeA)
			nodestring = nodestring + str(nodeA)
		if nodeB not in nodelist:
			nodelist.append(nodeB)
			nodestring = nodestring + str(nodeB)
	return nodelist, edgelist, nodestring

def count_FFL(nodelist,edgelist):
	InDict = {} #list of in neighbors of a node
	OutDict = {} #list of out neighbors of a node
	AllNDict = {} #dictionary of all in and out neighbors of a node

	for node in nodelist:
		OutList = []
		InList = []
		for edge in edgelist:
			if edge[0] == node: #if its the starting node
				OutList.append(edge[1]) #adding second node to list of out neighbors
			elif edge[1] == node: #if its the receiving node
				InList.append(edge[0]) #adding first node to list of in neighbors
		InDict[node] = InList
		OutDict[node] = OutList
		AllNList = InList + OutList
		AllNDict[node] = AllNList

	triples = itertools.combinations(nodelist,3)
	FFL_Count = 0
	for item in triples:
		A = item[0]
		AN = AllNDict[A] #all neighbors of A
		B = item[1]
		BN = AllNDict[B] #all neighbors of B
		C = item[2]
		CN = AllNDict[C] #all neighbors of C
		if B in AN and C in AN and B in CN: #if B and C are neighbors of A and B is a neighbor of C (ie edges exist between all)
			if C in OutDict[A]: #if there exists A->C
				if B in OutDict[A]: #if there exists A->B
					FFL_Count = FFL_Count + 1
			elif C in OutDict[B]: #if there exists B->C
				if A in OutDict[B]: #if there exists B->A
					FFL_Count = FFL_Count + 1
			elif B in OutDict[C]: #if there exists C->B
				if A in OutDict[C]: #if there exists C->A
					FFL_Count = FFL_Count + 1
	return FFL_Count

def count_AR(nodelist,edgelist):
	count = 0
	AR_List = []
	for node in nodelist:
		selfedge = (node,node)
		if selfedge in edgelist:
			count = count +1
			AR_List.append(selfedge)
	return AR_List, count

def count_MCL(nodelist,edgelist):
	pairs = itertools.combinations(nodelist,2)
	MCL_List =[]
	MCL_count = 0
	for pair in pairs:
		A = pair[0]
		B = pair[1]
		if (A,B) in edgelist and (B,A) in edgelist:
			MCL_count = MCL_count + 1
			MCL_List.append((A,B))
			MCL_List.append((B,A))
	return MCL_List, MCL_count

def rewire_graph(nodelist,edgelist,r):
	new_edges = []
	edgeset = set(edgelist)
	seenedges = []
	for i in range(r):
		pair = random.sample(edgeset,2) #two random edges (ie [(A,B),(C,D)])
		seenedges.append(pair[0])
		seenedges.append(pair[1])
		A = pair[0][0]
		B = pair[0][1]
		C = pair[1][0]
		D = pair[1][1]
		newedgeA = (A,D) #makes new edge (A,D)
		newedgeB = (C,B) #makes new edge (C, B)
		if newedgeA not in new_edges and newedgeB not in new_edges:
			new_edges.append(newedgeA)
			new_edges.append(newedgeB)
		for edge in edgelist:
			if edge not in seenedges:
				new_edges.append(edge)
	return new_edges

def p_function(nodelist, edgelist, r, t):
	OAR, OriginalAR = count_AR(nodelist,edgelist)
	OMCL, OriginalMCL = count_MCL(nodelist,edgelist)
	OriginalFFL = count_FFL(nodelist,edgelist)

	print("ORIGINAL VALUES. AR: ", OriginalAR, "MCL: ", OriginalMCL, "FFL: ", OriginalFFL)
	Master_Graph_List = [] #list containing all the new graphs (only the edges)
	FFL_CountList = []
	AR_CountList = []
	MCL_CountList = []
	for i in range(t):
		new_edges = rewire_graph(nodelist,edgelist,r)
		Master_Graph_List.append(new_edges)
	for i in range(len(Master_Graph_List)):
		FFL = count_FFL(nodelist,Master_Graph_List[i])
		FFL_CountList.append([FFL])
		MList, MCL = count_MCL(nodelist,Master_Graph_List[i])
		MCL_CountList.append([MCL])
		AList, AR = count_AR(nodelist,Master_Graph_List[i])
		AR_CountList.append([AR])
	
	print("NEW AR Values: ", AR_CountList)
	print("NEW FFL Values: ", FFL_CountList)
	print("NEW MCL Values: ", MCL_CountList)
	
	ARGreater = 0
	MCLGreater = 0
	FFLGreater = 0
	for i in range(t):
		AR = int(AR_CountList[i][0])
		MCL = int(MCL_CountList[i][0])
		FFL = int(FFL_CountList[i][0])
		if AR >= OriginalAR:
			ARGreater = ARGreater +1
		if MCL >= OriginalMCL:
			MCLGreater = MCLGreater +1
		if FFL >= OriginalFFL:
			FFLGreater = FFLGreater + 1


	print("AR COUNT", ARGreater, "FFL COUNT", FFLGreater, "MCL COUNT", MCLGreater)

	p_AR = ARGreater/t
	p_MCL = MCLGreater/t
	p_FFL = FFLGreater/t

	print(p_AR, p_MCL, p_FFL)

	return p_AR, p_MCL, p_FFL


def Graph_Maker(nodelist, edgelist, MCL_List, p_AR, p_MCL, p_FFL, MCL_count):
	G = GSGraph()
	G.set_name('Sol_YEAST_MCL_FULL_Graph')
	G.set_tags(['Lab 6'])
	desc_string = 'Fi MCL = ' + str(MCL_count) + ';r=20; t=20; p(AR) = ' + str(p_AR) + ' p(MCL) =' + str(p_MCL) + ' p(FFL) = ' + str(p_FFL)
	G.set_data(data={'description': desc_string})
	for n in nodelist:
		G.add_node(n, label=n)
	for e in edgelist:
		G.add_edge(e[0],e[1])
		if e in MCL_List:
			G.add_edge_style(e[0],e[1], width = 2,directed=True, color='yellow')
		else:
			G.add_edge_style(e[0],e[1], width = 2, directed=True)
	graphspace.post_graph(G)
	print("Graph Updated")

main()