'''
Bio 331 HW 3: Random Walks
Author: Sol Taylor-Brill
Date: 10/03/19
'''

import matplotlib.pyplot as plt
import random
import math
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session

def main():
	foutput = FileReader("EGFR1-edges.txt")
	nodes = foutput[0]
	edges = foutput[1]
	interactDict = foutput[2]

	s="EGF" #the starting node
	Tprob = 1000
	Tsim = 100000
	q = 0.95

	#run either add_edges OR add_self_loops to get either out connections to source node or self loops
	Gin, Gout = Neighbor_Dicts(nodes,edges)
	#G_in, G_out = add_edges(Gin, Gout, nodes,edges,s)
	G_in, G_out = add_self_loops(Gin, Gout, nodes,edges) 

	X, TDict = rw_probs(nodes,G_in,G_out,s,Tprob,q)
	TimeDict = rw_simulate(nodes,G_out,s,Tsim,q)
	Graph_Maker(nodes,edges,TimeDict,interactDict)

	'''
	TpList = [10,100,1000,10000]
	TsList = [10,100,1000, 10000]
	ErrorList = []
	for i in range(len(TpList)):
		X, TDict = rw_probs(nodes,G_in,G_out,s,TpList[i],q)
		TimeDict = rw_simulate(nodes,G_out,s,TsList[i],q)
		NormTimeDict = Normalizer(TimeDict)
		Error = Error_Calc(nodes,NormTimeDict,TDict,TpList[i])
		ErrorList.append(Error)
	print(ErrorList)
	'''
	return

#This function just reads the edgefile. It's based on a function from HW2 but it uses lists instead of sets
#and doesn't include the edge type information
def FileReader(f):
	print("reading " + f)
	EdgeList = open(f,'r')
	nodes = []
	edges = []
	InteractDict = {}

	#The following for loop goes through all the lines in the file and adds nodes and edges to their lists
	#It also makes an edge->interaction dictionary
	for l in EdgeList:
		line = str(l)
		nodesinline = line.split()
		n1 = nodesinline[0]
		n2 = nodesinline[1]
		file_edge = [n1,n2] #doesn't include edge type
		callstring = str(n1 + n2)
		Interaction = nodesinline[2]
		InteractDict[callstring] = Interaction
		if file_edge not in edges:
			edges.append(file_edge)
		if n1 not in nodes:
			nodes.append(n1)
		if n2 not in nodes:
			nodes.append(n2)
	return nodes, edges, InteractDict

#makes Gin and Gout (in neighbors and out neighbors) dictionaries
#This function was copied from Lab 5
def Neighbor_Dicts(nodes,edges):
	Gin = {} #dictionary node->in neighbors [0]
	Gout = {} #dictionary node->out neighbors [1]
	for n in nodes:
		Gin[n] = []
		Gout[n] =[]
	for e in edges:
		NinList = Gin[e[1]] #calling all the in neighbors of the out node (second node)
		NoutList = Gout[e[0]] #calling all the out neighbors of the in node (first node)
		NinList.append(e[0]) #adding the first node (in node) to the list of in neighbors of the out node
		NoutList.append(e[1]) #adding the second node (out node) to the in nodes list of out neighbors
		Gin[e[1]] = NinList
		Gout[e[0]] = NoutList
	return Gin, Gout

#This functions adds extra edges to the graph so that any node that has no outgoing edges will have an edge to s
def add_edges(Gin,Gout,nodes,edges,s):
	for v in nodes: #goes through all the nodes in the graph
		if len(Gout[v]) == 0: #if the node has no outgoing neighbors
			edges.append((v,s))
			Gout[v] = [s] #the list of out-neighbors of v = s
			inputval = Gin[s] + [v] #adds v to the list of in-neighbors of s
			Gin[s] = inputval #updates in-neighbors of s
	return Gin, Gout

def add_self_loops(Gin,Gout,nodes,edges):
	for v in nodes: #goes through all the nodes in the graph
		if len(Gout[v]) == 0: #if the node has no outgoing neighbors
			edges.append((v,v))
			Gout[v] = [v] #the list of out-neighbors of v = v
			inputval = Gin[v] + [v] #adds v to the list of in-neighbors of v
	return Gin, Gout

# This function intakes a list of nodes, the Gin/Gout dictionaries, a starting node, s, timesteps, T, and probability value, q.
#This function was copied from Lab 5
def rw_probs(nodes,Gin,Gout,s,T,q):
	table_Dict = {} #a dictionary that will contain the dictionary of values for each time step
	table_list = [] #a list of lists

	#this part just initializes a dictionary of dictionaries so that it will run more easily.
	for t in range(0,(T+1)):
		nodedict = {}
		for n in nodes:
			nodedict[n] = []
		table_Dict[t]= nodedict


	for node in nodes:
		if node == s:
			table_Dict[0][node]=1
		else:
			table_Dict[0][node]=0
	for t in range(1,(T+1)):
		for n in nodes: #goes through all nodes
			val = 0
			for x in Gin[n]: #goes through all in-neighbors
				inputval = table_Dict[(t-1)][x]/len(Gout[x])
				val = val + inputval#making a sum of all p(u)^t-1
			eq_val = q*val + ((1-q)*table_Dict[0][n]) #probability value for each node in time, t
			table_Dict[t][n] = eq_val #adding the value to the table_Dict
	
	for call in table_Dict:
		tlist = []
		for val in table_Dict[call]:
			tlist.append(table_Dict[call][val])
		table_list.append(tlist) #making a list of lists
	return table_list, table_Dict

#This function simulates a random walker moving through a graph and 
#outputs a dictionary of times a node was visited by a random walker
#inputs = nodes, G_out (dictionary node->Nout), s (source node), T (time steps), q (probability [0-1])
def rw_simulate(nodes,G_out,s,T,q):
	TimeDict = {} #Dictionary node->times walker visited
	for n in nodes: #initializing TimeDict w/ all nodes except s at 0
		if n == s:
			TimeDict[n] = 1 #because it starts at 1
		else:
			TimeDict[n] = 0 #all other nodes initialized at 0

	#Moving through time steps
	CurrentNode = s #where the traveler is at. 
	for t in range(1,T):
		start = CurrentNode
		n_out = G_out[start] #a list of all outgoing nodes of the current node
		num = random.random() #generates random number between 0 and 1
		if num > q: 
			CurrentNode = s #redefining the current node. Teleport back to source
		else:
			newnode = random.choice(n_out) #choose random outgoing neighbor as new current node
			CurrentNode = newnode
		numtimes = TimeDict[CurrentNode] + 1 #adds 1 to the node that was just moved to
		TimeDict[CurrentNode] = numtimes #redefines node that just got moved to
	return TimeDict

#Normalizes the TimeDict values so that they are between 0 and 1
def Normalizer(TimeDict):
	TimeList = [] #a list of the number of time each node has been visited
	NormTimeDict = {}
	for node in TimeDict:
		TimeList.append(TimeDict[node])
	for n in TimeDict:
		NormTimeDict[n] = TimeDict[n]/max(TimeList)
	return NormTimeDict

#calculates the error (sum for all the nodes of the difference between the simulated value and calculated value)
def Error_Calc(nodes, NormTimeDict, TableDict, Tprob):
	sumval= 0
	for n in nodes:
		errorval = (NormTimeDict[n] - TableDict[Tprob][n])**2 #squared difference
		sumval = sumval + errorval #sum of squared differences
	print(sumval)
	return sumval

#Makes the graph.
#Taken from Lab 4
def Graph_Maker(node_list, edge_list, TimeDict,interactDict):
	G = GSGraph()
	G.set_name('Sol_Self-Loop_Walk_Graph')
	G.set_tags(['HW3'])

	#This section transforms and normalizes the data
	Timelist = [] #I'm just making a list of counts to use to calculate node color
	logList = []
	NewTimeDict = {}
	for noden in TimeDict: #just making a list of values in the dictionary (useful for nodecolor)
		newval = TimeDict[noden] + 1 #adds 1 to each countso there are no zeros
		NewTimeDict[noden] = newval #adds the value
		if NewTimeDict[noden] not in Timelist:
			Timelist.append(NewTimeDict[noden]) #adds the value to a list
	for val in Timelist: #to figure out the max log transformed value
		logv = math.log(val)
		logList.append(logv)
	for n in NewTimeDict: #log transforms values and then normalizes
		logval = math.log(NewTimeDict[n])
		NewTimeDict[n] = logval/max(logList) #new value is the log of the count normalized by dividing by the max log value
	
	#This section makes the graph and updates on GraphSpace
	seenedges = [] # a list of edges that have gone over so there aren't multiple edges between two nodes
	for n in node_list:
		poplabel = "Count: " + str(TimeDict[n]) + ";Normalized: " + str(NewTimeDict[n])
		G.add_node(n,popup=poplabel, label=n)
		blue = NewTimeDict[n] #a value between 0 and 1 (most visited = most blue)
		red = 1-blue #opposite of blue
		nodecolor = rgb_to_hex(red,0,blue)
		G.add_node_style(n, width = 20, height = 20, color=nodecolor,)
	for e in edge_list:
		namecall = str(e[0] + e[1]) #to call up interaction type for the popup
		if namecall in interactDict:
			plabel = interactDict[namecall]
		else:
			plabel = "N/A"
		if [e[1],e[0]] in edge_list: #if the inverse exists
			if e not in seenedges:
				G.add_edge(e[0],e[1], popup = plabel, directed=False)
			seenedges.append(e)
			seenedges.append([e[1],e[0]])
		else: #if the graph isn't bidirectional add an arrow
			G.add_edge(e[0],e[1],popup = plabel, directed=True)
			G.add_edge_style(e[0],e[1], width = 2,directed=True)
	graphspace.post_graph(G)
	print("Graph Updated")

'''
This function just outputs a hex color code from inputs of the amount of red and blue calculated from layer in Graph_Maker()
This function was just copied from the lab 4 handout.
'''
def rgb_to_hex(red,green,blue):
	maxHexValue = 255
	r = int(red*maxHexValue)
	g = int(green*maxHexValue)
	b = int(blue*maxHexValue)
	RR = format(r,'02x')
	GG = format(g,'02x')
	BB = format(b,'02x')
	colorname = '#' + RR + GG + BB
	return colorname

main()