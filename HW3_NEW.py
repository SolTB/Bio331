'''
Bio 331 HW 3: Random Walks
Author: Sol Taylor-Brill
Date: 10/03/19
'''

import matplotlib.pyplot as plt

def main():
	foutput = FileReader("test-edges1.txt")
	nodes = foutput[0]
	edges = foutput[1]
	s="A" #the starting node
	T = 100
	q = 0.5
	Gin, Gout = Neighbor_Dicts(nodes,edges)
	G_in, G_out = add_edges(Gin, Gout, nodes,edges,s)
	X = rw_probs(nodes,G_in,G_out,s,T,q)
	return

#This function just reads the edgefile. It's based on a function from HW2 but it uses lists instead of sets
def FileReader(f):
	print("reading " + f)
	EdgeList = open(f,'r')
	nodes = []
	edges = []

	#The following for loop goes through all the lines in the file and adds nodes and edges to their sets
	for l in EdgeList:
		line = str(l)
		nodesinline = line.split()
		if nodesinline not in edges:
			edges.append(nodesinline)
		for n in nodesinline:
			if n not in nodes:
				nodes.append(n)
	return nodes, edges

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

# This function intakes a list of nodes, the Gin/Gout dictionaries, a starting node, s, timesteps, T, and probability value, q.
#This function was copied from Lab 5
def rw_probs(nodes,Gin,Gout,s,T,q):
	table_Dict = {} #a dictionary that will contain the dictionary of values for each time step
	table_list = [] #a list of lists

	#this part just initializes a dictionary of dictionaries so that it will run more easily.
	for t in range(0,T):
		nodedict = {}
		for n in nodes:
			nodedict[n] = []
		table_Dict[t]= nodedict


	for node in nodes:
		if node == s:
			table_Dict[0][node]=1
		else:
			table_Dict[0][node]=0
	for t in range(1,T):
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
	return table_list

def rw_simulate():
	return

main()