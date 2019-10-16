'''
Bio 331: HW 4
Author: Sol Taylor-Brill
Date: 10/16/19

The purpose of this code is to implement Kruskal's algorithm and post a graph to GraphSpace. 
'''

from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

def main():
	test_file = "weighted-graph.txt" #example graph
	real_file = "weighted-EGFR.txt" #EGFR graph
	
	edge_list = read_edge_file(real_file) #reads edgefile and makes edge_list (list of lists)
	T = kruskal(edge_list) #calculates minimum spanning tree T (a set of tuples)
	post_graph(edge_list,T) #Posts graph to graphspace
	return

#Input: file, a 3 column file (node A, node B, weight). 
#output: edge_list, a list of lists containing weighted edges.
def read_edge_file(file):
	f = open(file, 'r')
	edge_list = []
	for line in f:
		items = line.split()
		w = float(items[2])
		edge = [items[0],items[1],w]
		edge_list.append(edge)
	return edge_list

#input = a list of lists containing the two nodes in an edge and the weight of the edge
#output = T a set of tuples that are edges that make up the minimum spanning tree
def kruskal(edge_list):
	C = []
	T = set()
	for e in edge_list:
		nodeA = [e[0]] #single element list
		if nodeA not in C:
			C.append(nodeA)
		nodeB = [e[1]] #single element set
		if nodeB not in C:
			C.append(nodeB)
	Ordered = sorted(edge_list,key = lambda x: x[2]) #ordered by weight of the edge
	for e in Ordered:
		u = e[0] #first node in edge
		Cu = in_tree(u,C) #subtree that u is in
		v = e[1] #second node in edge
		Cv = in_tree(v,C) #subtree that v is in
		if Cu != Cv: #if the two trees are different
			T.add((u,v)) #adds a tuple to the set of edges that make up the tree
			Cu_set = set(Cu)
			Cv_set = set(Cv)
			newC = Cu_set.union(Cv_set) #combines the two sets
			listnewC = list(newC) #this just makes it easier to work with
			C.remove(Cu)
			C.remove(Cv)
			C.append(newC)
	return T

#just returns the subtree that a node is in
#input = v (a node) and C (a list of lists containing subtrees)
#output = the subtree that v is in
def in_tree(v,C):
	for subtree in C:
		if v in subtree:
			return subtree

#This function posts a graph to graphspace with the MST edges in blue
def post_graph(edge_list,T):
	graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session
	G = GSGraph()
	G.set_name('Sol_MST_EGFR_Graph')
	G.set_tags(['HW 4'])

	seennodes = set() #keeps track of nodes which have been added since nodes can only be added once
	for e in edge_list:
		nodeA = e[0]
		nodeB = e[1]
		if nodeA not in seennodes:
			G.add_node(nodeA, label=nodeA) #adds first node
			seennodes.add(nodeA)
		if nodeB not in seennodes:
			G.add_node(nodeB, label=nodeB) #adds second node
			seennodes.add(nodeB)

		G.add_edge(nodeA,nodeB) #adds edge
		edge_tup = (nodeA, nodeB)
		if edge_tup in T: 
			G.add_edge_style(nodeA, nodeB, color='blue')

	graphspace.post_graph(G)
	print("Graph Posted")
		

main()