import statistics 
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import time
from collections import Counter #I read about this package at https://realpython.com/python-histograms/

def main():
	positive_f = 'labeled_nodes.txt'
	All_f = 'interactome-flybase-collapsed-weighted.txt' 
	Outputf = "Clustering_Candidates.txt"
	NeighborOutputf = "Common_Neighbor_Candidates.txt"
	WeightedOutputf = "Weighted_Neighbor_Candidates.txt"

	Pos_Name_set = PositiveReader(positive_f) #makes set of all known positives
	nodes, edges, simpedgeset = InteractomeReader(All_f) #makes edgeset and nodeset from interactome file

	DDict, NDict = DegreeCalc(nodes,edges) #makes dictionaries of neighbors and degree for each node
	ANDDict = NeighborD(nodes, NDict, DDict) #dictionary of the average neighbor degree of each node.
	kANDDict = KeyAnd(ANDDict,DDict)[0] #average AND for all the nodes with degree k
	knodeDict = KeyAnd(ANDDict, DDict)[1] #dictionary of degree->nodes with that degree

	#This section selects candidates based on the number of neighbors they share with known positives
	#It does this either by just counting the number of neighbors OR by weighting neighbors by the # of times
	#they are connected to known positives
	known_neighborset, ranked_nDict = Pos_Neighbor(Pos_Name_set, NDict) #returns set of all neighbors of positives + dictionary of neighbor -> how many positives the neighbor is attached to
	weighted_candidates, weighted_weightD, simple_candidates, simple_weightD = NeighborC_Picker(nodes, NDict, known_neighborset, ranked_nDict, Pos_Name_set)

	''' #This section produces candidates based on clustering coefficient
	startclust = time.time()
	kClustDict,k_list, AvgCv_list, ClustDict = ClusterCo(nodes, edges, NDict, DDict, knodeDict, simpedgeset)
	endclust = time.time()
	changeclust = endclust - startclust
	print("Cluster Time =" + str(changeclust))
	rangemin, rangemax, mean = ClustOutput(ClustDict, Pos_Name_set)
	candidates, candidate_weights = node_picker(nodes, Pos_Name_set, rangemin, rangemax, ClustDict, mean)
	'''

	output_maker(simple_candidates, simple_weightD, NeighborOutputf)
	output_maker(weighted_candidates, weighted_weightD, WeightedOutputf)
	return

def PositiveReader(positive_f):
	Pos_Name_set = set()
	f= open(positive_f,'r') #opens the graph file
	print("Opening " + positive_f)

	for line in f:
		l = str(line)
		l_sublist = l.split()
		if l_sublist[2] == 'Positive':
			Pos_Name_set.add(l_sublist[0])
	return Pos_Name_set

def InteractomeReader(All_f):
	nodeset = set()
	edgeset = set()
	simpedgeset = set()
	f= open(All_f,'r') #opens the graph file
	print("Opening " + All_f)

	for line in f:
		l = str(line)
		l_sublist = l.split()
		node1 = l_sublist[0]
		node2 = l_sublist[1]
		weight = l_sublist[2]
		edge = (node1, node2, weight) #a tuple
		simpedge = (node1, node2)

		if node1 not in nodeset:
			nodeset.add(node1)
		if node2 not in nodeset:
			nodeset.add(node2)

		if edge not in edgeset and node1 != '#symbol1':
			edgeset.add(edge)

		if simpedge not in simpedgeset and node1 != '#symbol1':
			simpedgeset.add(simpedge)
	
	nodeset.discard('#symbol1')
	nodeset.discard('symbol2')

	return nodeset, edgeset, simpedgeset


'''
MOSTLY COPIED FROM HW2
This function takes a node set and edge set as inputs and outputs three things for each interactome:
DDict, a dictionary of the degree of each node, NDict, a dictionary of the neighbors of each node 
(a dictionary of sets)
'''
def DegreeCalc(nodes,edges):
	DDict = {} #dictionary of degrees
	NDict = {} #dictionary of neighbor lists

	for n in nodes:
		DDict[n] = 0 #initializes the degree of each node at 0
		NDict[n] = set() #initializes an empty set for each node
	for e in edges:
		Nset0 = NDict[e[0]]
		Nset1 = NDict[e[1]]
		Nset0.add(e[1]) #adds the second node to the neighbor list of the first node
		Nset1.add(e[0]) #adds the first node to the neighborlist of the second node
		NDict[e[0]] = Nset0 #changes the dictionary content for nodes
		NDict[e[1]] = Nset1
		for i in e[0:2]: 
			newn = int(DDict[i]) + 1 #adds 1 to the previous degree of the node (initialized at 0)
			DDict[i] = newn #changes DDictionary degree of the node
	return DDict, NDict

'''
COPIED FROM HW2
This file takes a node set, the neighbor dictionary (NDict) and the degree dictionary (DDict)
as input and outputs DnDict, a dictionary of the average neighbor degree of each node.
'''
def NeighborD(nodes, NDict, DDict):
	ANDDict = {} #dictionary of neighbor degrees of all the nodes
	for n in nodes: #for every node
		nset = NDict[n]
		dnsum = 0 #sum of the degrees of the neighbors of a single node
		for i in nset: #for every neighbor of the node
			dnsum = dnsum + int(DDict[i]) #adding the degree of each neighbor to the sum
		dnavg = dnsum/len(nset) #the avg neighbor degree of the node n
		ANDDict[n] = dnavg
	return ANDDict

'''
COPIED FROM HW2
This function takes ANDDict (Dictionary node->AND) and DDict(Dictionary node->degree(k)) and outputs kANDDict (k->AND)
which has the average AND for all the nodes with degree k
'''
def KeyAnd(ANDDict, DDict):
	kANDDict = {} #dictionary of the average AND for each k (k->average AND)
	knodeDict = {} #dictionary of a list of nodes of every degree (k->nodelist)
	DictInput = []
	listofk = [] #a list of all the degrees

	#this part is just making a list of all the degrees (each degree will only appear once)
	for node in DDict:
		k = DDict[node]
		if k not in listofk:
			listofk = listofk + [k]

	#Runs through all of the degrees, k, in the degree list and initializes 2 dictionaries
	for k0 in listofk:
		knodeDict[k0] = []#initializing a dictionary that will have: degree->nodes with that degree
		kANDDict[k0] = 0 #initializing a dictionary that will have: degree(k)->average AND for nodes of degree k

	#Runs through all the nodes, finds their degree and adds them to the proper dictionaries.
	for node in DDict:
		k1 = DDict[node] #finding the degree of the node
		kinput = knodeDict[k1] + [node] #adding it to a list of nodes with that degree
		knodeDict[k1] = kinput #updating the list of nodes with degree k1 in the dictionary of degree->nodes with degree

	#goes through all the degrees and finds the average AND of each and adds that value to kANDDict
	for k2 in knodeDict:
		kAND = 0 #average AND for k (reinitialized for each degree)
		for node in knodeDict[k2]: #goes through all the nodes that have a certain degree, k2.
			kAND = kAND + ANDDict[node] #initial sum + the AND of each node in the list of nodes 
		kANDaverage = kAND/int(len(knodeDict[k2])) #sum of the AND for all nodes with degree=k2/#nodes with degree=k2
		kANDDict[k2] = kANDaverage #adds the average AND of nodes with degree=k2 to dictionary.
	return kANDDict, knodeDict

def Pos_Neighbor(Pos_Name_set, NDict):
	known_neighborset = set()
	ranked_nDict = {}
	seenneighbors = set()

	for positive in Pos_Name_set:
		neighbors = NDict[positive] #all of the neighbors of a single positive node
		for n in neighbors: #goes through all of positive's neighbors
			known_neighborset.add(n) #adds all of the neighbors of positive to the list of known neighbors
			if n in seenneighbors: #if it's already been seen
				ranked_nDict[n] = ranked_nDict[n] + 1 #take the previous score +1 
			else:
				ranked_nDict[n] = 1 #if it just got added to the dictionary, value = 1
			seenneighbors.add(n)

	return known_neighborset, ranked_nDict

#goes through all the neighbors of the nodes and ranks weights by either the number of neighbors they share with known positives (simple)
#OR by the cumulative sum of the #positives each neighbor is attached to
#Outputs: two candidate lists,simple_candidates and weighted_candidates, and their respective node weight dictionaries
def NeighborC_Picker(nodes, NDict, known_neighborset, ranked_nDict, Pos_Name_set):
	simple_list = []
	weighted_list = []
	simple_candidates = []
	simple_weightD = {}
	weighted_candidates = []
	weighted_weightD = {}

	sweightlist = [] #just a list of the weights unordered
	#this part goes through simply how many neighbors every node has 
	for node in nodes:
		simp_neighbor_sum = 0 #how many neighbors node has that are neighbors of known positives
		neighbors = NDict[node]
		for n in neighbors: #goes through all the neighbors of the node
			if n in known_neighborset: #if a neighbor is a known neighbor
				simp_neighbor_sum = simp_neighbor_sum + 1 #add 1 to previous number of known neighbors
		if node not in Pos_Name_set: #only nodes that aren't known positives included
			sweightlist.append(simp_neighbor_sum) #this is important to find the max to normalize
			simple_list.append((simp_neighbor_sum, node)) #adding (#neighbors, node) to list
	
	simple_list.sort(reverse = True)

	for x in simple_list[0:100]: #the top 100 candidates based on weight
		simple_candidates.append(x[1]) #adds node to list of top 100 candidates
		simple_weightD[x[1]] = x[0]/max(sweightlist) #adds normalized weight of the node

	wweightlist = [] #just a list of weighted sums
	for node in nodes:
		weightednsum = 0 #sum of weightos of neighbors of the node has that are neighbors of known positives
		neighbors = NDict[node]
		for n in neighbors:
			if n in known_neighborset:
				weightednsum = weightednsum + ranked_nDict[n] #previous sum + weight(#positives it's attached to) of neighbor
		if node not in Pos_Name_set:
			wweightlist.append(weightednsum)
			weighted_list.append((weightednsum, node)) #adding (#neighbors, node) to list

	weighted_list.sort(reverse = True)

	for x in weighted_list[0:100]: #the top 100 candidates based on weight
		weighted_candidates.append(x[1]) #adds node to list of top 100 candidates
		weighted_weightD[x[1]] = x[0]/max(wweightlist) #adds normalized weight of the node

	print(weighted_weightD, simple_weightD)
	return weighted_candidates, weighted_weightD, simple_candidates, simple_weightD


'''
Code copied from Lab 4
Inputs = nodelist, edgelist, adj_list (a dictionary node->neighbors) and node_s, a starting node.
Calculates the shortest path from node_s to every other node.
Output = D, dictionary of pathlength of every node from node_s
'''
def shortest_paths(node_list, edge_list, adj_list, node_s):
	D={} #a dictionary of shortest paths from a node,s, to every other node.
	for node in node_list:
		D[node] = 100 #sets the layer to 1000
	D[node_s] = 0 #sets the layer of node_s to 0 since it is the starting node
	Q = [node_s] #the queue starts off with node_s
	while len(Q)>0: #if there are items in Q
		w = Q[0] #the node to be explored is the first item on the list
		Q.remove(w) #first item on the list is removed
		for neighbor in adj_list[w]:
			if D[neighbor] == 100:
				D[neighbor] = D[w] + 1 #layer set to whatever the previous layer was + 1
				Q.append(neighbor)#adds the neighbor to the end of Q
	return D


'''
COPIED FROM HW2
This is a function to calculate the clustering coefficient
Inputs = nodes (a list), NDict (a dictionary of node->set of neighbors,), DDict (dictionary: node->degree), and knodDict (dictionary k->nodes w/ degree k)
Outputs kClustDict (dictionary of k->avg Cv), k_list (list of all the k in a dataset) and AvgCv_list (list of all the average Cvs in the same orders as the ks)
'''
def ClusterCo(nodes, edges, NDict, DDict, knodeDict, simpedge):
	ClustDict = {} #initializing a dictionary node->clustering coefficient
	kClustDict = {} #initializng a dictionary k-> average clustering coefficient
	k_list = []
	AvgCv_list = []

	#Runs through all the nodes and constructs ClustDict, a node -> Cv dictionary.
	for node in nodes:
		neighbors = NDict[node] # a list of all the neighbors of a node
		seen = set() #just checking I'm not counting edges between neighbors multiple times for one node
		edgesum = 0 #sum of edges between neighbors reinitialized for each node
		for u in neighbors: #goes through each neighbor in the list of neighbors
			for w in neighbors: #and checks if its connected to every other neighbor.
				if (u,w) in simpedge and (u,w) not in seen:
					edgesum = edgesum + 1 #adds 1 for every edge between neighbors
					seen.add((u,w))#to make sure that edges aren't counted multiple times
		denom = (DDict[node] * (DDict[node]-1))
		if denom == 0: #this is just to prevent a divide by 0 error. 
			Cv = 0 #technically it might be more correct for Cv to equal 1 (because if denom = 0 then there is only 1 neighbor and by default 0 edges and 0/0 =1).
			#  However I feel that that is just conceptually bad even if it's technically correct because an interactome with many single nodes would appear to be highly clustered.
		else:
			Cv = edgesum/denom
		ClustDict[node] = Cv #Adding Cv to node->Cv dictionary

	#Runs through all the degrees and constructs kClustDict, a dictionary of k-> average Cv.
	for k in knodeDict: #goes through all the nodes
		CvSum = 0
		for node in knodeDict[k]: #all the nodes that have degree k
			CvSum = CvSum + ClustDict[node] #adds the Cv of every node with degree k
		CvAvg = CvSum/len(knodeDict[k]) #the sum of all the Cvs/number of nodes with degree k
		kClustDict[k] = CvAvg #adding to k->avgCv Dictionary

	#Makes lists of k and avg Cv (same order) for use to make a plot
	for k2 in kClustDict:
		k_list.append(k2) 
	k_list.sort() #this makes the plot look nicer
	for k in k_list: #this is separated so that it will have the same order as the sorted k_list
		AvgCv_list.append(kClustDict[k])

	return kClustDict,k_list, AvgCv_list, ClustDict

def ClustOutput (ClustDict, Pos_Name_set):
	clustlist = []
	for name in Pos_Name_set:
		Co = ClustDict[name]
		clustlist.append(Co)
		
	clustlist.sort()

	minimum = min(clustlist)
	maximum = max(clustlist)
	print("minimum clustering coefficient: ", minimum, "maximum clustering coefficient: ", maximum)
	std = statistics.stdev(clustlist)
	mean = statistics.mean(clustlist)
	print("standard deviation: ", std, "mean: ", mean)
	rangemax = float(mean) + (float(std)/2)
	rangemin = float(mean) - (float(std)/2)
	print(rangemin,rangemax)
	return rangemin, rangemax, mean

#goes through all the nodes and returns a list of nodes that have a clustering coefficient that is closest to the mean of the known positives
def node_picker(nodeset, Pos_Name_set, rangemin, rangemax, ClustDict, mean):
	candidates = []
	candidate_weights = {}
	rawweightlist = []
	tuplelist = []

	for node in nodeset:
		if node not in Pos_Name_set:
			Co = ClustDict[node]
			diff = mean - Co
			AbsoluteDiff = abs(diff)
			tuplelist.append((AbsoluteDiff,node))
	tuplelist.sort()
	for tup in tuplelist[0:100]:
		node = tup[1]
		rweight = tup[0]
		rawweightlist.append(rweight)
		candidates.append(node)
	
	for tup2 in tuplelist[0:100]:
		node = tup2[1]
		rawweight = tup2[0]
		normalizedweight = 1 - (rawweight/max(rawweightlist))
		candidate_weights[node] = normalizedweight
	
	return candidates, candidate_weights

#Just makes an output file with the top candidates
def output_maker(candidates, candidate_weights, outputf):
	f= open(outputf,'w')
	print("Opening " + outputf)

	for candidate in candidates:
		f.write(candidate + '\t' + str(candidate_weights[candidate]) + '\n')

	print("printed to Output file. Done")

main()
