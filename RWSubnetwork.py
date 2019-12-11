'''
This program runs a random walk simulation on the subnetwork of candidates using either each node as a starting node
OR using known positives as starting nodes (potentially better ways to do this).
Author: Sol Taylor-Brill
Last Edit: 12/10/19
'''

import math
import random
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session

def main():
	## READ INTERACTOME AND KNOWN POSITIVE FILES ##
	positive_f = 'labeled_nodes.txt' #known-positives
	interactome_f = 'interactome-flybase-collapsed-weighted.txt'#whole interactome

	Pos_Name_set = PositiveReader(positive_f) #set of known positives
	nodeset, edgeset, simpedgeset = InteractomeReader(interactome_f) #set of ALL nodes, ALL edges + weights, ALL edges

	## CANDIDATE FILES ##
	Maddy_f = "Maddy_candidates_10.txt" #unweighted, 10
	CN_f = "Common_Neighbor_Candidates.txt" #weighted, 100
	WN_f = "Weighted_Neighbor_Candidates.txt" #weighted, 100
	Co_F = "Clustering_Candidates.txt" #weighted, 100

	All_Candidates_List = (Maddy_f, CN_f, WN_f, Co_F) #list of all candidate files to go through
	All_Candidate_Nodes = set() #a (unranked) set of all nodes ranked as top 10 by different methods (NOTE: RANK NOT PRESERVED)

	## PRODUCES SUBNETWORK ##
	#Runs through all candidate files and produces set of all top 10 nodes produced by different methods
	for f in All_Candidates_List:
		candidates, RankDict = weighted_candidate_reader(f)
		for node in candidates: #goes through all 10 nodes
			All_Candidate_Nodes.add(node) #Adds node to set.


	subedge_set, sub_simp_edge_set = get_subnetwork(All_Candidate_Nodes, edgeset) #Finds edges between candidates
	CandidatesWITHPositives = list(All_Candidate_Nodes) + list(Pos_Name_set)
	posedge_set, pos_simp_edge_set = get_subnetwork(CandidatesWITHPositives, edgeset) #Finds edges between candidates + Positives

	print('\n', "#Nodes in subnetwork = ", len(All_Candidate_Nodes))
	print("#Edges in subnetwork = ", len(subedge_set), '\n')

	## CALCULATE NEIGHBORS AND DEGREE ##
	DDict, NDict, TopDList = DegreeCalc(All_Candidate_Nodes, sub_simp_edge_set) #DDict = node-> degree dict; NDict = node->neighbors Dict; TopDList = 10 nodes with highest degree
	DPosDict, NPosDict, TopDList = DegreeCalc(CandidatesWITHPositives, pos_simp_edge_set) #DDict = node-> degree dict; NDict = node->neighbors Dict; TopDList = 10 nodes with highest degree
	
	## RANDOM WALKS ##
	T = 1000 #number of time steps for random walks simulation
	q = 0.95 #teleportation probability for random walks

	MasterTimeDict, Top10RW = run_rw(All_Candidate_Nodes, NDict, T, q)
	PositiveTimeDict, Top10wPos = rw_with_positives(CandidatesWITHPositives, NPosDict, T, q, Pos_Name_set)

	print("Top 10 nodes by random walks with restarts. Starting from every node", Top10RW, '\n')
	print("Top 10 nodes by random walks with restarts. Starting from known positives", Top10wPos)

	## OUTPUTS ##
	#MAKE A FUNCTION TO OUTPUT A RANKED NODE LIST FOR SOME SUBSECTION OF THESE

	#Posting RW file GraphSpace
	#Graph_Maker(All_Candidate_Nodes, sub_simp_edge_set, MasterTimeDict)
	#Graph_Maker(CandidatesWITHPositives, pos_simp_edge_set, PositiveTimeDict)

	return

## READ FILES AND GET SUBNETWORK SECTION ##


#Input: text file of known positives
#Output: Set containing all known positives
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

#Input: text file containing edges with their weight and other info
#Output: nodeset (ALL nodes in the graph), edgeset (ALL edges in the graph (node1, node2, weight))
# and simpedgeset (ALL edges in the graph WITHOUT the weight ie (node1, node2))
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

#returns list of candidates
def weighted_candidate_reader(candidate_f):
	candidates = [] #list of candidates from first to last
	rankDict = {}
	count = 1 #FOR NOW ONLY TAKING 10 or FEWER CANDIDATES

	f= open(candidate_f,'r') #opens the graph file
	print("Opening " + candidate_f)

	for line in f: #goes through every gene in the file
		if count < 11: #Only goes through first 10 lines in the graph (if less than 10 candidates it will take all of them)
			l = str(line) #takes the line and makes it a string
			l_sublist = l.split() #splits items in string by tabs
			candidates.append(l_sublist[0]) #Add gene to list of candidates
			rankDict[l_sublist[0]] = 1 - count/10 #node-> ranking (rank/10)
			count = count + 1
	return candidates, rankDict

#Input: All_Candidate_Nodes (set of all top 10 nodes predicted by different methods) & edgeset (ALL edges in interactome)
#Output: subedge_set (set of all edges connecting any two candidate nodes w/ weights), subsimpedge_set (basically subedge_set w/o weights)
def get_subnetwork(All_Candidate_Nodes, edgeset):
	subedge_set = set() #set of all edges connecting any two nodes in subnetwork (w/ weights)
	subsimpedge_set = set() #set of all edges connecting any two nodes in subnetwork (w/o weights)

	for edge in edgeset:
		if edge[0] in All_Candidate_Nodes and edge[1] in All_Candidate_Nodes: #if both nodes in the edge are candidates
			subedge_set.add(edge) #adds the whole edge to edgeset
			subsimpedge_set.add((edge[0],edge[1])) #adds edge w/o weight to simplified edge set

	return subedge_set, subsimpedge_set


## CALCULATE DEGREE AND NEIGHBORS ##

'''
MOSTLY COPIED FROM HW2
This function takes a node set and edge set as inputs and outputs three things for each interactome:
DDict, a dictionary of the degree of each node, NDict, a dictionary of the neighbors of each node 
(a dictionary of sets), and TopDList (list of top 10 nodes with highest subnetwork degree)
'''
def DegreeCalc(nodes,edges):
	DDict = {} #dictionary of degrees
	DList = [] #(degree,node) list
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

	for x in DDict:
		DList.append((DDict[x],x))
	DList.sort(reverse=True)
	TopDList = DList[0:10] #top 10 nodes ranked by degree

	return DDict, NDict, TopDList


## RANDOM WALKS SIMULATION ##


def rw_with_positives(nodes, NDict, T, q, positives):
	MasterTimeDict = {}
	DictionaryList = []
	ListToSort = []


	#goes through the simulation for each node as the starting node
	for s in positives: #only positives as start nodes.
		TimeDict = rw_simulate(nodes,NDict,s,T,q)
		DictionaryList.append(TimeDict)

	#goes through all the nodes and finds the average from the simulations
	for n in nodes:
		sumtime = 0 #sum from every dictionary
		for Dict in DictionaryList: #goes through all dictionaries
			sumtime = sumtime + Dict[n] #adds dictionary value to sum
		avgsum = sumtime/len(DictionaryList)
		MasterTimeDict[n] = avgsum #node value assigned as the mean from the simulations
		ListToSort.append((avgsum, n))


	ListToSort.sort(reverse=True) #sort node list by random walk probability
	Top10RW = ListToSort[0:10] #top 10 nodes by random walk probability

	return MasterTimeDict, Top10RW

#this function runs random walks with each node taking a turn as the start node and then averages the results
# NOTE: It might be better to just pick a known important node such as NMII and have that be the start node. Ask on Wed.
#Input: nodes (set of nodes in subnetwork); NDict (dictionary of node-> neighbors (NOT DIRECTED)); 
#T (#times to run simulation); & q (probability of the node teleporting to the start)
#Output: Dictionary with averaged times from simulations from each starting node and Top10RW (top 10 nodes by random walk prob)
def run_rw(nodes, NDict, T, q):
	MasterTimeDict = {}
	DictionaryList = []
	ListToSort = []


	#goes through the simulation for each node as the starting node
	for s in nodes:
		TimeDict = rw_simulate(nodes,NDict,s,T,q)
		DictionaryList.append(TimeDict)

	#goes through all the nodes and finds the average from the simulations
	for n in nodes:
		sumtime = 0 #sum from every dictionary
		for Dict in DictionaryList: #goes through all dictionaries
			sumtime = sumtime + Dict[n] #adds dictionary value to sum
		avgsum = sumtime/len(DictionaryList)
		MasterTimeDict[n] = avgsum #node value assigned as the mean from the simulations
		ListToSort.append((avgsum, n))


	ListToSort.sort(reverse=True) #sort node list by random walk probability
	Top10RW = ListToSort[0:10] #top 10 nodes by random walk probability

	return MasterTimeDict, Top10RW


#COPIED FROM HW3
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
	if len(G_out[CurrentNode]) == 0: #if there are no neighbors in the subnetwork (ie you start on a lone node)
		return TimeDict #just return a dictionary full of 0s

	for t in range(1,T):
		start = CurrentNode
		n_out = list(G_out[start]) #a list of all outgoing nodes of the current node
		num = random.random() #generates random number between 0 and 1
		if num > q: 
			CurrentNode = s #redefining the current node. Teleport back to source
		else:
			newnode = random.choice(n_out) #choose random outgoing neighbor as new current node
			CurrentNode = newnode
		numtimes = TimeDict[CurrentNode] + 1 #adds 1 to the node that was just moved to
		TimeDict[CurrentNode] = numtimes #redefines node that just got moved to
	return TimeDict



## OUTPUTS ##


#Makes graph from Random walks output and posts to Graphspace
#Input: nodes in subnetwork, set of edges in subnetwork, TimeDict (master time dictionary)
#Output: Posts graph to graphspace
def Graph_Maker(node_list, edge_list, TimeDict):
	G = GSGraph()
	G.set_name('GroupProject_Test_w/positives')
	G.set_tags(['GroupProject'])

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
	for n in node_list:
		poplabel = "Count: " + str(TimeDict[n]) + ";Normalized: " + str(NewTimeDict[n])
		G.add_node(n,popup=poplabel, label=n)
		blue = NewTimeDict[n] #a value between 0 and 1 (most visited = most blue)
		red = 1-blue #opposite of blue
		nodecolor = rgb_to_hex(red,0,blue)
		G.add_node_style(n, width = 20, height = 20, color=nodecolor,)
	for e in edge_list:
		G.add_edge(e[0],e[1])
		G.add_edge_style(e[0],e[1], width = 2)
	
	graphspace.post_graph(G)
	print("Graph Updated")

#COPIED FROM HW3
#Normalizes the TimeDict values so that they are between 0 and 1
def Normalizer(TimeDict):
	TimeList = [] #a list of the number of time each node has been visited
	NormTimeDict = {}
	for node in TimeDict:
		TimeList.append(TimeDict[node])
	for n in TimeDict:
		NormTimeDict[n] = TimeDict[n]/max(TimeList)
	return NormTimeDict

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

	