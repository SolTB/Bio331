'''
This program will form a subnetwork with only the top 10 nodes that were suggested by the other methods 
(and the edges that connect those nodes) and then will use network analyses (ie Random Walks)
to examine key nodes within the subnetwork. I will also look at various centrality measures 
for the nodes (ie degree of the node in the subnetwork, the clustering coefficient and shortest path betweenness)

QUESTION: SHOULD I INCLUDE THE KNOWN-POSITIVES IN THE SUBNETWORK??
QUESTIONS: DIJKSTRA: Between candidates or between positives?; Also use edge weights or use ranking?; 
RANDOM WALKS: from any candidate, from all positives, or from a specific node (ie NMII)?

Author: Sol Taylor-Brill
Date of last edit: 12/10/19
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
	MasterRankDict = {}
	RankDictList = []



	## MAKES SUBNETWORK ##

	#Runs through all candidate files and produces set of all top 10 nodes produced by different methods
	'''
	for f in All_Candidates_List:
		candidates = unweighted_candidate_reader(f)
		for node in candidates: #goes through all 10 nodes
			All_Candidate_Nodes.add(node) #Adds node to set. RANK NOT PRESERVED

	'''

	#Runs through all candidate files and produces set of all top 10 nodes produced by different methods
	for f in All_Candidates_List:
		candidates, RankDict = weighted_candidate_reader(f)
		for node in candidates: #goes through all 10 nodes
			All_Candidate_Nodes.add(node) #Adds node to set.
			RankDictList.append(RankDict)

	#Creates a master dictionary with an average rank (sum/total methods) for each node
	for node in All_Candidate_Nodes:
		nodesum = 0
		for Dict in RankDictList:
			if node in Dict:
				nodesum = nodesum + Dict[node] #rank already normalized between 0-1
		ranksum = nodesum/len(RankDictList) #sum of the rankings (0-1) divided by total #methods
		MasterRankDict[node] = ranksum

	subedge_set, sub_simp_edge_set = get_subnetwork(All_Candidate_Nodes, edgeset) #Finds edges between candidates
	sub_edge_cost = weights_to_costs(subedge_set) #set of edges with costs instead of weights (for use in dijkstra)
	rank_edge_cost, rank_edge_weight = noderank_to_edgecost(sub_simp_edge_set, MasterRankDict)
	print(rank_edge_weight)

	print('\n', "#Nodes in subnetwork = ", len(All_Candidate_Nodes))
	print("#Edges in subnetwork = ", len(subedge_set), '\n')

	CandidatesWITHPositives = list(All_Candidate_Nodes) + list(Pos_Name_set)
	posedge_set, pos_simp_edge_set = get_subnetwork(CandidatesWITHPositives, edgeset) #Finds edges between candidates + Posi
	



	## ANALYSES ##

	#Simple Analyses (Degree, Clustering, Shortest Paths):
	DDict, NDict, TopDList = DegreeCalc(All_Candidate_Nodes, sub_simp_edge_set) #DDict = node-> degree dict; NDict = node->neighbors Dict; TopDList = 10 nodes with highest degree
	DPosDict, NPosDict, TopDList = DegreeCalc(CandidatesWITHPositives, pos_simp_edge_set) #W/POSITIVES INCLUDED
	ClustDict, Top10Cv = ClusterCo(All_Candidate_Nodes, subedge_set, NDict, DDict, sub_simp_edge_set) #ClustDict (node->Cv Dict); Top10Cv = list of 10 nodes with highest Cv
	allpathslist, st_pathDict = st_paths(All_Candidate_Nodes,sub_edge_cost) #allpathslist (list of shortest paths between all nodes); st_pathDict ((s,t) -> shortest paths dictionary)
	NodePathDict, topShortestPathNodes = in_shortest_paths(allpathslist, All_Candidate_Nodes) #NodePathDict (node->#paths), topShortestPathNodes (10 nodes the appear in the highest number of shortest path)

	print("Top 10 Nodes by Degree: ", TopDList, '\n')
	print("Top 10 Nodes by Clustering Coefficient", Top10Cv, '\n')
	print("Top 10 Nodes by #Shortest Paths", topShortestPathNodes, '\n')


	#Random Walks
	DPosDict, NPosDict, TopDList = DegreeCalc(CandidatesWITHPositives, pos_simp_edge_set) #DDict = node-> degree dict; NDict = node->neighbors Dict; TopDList = 10 nodes with highest degree
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


#COPIED from "SimpleAnalyses.py" (a different file I made)
#Input: txt file containing unweighted ranked list of candidates
#Output: ordered list of candidates
#NOTE: LATER I MIGHT WANT TO MODIFY THIS TO SOMEHOW PRESERVE RANKS IN A DICTIONARY OR SOMETHING
def unweighted_candidate_reader(candidate_f):
	candidates = [] #list of candidates from first to last
	count = 0 #FOR NOW ONLY TAKING 10 or FEWER CANDIDATES

	f= open(candidate_f,'r') #opens the graph file
	print("Opening " + candidate_f)

	for line in f: #goes through every gene in the file
		if count < 10: #Only goes through first 10 lines in the graph (if less than 10 candidates it will take all of them)
			l = str(line) #takes the line and makes it a string
			l_sublist = l.split() #splits items in string by tabs
			candidates.append(l_sublist[0]) #Add gene to list of candidates
			count = count + 1

	return candidates

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

#Input: edgeset of subgraph (node1,node2) and dictionary with averaged rank of each node
#Output: returns edge set with cost (1 - (rank1 + rank2)/2) and edgeset with weight ((rank1 + rank2)/2)
def noderank_to_edgecost(sub_simp_edge_set, MasterRankDict):
	rank_edge_cost = set()
	rank_edge_weight = set()

	for edge in sub_simp_edge_set:
		rank1 = MasterRankDict[edge[0]] #avg rank of first node
		rank2 = MasterRankDict[edge[1]] #avg rank of second node
		edgeweight = (rank1 + rank2)/2
		rank_edge_weight.add((edge[0],edge[1],edgeweight))
		edgecost = 1 - edgeweight
		rank_edge_cost.add((edge[0],edge[1],edgecost))

	return rank_edge_cost, rank_edge_weight

#input: edge set of subnetwork with weights (node1, node2, weight)
#output: edge set of subnetwork with costs (node1, node2, cost (1-weight))
def weights_to_costs(subedge_set):
	sub_edge_cost = set()
	for edge in subedge_set:
		weight = edge[2]
		cost = 1-float(weight)
		sub_edge_cost.add((edge[0],edge[1],cost))
	return sub_edge_cost


## BASIC ANALYSES SECTION ##


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

'''
COPIED FROM HW2 and modified slightly
This is a function to calculate the clustering coefficient
Inputs = nodes (a list), NDict (a dictionary of node->set of neighbors,), DDict (dictionary: node->degree)
Outputs: ClustDict (Node -> Cv dictionary) and Top10Cv (list of 10 nodes with highest Cv)
'''
def ClusterCo(nodes, edges, NDict, DDict, simpedge):
	ClustDict = {} #initializing a dictionary node->clustering coefficient
	NodeCvList = []

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
		NodeCvList.append((Cv,node)) #adds (Cv,node) tuple to list

	NodeCvList.sort(reverse = True) #sorts by the highest Cv
	Top10Cv = NodeCvList[0:10]
	return ClustDict, Top10Cv


## SHORTEST PATHS SECTION ##


#MODIFIED FROM dijkstra() from "hw_functions.py" by Anna Ritz
#Returns ALL shortest paths not just 1
## Run Dijkstra's in the weighted, undirected graph.
## INPUT: set of nodes, 3-element list of edges [node1,node2,weight], source s
## OUTPUT: Dictionary of distances (D), Dictionary of predecessors (pi) 
def dijkstra_ALL(V,E,s):
	adj_list = get_adj_list_with_weights(E) #dictionary of dictionaries

	D = {n:10000 for n in V} #dictionary of distances from n->s. All distances start at 10000 (arbitrary large number)
	D[s] = 0 #distance from s->s is 0 

	pi = {n:None for n in V} #dictionary of predecessors. All nodes start at None

	Q = {n:D[n] for n in V} #dictionary of nodes to distance to that node. Basically a priority queue

	while len(Q) > 0: ## While the queue is not em[ty]
		## Find the node with the minimum weight.
		w = None #"winning" node starts at "None"
		for n in Q: #goes through all nodes in the queue
			if w == None or Q[n] < Q[w]: #if at the very beginning or if n is better than the previous w
				w = n #replace w with n

		del Q[w] #remove w from queue so that the while loop works correctly

		if w in adj_list:
			for x in adj_list[w]: #for all neighbors of w
				newpath = float(D[w]) + float(adj_list[w][x]) #the length of the current shortest path
				if D[x] > newpath:  #if the distance from x is larger than the shortest path so far
					D[x] = newpath # update the distance
					pi[x] = [w] # overwrite the predecessor with LIST containing w
					Q[x] = D[x] #update position in the queue

			#this is where the code diverges from a typical dijkstra
				elif D[x] == newpath: #if the path is just as good as a previous short path (so all shortest paths are found)
					newlist = pi[x] + [w]
					pi[x] = newlist
	return D,pi

##COPIED FROM get_all_paths.py. Written by Anna Ritz.
## Given a predecessor dictionary (e.g, pis from dijkstra()) and 
## a node, returns ALL paths from the node (t) to the source (s) used for 
## the predecessor dictionary. This is a recursive algorithm and expects
## that the predecessor dicationry has ALL potential predecessors for a node.
def get_all_paths(pis,node):
	## while there IS a predecessor...
	while pis[node] != None:  
		all_paths = [] # make a list of lists for current paths.
		## get all paths ending at each predecessor
		for node2 in pis[node]:
			next = get_all_paths(pis,node2)
			## append THIS node to the end of the list
			for i in range(len(next)):
				next[i].append(node)			
			## upate all paths
			all_paths += next
		## return all paths calculated at this level.
		return all_paths

	## If we're at the soure, return the node (base case)
	return [[node]] 

## THIS FUNCTION WAS COPIED FROM "hw_functions.py" by Anna Ritz
## Make an adjacency list that contains the weights of each edge.
## e.g., for edge (u,v), you can access the weight of that edge
## with adj_list[u][v] OR adj_list[v][u]
## INPUT: 3-element list of edges [node1,node2,weight]
## OUTPUT: dictionary of dictionaries
def get_adj_list_with_weights(edges):
	adj_list = {}
	## loop over all edges.
	for u,v,w in edges: ## another way to specify elements of key

		## We want to add the key-value pair (v,w) to adj_list[u].
		## First see if u is a key in adj_list.
		if u not in adj_list:
			adj_list[u] = {}  ## add the key (value is a DICTIONARY)
		## Add the key-value pair (v,w) to adj_list[u]
		adj_list[u][v] = w

		## We want to add the key-value pair (u,w) to adj_list[v].
		## First see if v is a key in adj_list.
		if v not in adj_list:
			adj_list[v] = {}  ## add the key (value is a DICTIONARY)
		## Add the key-value pair (u,w) to adj_list[v]
		adj_list[v][u] = w
		
	return adj_list

#COPIED FROM HW6
#input = V (nodeset) and E (edgelist)
#output = dictionary of all pi dicts for every node in V (node->pi dict with node as s)
def all_pi(V,E):
	ALL_pi = {} #dictionary of all pi dicts from any source node s in V
	for s in V: #goes through all nodes
		D, pi = dijkstra_ALL(V,E,s)
		ALL_pi[s] = pi #adds pi dict to dictionary of pi dicts corresponding to source node
	return ALL_pi

#MODIFIED FROM HW6
#input = V (nodeset) and E (edgelist)
#output = allpathslist (list of nodes in every shortest path between any two nodes in the graph)
#& st_pathDict (dictionary (s,t)-> list of all shortest paths between s and t)
def st_paths(V,E):
	ALL_pi = all_pi(V,E) #dicitionary of all pi dicts for every node in V
	st_pathDict = {} #initializing dictionary of (s,t) -> list of all shortest paths between s and t
	allpathslist = []

	for s in V: #for all source nodes
		for t in V: # for all terminal nodes
			if s!=t: #if the source is NOT the terminal
				allpaths = get_all_paths(ALL_pi[s],t) #for pi dict for source s and terminal t
				paths = Path_Output(allpaths)
				allpathslist.append(allpaths)
				st_pathDict[(s,t)] = paths #(s,t) as a tuple calls all the shortest paths between them
	return allpathslist, st_pathDict

#I just modified this function from HW5 to give me a list of edges in a path rather than an ordered list of nodes
#input: Paths (list of ordered lists of nodes (which are paths))
#output: list of lists of tuples (which are edges)
def Path_Output(Paths):
	P_edgelist = [] #list of edges (tuples) in T
	for sublist in Paths:
		subedgelist = [] #list of edges in 1 path
		for i in range(len(sublist)-1):
			edge = (sublist[i], sublist[i+1])
			subedgelist.append(edge)
		P_edgelist.append(subedgelist)
	return P_edgelist

#Input: allpathslist (list of shortest paths between any two nodes in the graph([u,v,x][u,y,x]])
#Output: NodePathDict (node-> # of paths), topShortestPathNodes (list of nodes that are in the highest number of shortest paths)
def in_shortest_paths(allpathslist, All_Candidate_Nodes):
	NodePathDict = {} #node->#paths dict
	AllShortestPathNodes = [] # [(#paths, node)] list

	#for every candidate node goes through all the lists and counts how many times it appears in a shortest path
	for node in All_Candidate_Nodes:
		pathnum = 0
		for l in allpathslist: #lists in all paths (all of the paths between 2 nodes)
			for sub_l in l: #All sublists in the list: (a path in the list of all paths between 2 nodes)
				if node in sub_l: #if the node is in a sublist (ie a path)
					pathnum = pathnum + 1
		NodePathDict[node] = pathnum
		AllShortestPathNodes.append((pathnum,node))

	AllShortestPathNodes.sort(reverse = True)
	topShortestPathNodes = AllShortestPathNodes[0:10]
	return NodePathDict, topShortestPathNodes



## RANDOM WALK SECTION ##
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