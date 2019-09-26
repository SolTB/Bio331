
'''
Bio 331: HW2
Author: Sol Taylor-Brill

This code does several things: it takes edges files and extracts edges and nodes, calculates degree of each node,
calculates AND of each node, calculates average AND for all the nodes with the same degree, calculates clustering coefficient, 
calculates pathlength and outputs several histograms and plots showing these statistics for each dataset.

Inputs = filenames and the names of the datasets

Parts of the code are copied from the HW2 Handout(noted in comments), or copied from Lab2
Additional information about formatting histograms was found at: https://matplotlib.org/3.1.1/gallery/statistics/histogram_multihist.html
Information about using counter to count the number of times a item appears in a list was found at: https://realpython.com/python-histograms/
'''
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import time
from collections import Counter #I read about this package at https://realpython.com/python-histograms/

def main():
	names = ['Yeast-APMS','Yeast-LC','Yeast-Y2H','Fly','HIPPIE'] #This is from the HW2 Handout
	files = ['Yeast_Combined_APMS.txt','Yeast_LC_Multiple.txt','Yeast_Y2H_Union.txt','Fly_Unpublished.txt',
	'HIPPIE_Unweighted.txt']
	master_node_list = []
	master_DDict_list = []
	kAND_list = []
	PathDict_list = []
	kNodeList = []
	kList_List= []
	AvgCvList_List = []

	#Runs through all five of the datasets named in the files list
	for i in range(len(names)): #from HW2 Handout
		print('DATASET:',names[i])
		print('READING FILE:',files[i])

		#The following section just goes through the file and gets all the information about nodes, edges, degrees and neighbors out
		output = read_edge_file(files[i], names[i])
		edges = output[0]
		nodes = output[1]
		DDict = output[2]
		NDict = output[3]

		#the following lists of lists store the information about each of the datasets so they can be run through later
		master_node_list = master_node_list + [nodes]
		master_DDict_list = master_DDict_list + [DDict]
		ANDDict = NeighborD(nodes, NDict, DDict)
		kANDDict = KeyAnd(ANDDict,DDict)[0]
		knodeDict = KeyAnd(ANDDict, DDict)[1]
		kAND_list.append(kANDDict) 
		kNodeList.append(knodeDict)
		startclust = time.time()
		ClusterOutput = ClusterCo(nodes, edges, NDict, DDict, knodeDict)
		endclust = time.time()
		changeclust = endclust - startclust
		print("Cluster Time =" + str(changeclust))
		k_list = ClusterOutput[1]
		kList_List.append(k_list)#adds the list of k to the end of the list of list for each dataset
		AvgCv_list = ClusterOutput[2]
		AvgCvList_List.append(AvgCv_list)
		AvgCvList_List #appends AvgCv_list to a list containing lists of all five datasets.


		'''
		#this part calculates the first 100,000 shortest paths for the Fly/HIPPIE databases and all the yeast paths
		#I've commented it out so that I can work on the rest of my code without it running unnecessarily.
		start1 = time.time() #Just so I know how long it takes to run
		count = 0
		pathDict = {} #dictionary of pathlengths that gets reinitialized for every dataset
		for u in nodes: #goes through all the nodes to calculate shortest path to every other node from u
			if names[i] == "Fly" or names[i]== 'HIPPIE':
				if (count*(len(nodes)-1)) > 100000: #not calculating if over 100,000
					print("BROKE AT u=", str(count))
					break 
				else:
					D = shortest_paths(nodes,edges, NDict, u)
					pathDict[u] = D #adding to a dictionary of dictionaries
			else:
				D = shortest_paths(nodes,edges, NDict, u)
				pathDict[u] = D
			count = count + 1
		PathDict_list.append(pathDict) #appending a list of dictionaries of dictionaries
		end1 = time.time()
		change1 = end1 - start1
		print("PATHLENGTH TIME:", change1)
		'''

	'''
	This section runs functions that go through the big lists that have information for all five datasets
	Most of the functions are commented out because I already used them to get what I needed. 
	'''
	plot_k(kList_List, AvgCvList_List, names)
	#pathlists_list = HistoFixo(PathDict_list,names)
	#HistoInput(pathlists_list, names) #this function was used to construct histograms of shortest pathlengths
	#SimplePlot(kAND_list, names) #was used to make plot of average AND
	#listofdlists = HistoInput(master_node_list, master_DDict_list, names) #was used to make histograms
	#plot_some_numbers(listofdlists, names) #was used to make plot of log(degree) vs. log(#nodes)
	#ClusterCo(nodes, NDict, DDict)
	return

'''
This is based on plot_some_numbers which was based on the code given by Anna
This function just plots the average AND of each degree for each dataset
inputs = kAND_list, a list of dictionaries of degree to average AND for each data set, and names, a list of names of each dataset.
'''
def SimplePlot(kAND_list,names):
	node_style_list = ['*r-', 'og-', '.b-', '+y-', 'dm-']
	fig = plt.figure(figsize=(6.5,4)) # make a 6.5" wide by 4" tall figure.
	for i in range(len(kAND_list)): #for each dictionary in the list of dictionaries
		klist = [] #reinitiated for each dataset
		ANDlist = [] #reinitiated for each dataset
		for k in kAND_list[i]: #for each degree in the dictionary
			klist.append(k) #the list of degrees + the degree
		klist.sort() #sorts the list so it will look nice on the plot
		for k in klist: #goes through the list of degrees and constructs a list of corresponding average ANDs
			ANDlist.append(kAND_list[i][k]) #appends the average AND of the node
		plt.plot(klist,ANDlist,node_style_list[i], label=names[i])
	plt.xlim(0,100)
	plt.legend(loc='upper right')
	plt.xlabel('k')
	plt.ylabel('average AND')
	plt.title('Degree (k) vs. avg. AND for k-degree nodes')
	plt.tight_layout()
	plt.savefig('kAND1.png')
	return

'''
Plots the log of degree vs. the log of the number of nodes w/ that degree for each dataset onto the same plot
inputs = listofdlists = a list of lists of degrees for each dataset and names, the name of each dataset
This is based on the code that was given to us by Anna Ritz in the HW2 Instructions
'''
def plot_some_numbers(listofdlists, names): #From HW2 instructions
	node_style_list = ['*r-', 'og-', '.b-', '+y-', 'dm-']
	fig = plt.figure(figsize=(6.5,4)) # make a 6.5" wide by 4" tall figure.
	for i in range(len(listofdlists)):
		listofdlists[i].sort()
		counted = Counter(listofdlists[i]) #dictionary of node degree to number of nodes
		listx = []
		listy = []
		for n in counted:
			if counted[n] != 0:
				lognumbernodes = math.log(counted[n])
				logneighbors = math.log(n)
				listx.append(logneighbors)
				listy.append(lognumbernodes)
		#print(listx,listy)
		plt.plot(listx,listy,node_style_list[i], label=names[i])
	plt.legend(loc='upper right')
	plt.xlabel('log(k)')
	plt.ylabel('log(#nodes)')
	plt.title('log(degree) vs. log(#nodes)')
	plt.tight_layout()
	plt.savefig('number1.png')
	print('wrote to numbers.png')
	return

'''
This is basically copied from plot_some_numbers() and then slightly modified so it would work better with the k/Clustering coefficient data.
Inputs: kList_list (list of lists of degrees for all five datasets), AvgCvList_list (list of lists of avgCvs (in the same order as the klists)), and names
Outputs: png file that plots k vs. avg Cv for all five datasets
'''
def plot_k(kList_List, AvgCvList_List, names): #From HW2 instructions
	node_style_list = ['*r-', 'og-', '.b-', '+y-', 'dm-']
	fig = plt.figure(figsize=(6.5,4)) # make a 6.5" wide by 4" tall figure
	for i in range(len(kList_List)):
		plt.plot(kList_List[i], AvgCvList_List[i],node_style_list[i], label=names[i])
	plt.xlim(0,100)
	plt.legend(loc='upper right')
	plt.xlabel('k')
	plt.ylabel('Avg C(v)')
	plt.title('Degree vs. Avg Clustering Coefficient')
	plt.tight_layout()
	plt.savefig('Cluster.png')
	print('wrote to Cluster.png')
	return
'''
Inputs = file name and the dataset name
Reads the file and returns edges, a list of edges, nodes, a list of nodes, DDict, a dictionary of node->degree,
and NDict, a dictionary of node->node's neighbors
'''
def read_edge_file(f,name): #copied partially from lab2
	EdgeList = open(f,'r')
	print("Opened ", f)
	nodes = set()
	edges = set()

	#The following for loop goes through all the lines in the file and adds nodes and edges to their sets
	for l in EdgeList:
		line = str(l)
		nodesinline = line.split()
		nodesinline = tuple(nodesinline)
		if nodesinline not in edges:
			edges.add(nodesinline)
		for n in nodesinline:
			if n not in nodes:
				nodes.add(n)

	#The following section takes the information that was taken from the file and plugs them into the DegreeCalc()
	#function to get the degrees and neighbors of each node.
	startD = time.time()
	DegreeCOutput = DegreeCalc(nodes,edges)
	endD = time.time()
	changeD = endD - startD
	print("DEGREECALC TIME:", changeD)
	DDict = DegreeCOutput[0]
	NDict = DegreeCOutput[1]
	davg = DegreeCOutput[2]
	return edges, nodes, DDict, NDict

'''
This function takes a node set and edge set as inputs and outputs three things for each interactome:
DDict, a dictionary of the degree of each node, NDict, a dictionary of the neighbors of each node 
(a dictionary of sets), and davg, the avg node degree of the interactome
'''
def DegreeCalc(nodes,edges):
	DDict = {} #dictionary of degrees
	NDict = {} #dictionary of neighbor lists
	dsum = 0 #sum of the degree of every node
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
		for i in e: 
			newn = int(DDict[i]) + 1 #adds 1 to the previous degree of the node (initialized at 0)
			DDict[i] = newn #changes DDictionary degree of the node
	for node in DDict:
		dsum = dsum + DDict[node]
	davg = dsum/len(nodes)
	return DDict, NDict, davg

'''
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

'''
This is a function to calculate the clustering coefficient
Inputs = nodes (a list), NDict (a dictionary of node->set of neighbors,), DDict (dictionary: node->degree), and knodDict (dictionary k->nodes w/ degree k)
Outputs kClustDict (dictionary of k->avg Cv), k_list (list of all the k in a dataset) and AvgCv_list (list of all the average Cvs in the same orders as the ks)
'''
def ClusterCo(nodes, edges, NDict, DDict, knodeDict):
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
				if (u,w) in edges and (u,w) not in seen:
					edgesum = edgesum + 1 #adds 1 for every edge between neighbors
					seen.add((u,w))#to make sure that edges aren't counted multiple times
		denom = (DDict[node] * (DDict[node]-1))
		if denom == 0:
			Cv = 0
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

	return kClustDict,k_list, AvgCv_list

'''
Code copied from Lab 4
Inputs = nodelist, edgelist, adj_list (a dictionary node->neighbors) and node_s, a starting node.
Calculates the shortest path from node_s to every other node.
Output = D, dictionary of pathlength of every node from node_s
'''
def shortest_paths(node_list, edge_list, adj_list, node_s):
	D={}
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
Inputs = nodelist (a list of lists of nodes) and DDict (a list of degree dictionaries))
Outputs a figure with six histograms: 1 for each dataset + 1 with all 5 overlaid
returns listofdlists a list of lists of degrees (for use later on)
'''
def HistoInput(listoflists,names):
	facecolor_list = ['blue', 'red', 'cyan', 'magenta', 'green']
	fig, axes = plt.subplots(nrows=2, ncols=3) #I got information for how to make six plots in a grid from https://matplotlib.org/3.1.1/gallery/statistics/histogram_multihist.html
	for j in range(len(listoflists)): #this just makes a separate histogram for each dataset
		num = j + 1
		plt.subplot(2,3,num) #I copied this general plot format from the plot_some_numbers() funtion that was provided in the HW2 Handout
		plt.hist(listoflists[j],bins=25, range= (0,100), color=facecolor_list[j])
		plt.xlabel('Pathlength')
		plt.ylabel('#paths')
		plt.title(names[j])
	plt.subplot(2,3,6) #this is the extra combined histogram
	plt.hist(listoflists,bins=25,range = (0,100), color=facecolor_list)
	plt.xlabel('Pathlength')
	plt.ylabel('#paths')
	plt.title('All')
	print("Finished")
	Histoname = "Histogram" + "Revised9" + '.png'
	plt.xlim(0,100)
	fig.tight_layout()
	plt.savefig(Histoname)
	print('wrote to ', Histoname)
	return 

def DList_List_Maker(nodelist, DDict):
	for n in nodelist[i]: #for all the nodes in a singl dataset
		value = DDict[i][n] #the degree of that node
		Dlist.append(value) #adding the degree to the list
	listofdlists = listofdlists + [Dlist] #adding list of degrees to a list of lists

'''
Input = PathDict_list (a list of dictionaries of dictionaries) and names (a list of dataset names)
Output = a list of lists of shortest paths for all five data sets to use to make the histograms.
'''
def HistoFixo(PathDict_list, names):
	pathlists_list = [] #a list of lists of the first 100,000 pathlengths for each dataset
	for i in range(len(PathDict_list)): #should go through the dictionary of dictionaries for all five datasets
		pathlist = [] # a list of the first ~100,000 pathlengths for dataset i
		for u in PathDict_list[i]: #goes through all the sub_dictionaries (shortest paths from each node)
			for w in PathDict_list[i][u]: #for each neighbor node to node u in dataset i
				if names[i] == 'Fly' or names[i] == 'HIPPIE':
					if len(pathlist)>100000:
						break
					else:
						pathlist.append(PathDict_list[i][u][w]) #append path length to w from u in dataset i
				else:
					pathlist.append(PathDict_list[i][u][w])
		pathlists_list.append(pathlist)
	return pathlists_list


main()