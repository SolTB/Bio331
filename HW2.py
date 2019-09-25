import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import time
from collections import Counter #I read about this package at https://realpython.com/python-histograms/

'''
This code does several things: it takes edges files and extracts edges and nodes, calculates degree of each node,
calculates AND of each node, calculates average AND for all the nodes with the same degree, calculates clustering coefficient, 
calculates pathlength and outputs several histograms and plots showing these statistics for each dataset.

Inputs = filenames and the names of the datasets

Parts of the code are copied from the HW2 Handout(noted in comments), or copied from Lab2
Additional information about formatting histograms was found at: https://matplotlib.org/3.1.1/gallery/statistics/histogram_multihist.html
Information about using counter to count the number of times a item appears in a list was found at: https://realpython.com/python-histograms/
'''

def main():
	#plot_some_numbers()
	#filename = "example.txt"
	names = ['Yeast-APMS','Yeast-LC','Yeast-Y2H','Fly','HIPPIE'] #This is from the HW2 Handout
	files = ['Yeast_Combined_APMS.txt','Yeast_LC_Multiple.txt','Yeast_Y2H_Union.txt','Fly_Unpublished.txt',
	'HIPPIE_Unweighted.txt']
	master_node_list = []
	master_DDict_list = []
	kAND_list = []
	for i in range(len(names)): #from HW2 Handout
		print('DATASET:',names[i])
		print('READING FILE:',files[i])
		output = read_edge_file(files[i], names[i])
		nodes = output[1]
		DDict = output[2]
		NDict = output[3]
		master_node_list = master_node_list + [nodes]
		master_DDict_list = master_DDict_list + [DDict]
		ANDDict = NeighborD(nodes, NDict, DDict)
		kANDDict = KeyAnd(ANDDict,DDict)
		kAND_list = kAND_list + [kANDDict]
	#SimplePlot(kAND_list, names) #was used to make plot of average AND
	#listofdlists = HistoInput(master_node_list, master_DDict_list, names) #was used to make histograms
	#plot_some_numbers(listofdlists, names) #was used to make plot of log(degree) vs. log(#nodes)
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
	plt.xlim(0,1000)
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
Inputs = file name and the dataset name
Reads the file and returns edges, a list of edges, nodes, a list of nodes, DDict, a dictionary of node->degree,
and NDict, a dictionary of node->node's neighbors
'''
def read_edge_file(f,name): #copied partially from lab2
	EdgeList = open(f,'r')
	print("Opened ", f)
	nodes = set()
	edges = set()
	for l in EdgeList:
		line = str(l)
		nodesinline = line.split()
		nodesinline = tuple(nodesinline)
		if nodesinline not in edges:
			edges.add(nodesinline)
		for n in nodesinline:
			if n not in nodes:
				nodes.add(n)
	startD = time.time()
	DegreeCOutput = DegreeCalc(nodes,edges)
	endD = time.time()
	changeD = endD - startD
	print("DEGREECALC TIME:", changeD)
	DDict = DegreeCOutput[0]
	NDict = DegreeCOutput[1]
	davg = DegreeCOutput[2]
	DnDict = NeighborD(nodes, NDict, DDict)
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
	for node in DDict:
		k = DDict[node]
		if k not in listofk:
			listofk = listofk + [k]
	for k0 in listofk:
		knodeDict[k0] = []
		kANDDict[k0] = 0
	for node in DDict:
		k1 = DDict[node]
		kinput = knodeDict[k1] + [node]
		knodeDict[k1] = kinput
	for k2 in knodeDict:
		kAND = 0 #average AND for k (reinitialized for each degree)
		for node in knodeDict[k2]:
			kAND = kAND + ANDDict[node] #initial sum + the AND of each node in the list of nodes 
		kANDaverage = kAND/int(len(knodeDict[k2]))
		kANDDict[k2] = kANDaverage
	return kANDDict

'''
Outputs a figure with six histograms: 1 for each dataset + 1 with all 5 overlaid
inputs = nodelist (a list of lists of nodes) and DDict (a list of degree dictionaries))
returns listofdlists a list of lists of degrees (for use later on)
'''
def HistoInput(nodelist,DDict,names):
	facecolor_list = ['blue', 'red', 'cyan', 'magenta', 'green']
	fig, axes = plt.subplots(nrows=2, ncols=3) #I got information for how to make six plots in a grid from https://matplotlib.org/3.1.1/gallery/statistics/histogram_multihist.html
	listofdlists = []
	for i in range(len(nodelist)): #going through the nodes from each dataset
		Dlist = []
		color = facecolor_list[i]
		for n in nodelist[i]: #for all the nodes in a singl dataset
			value = DDict[i][n] #the degree of that node
			Dlist.append(value) #adding the degree to the list
		listofdlists = listofdlists + [Dlist] #adding list of degrees to a list of lists
	for j in range(len(nodelist)): #this just makes a separate histogram for each dataset
		num = j + 1
		plt.subplot(2,3,num) #I copied this general plot format from the plot_some_numbers() funtion that was provided in the HW2 Handout
		plt.hist(listofdlists[j],bins='auto',color=facecolor_list[j])
		plt.xlabel('k')
		plt.ylabel('#nodes')
		plt.title(names[j])
	plt.subplot(2,3,6) #this is the extra combined histogram
	plt.hist(listofdlists,bins='auto',color=facecolor_list)
	plt.xlabel('k')
	plt.ylabel('#nodes')
	plt.title('All')
	print("Finished")
	Histoname = "Histogram" + "Revised3" + '.png'
	plt.xlim(0,100)
	fig.tight_layout()
	plt.savefig(Histoname)
	print('wrote to ', Histoname)
	return listofdlists


main()