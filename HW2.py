import matplotlib
import matplotlib.pyplot as plt
import math
import time
from collections import Counter #I read about this package at https://realpython.com/python-histograms/

def main():
	#plot_some_numbers()
	#filename = "example.txt"
	names = ['Yeast-APMS','Yeast-LC','Yeast-Y2H','Fly','HIPPIE'] #This is from the HW2 Handout
	files = ['Yeast_Combined_APMS.txt','Yeast_LC_Multiple.txt','Yeast_Y2H_Union.txt','Fly_Unpublished.txt',
	'HIPPIE_Unweighted.txt']
	'''
	for i in range(len(names)): #from HW2 Handout
		print('DATASET:',names[i])
		print('READING FILE:',files[i])
		read_edge_file(files[i], names[i])
	'''
	read_edge_file('Fly_Unpublished.txt', 'example')
	return

def plot_some_numbers(): #From HW2 instructions
	fig = plt.figure(figsize=(6.5,4)) # make a 6.5" wide by 4" tall figure.
	x = list(range(1,10))
	y = [xval/1.5+2 for xval in x]
	logx = [math.log(a) for a in x]
	logy = [math.log(b) for b in y]
	plt.subplot(1,2,1)
	plt.plot(x,y,'o-r')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title('Some Numbers')
	plt.subplot(1,2,2)
	plt.plot(logx,logy,'s-b')
	plt.xlabel('log x')
	plt.ylabel('log y')
	plt.title('Some Numbers (log)')
	plt.tight_layout()
	plt.savefig('numbers.png')
	print('wrote to numbers.png')
	return

def read_edge_file(f,name): #copied from lab2
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
	HistoMaker(nodes, DDict)
	return edges

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
	DnDict = {} #dictionary of neighbor degrees of all the nodes
	for n in nodes: #for every node
		nset = NDict[n]
		dnsum = 0 #sum of the degrees of the neighbors of a single node
		for i in nset: #for every neighbor of the node
			dnsum = dnsum + int(DDict[i]) #adding the degree of each neighbor to the sum
		dnavg = dnsum/len(nset) #the avg neighbor degree of the node n
		DnDict[n] = dnavg
	return DnDict

def HistoMaker(nodes,DDict):
	DSet = []
	for n in nodes:
		value = DDict[n]
		DSet.append(value)
	counted = Counter(DSet)
	listx = []
	listy = []
	for n in counted:
		if counted[n] != 0:
			logval = math.log(counted[n])
			listx.append(n)
			listy.append(logval)
	plt.xlim([0,100])
	plt.plot(listx,listy,'*r')
	plt.savefig('NewHisto.png')
	print('wrote to numbers.png')
	return



main()