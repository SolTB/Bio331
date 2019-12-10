'''
Part C: Finding nodes with similar degrees to known positives
Author: Sol Taylor-Brill
Starting Date: 11/27/19
Last Edit Date: 11/27/19
'''
import statistics 

def main():
	positive_f = 'labeled_nodes.txt'
	All_f = 'interactome-flybase-collapsed-weighted.txt' 
	Outputf = "PartC_Candidates.txt"

	Pos_Name_set = PositiveReader(positive_f) #makes set of all known positives
	nodeset, edgeset = InteractomeReader(All_f) #makes edgeset and nodeset from interactome file

	DDict, NDict = DegreeCalc(nodeset,edgeset) #makes dictionaries of neighbors and degree for each node
	rangemin, rangemax = degree_of_positives(Pos_Name_set,DDict) #calculates a certain minimum and maximum degree value from degrees of known positives
	candidateset = node_picker(nodeset, Pos_Name_set, rangemin, rangemax, DDict) #finds a set of candidate genes whos degree is within the set values

	output_maker(candidateset, Outputf) #prints candidates to output file
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
	f= open(All_f,'r') #opens the graph file
	print("Opening " + All_f)

	for line in f:
		l = str(line)
		l_sublist = l.split()
		node1 = l_sublist[0]
		node2 = l_sublist[1]
		weight = l_sublist[2]
		edge = (node1, node2, weight) #a tuple

		if node1 not in nodeset:
			nodeset.add(node1)
		if node2 not in nodeset:
			nodeset.add(node2)

		if edge not in edgeset and node1 != '#symbol1':
			edgeset.add(edge)
	
	nodeset.discard('#symbol1')
	nodeset.discard('symbol2')

	return nodeset, edgeset

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
	return DDict, NDict,

#This function takes the list of positives and the dictionary of node-> degree as inputs
#Outputs the mean, min and max of node degrees of confirmed positive
def degree_of_positives(Pos_Name_set,DDict):
	degreelist = []
	for name in Pos_Name_set:
		degree = DDict[name]
		degreelist.append(degree)
		if degree == 0:
			print(name)
	
	degreelist.sort()
	print("degree list ", degreelist)

	minimum = min(degreelist)
	maximum = max(degreelist)
	print("minimum degree: ", minimum, "maximum degree: ", maximum)
	std = statistics.stdev(degreelist)
	med = statistics.mean(degreelist)
	print("standard deviation: ", std, "mean: ", med)
	rangemax = float(med) + (float(std)/10)
	rangemin = float(med) - (float(std)/10)
	return rangemin, rangemax

#goes through all the nodes and returns a list of nodes that have a degree within the range calculated in degree_of_positives
def node_picker(nodeset, Pos_Name_set, rangemin, rangemax, DDict):
	candidateset = set()

	for node in nodeset:
		if node not in Pos_Name_set:
			degree = DDict[node]
			if degree <= rangemax and degree >= rangemin:
				candidateset.add(node)
	print("Number of candidates: ", len(candidateset))
	return candidateset

def output_maker(candidateset, outputf):
	f= open(outputf,'w')
	print("Opening " + outputf)
	candidatelist = list(candidateset)
	cutcandidates = candidatelist[0:100]

	for candidate in cutcandidates:
		f.write(candidate + '\n')

	print("printed to Output file. Done")


main()
