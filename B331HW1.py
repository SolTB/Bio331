'''
Biology 331: HW 1: Badger Networks
Sol Taylor-Brill

This code constructs a graph of badger interactions on GraphSpace. The graph shows badgers grouped by social group (shown as color)
Females are shown as circles, males are rectangles. Size of node is equal to its degree. Edge thickness is proportional to contact time
between the two individuals
'''
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

import math

def main():
	matfile = "BadgerMatrix.txt" #Adjacency matrix text file
	matrixoutput = read_matrix(matfile) #txt file -> nodelist and adjacency matrix
	nodelist = matrixoutput[0]
	adjacency_matrix = matrixoutput[1]
	check_symmetric(adjacency_matrix) #checks if matrix is symmetrical
	edgeoutput = mat_to_edgelist(adjacency_matrix, nodelist) #adj_mat + nodelist -> edgelist + edgeDict
	edgelist = edgeoutput[0]
	edgeDict = edgeoutput[1]
	BInfo = parse_info("BadgerInfo.txt") #returns 3 dictionaries of badger sex, TB status and social group
	SexDict = BInfo[0]
	TBDict = BInfo[1]
	GroupDict = BInfo[2] 
	DegreeDict = degree_calculator(nodelist,edgelist) #edges + nodes -> Dictionary[node]-> degree
	post_graph(nodelist,edgelist, SexDict, TBDict, GroupDict, DegreeDict, edgeDict) #posts graph to GraphSpace (or updates current graph)
	return

'''
This function reads the matrix txt file
and returns a list of nodes (list of strings) & adjacency matrix (list of lists)
'''
def read_matrix(matfile): 
	matrix = open(matfile, 'r')
	lines = matrix.readlines() #makes a list of lists containing every line in the file
	listoflines= [] #initiates a list that will be a list of lists with everyline in the badger file.
	nodelist = [] #initiates node list
	adj_mat = [] #initiates adjacency matrix (list of lists)
	for line in lines: 
		fileline = line.split() #cleans up lines to remove /t and make every object an individual string in a list
		listoflines = listoflines + [fileline] #adds each line to the list of lists as a new list
	for n in range(len(listoflines[0])):
		if n > 0: #so that "badger" will not be in the nodelist
			nodelist = nodelist + [listoflines[0][n]] #makes node list from the header
	for n in range(len(listoflines)):
		if n >0: #so that the header isn't included
			adj_mat = adj_mat + [listoflines[n][1:]]#badger name which is the first item in the list (index 0) not included
	return nodelist, adj_mat

'''
This function checks if the adjacency matrix is symmetric 
(aka A[i][j] = A[j][i] for all pairs of indices i and j)
returns True if symmetric, False if not
'''
def check_symmetric(A):
	for i in range(len(A)): #for every line (a list) in the adjacency matrix
		for j in range(len(A[i])): #for every character in the line
			if A[i][j] != A[j][i]:
				print("No")
				return False #If there is an asymmetrical value return False
	print("Yes")
	return True #Otherwise, return True

'''
Takes adjacency matrix, A, and nodelist, L, and returns a list of edges (list of 2 element lists) and a dictionary of edge weights (contact time)
'''
def mat_to_edgelist(A,L):
	edgelist = [] #initiates list of lists
	EdgeDict = {}
	for i in range(len(A)):
		for j in range(len(A[i])):
			if int(A[i][j])>0:
				if [L[j],L[i]] not in edgelist:
					edgelist = edgelist + [[L[i],L[j]]]
					dictcall = L[i] + L[j] #concatenates 2 badger names into one string
					EdgeDict[dictcall] = int(A[i][j])
	return edgelist, EdgeDict

'''
This function takes "BadgerInfo.txt" as input and outputs three dictionaries: Name-> sex, TB status and social group
'''
def parse_info(InfoF):
	BadgerInfo = open(InfoF, 'r')
	listoflines = []
	lines = BadgerInfo.readlines() #list of lines in Badger info file
	SexDict = {} #Name -> Sex dictionary
	TBDict = {} #Name -> TB status dictionary
	GroupDict = {} #Name -> Social Group dictionary
	for line in lines: #This part is copied from read_matrix() above
		fileline = line.split() #cleans up lines to remove /t and make every object an individual string in a list
		listoflines = listoflines + [fileline] #adds each line to the list of lists as a new list
	for l in range(len(listoflines)):
		if l>0: #ignores header
			name = listoflines[l][0] #local variable
			SexDict[name] = listoflines[l][1]
			TBDict[name] = listoflines[l][2]
			GroupDict[name] = listoflines[l][3]
	return SexDict, TBDict, GroupDict

'''
This function takes nodelist and edgelist as inputs and outputs a dictionary of the degree of each node
I just didn't want this making post_graph() difficult to follow so I decided to make it its own function.
'''
def degree_calculator(nodes,edges):
	DegreeDict = {}
	for n in nodes: #goes through all the nodes
		degree = 0
		for e in edges: #checks if that node is in an edge by going through each edge
			if n in e:
				degree = degree + 1 #if there is an edge add 1 to the degree
		DegreeDict[n] = degree #creates dictionary of node to degree
	return DegreeDict
'''
Posts graph to GraphSpace
Takes nodelist(from read_matrix()), edgelist and edgeDict(from mat_to_edgelist()), SexDict, TBDict, and Group Dict (from parseinfo()) 
and DegreeDict (from degree_calculator()) as inputs.
'''
def post_graph(nodes,edges,SexDict,TBDict,GroupDict,DegreeDict,edgeDict):
	graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session
	G = GSGraph()
	G.set_name('Sol_Badger_Graph')
	G.set_tags(['HW 1'])
	G.set_data(data={'description': 'Male: Rectangle; Female: Ellipse<br> TB+: Dotted Border; TB-: Solid Border<br> Group1: pink; Group2: red; Group3: orange; Group4: yellow; Group5: green;<br> Group6: light blue; Group7: indigo; Group8: purple <br> View in "Social_Groups" Layout'})
	contactlist = [] #I'm going to use this list to come up with k values for the edges
	for item in edgeDict:
		if edgeDict[item] not in contactlist: #each contact time value only added once
			contactlist = contactlist + [edgeDict[item]] #making a list with all the contact times
	contactlist.sort(reverse=True) #sorting the list of contact times from largest to samllest
	for node in nodes: #adds all of the nodes
		sex = SexDict[node]
		TB = TBDict[node]
		Group = GroupDict[node]
		size = DegreeDict[node]*3
		popupstring = "<i>Sex</i>: " + sex + "<br> <i>TB Status</i>: " + TB + "<br> <i>Group</i>: " + Group
		G.add_node(node, label=node, popup=popupstring, k=0)
		if sex == "Male":
			nshape = 'rectangle'
		elif sex == "Female":
			nshape = 'ellipse'
		if TB == 'P':
			nborder = 'dotted'
		elif TB == 'N':
			nborder = 'solid'
		if Group == '1':
			ncolor = '#EF0FA6' #pink
		elif Group == '2':
			ncolor = '#EF0F16' #red
		elif Group == '3':
			ncolor = '#EF7A0F' #orange
		elif Group == '4':
			ncolor = '#EFDA0F' #yellow
		elif Group == '5':
			ncolor = '#0FEF38' #green
		elif Group == '6':
			ncolor = '#0FEFEC' #light blue
		elif Group == '7':
			ncolor = '#000ACB' #indigo
		elif Group == '8':
			ncolor = '#A614E2' #purple
		G.add_node_style(node, color= ncolor, style=nborder, border_width = 5, shape=nshape, height=size, width=size)
	for edge in edges: #adds all of the edges
		callkey = edge[0] + edge[1] #constructs the keystring for the dictionary
		thickness = edgeDict[callkey] #badger contact time
		kval = contactlist.index(thickness) #k value of edge = index of thickness in the list of contact time (sorted from biggest to smallest)
		logthick = math.log(thickness) #transforms it to log scale
		intlog = int(logthick) #easier to work with than a float
		popstring = "<i>Contact Time</i>: " + str(thickness) + "<br> <i>log(contact time)</i>: " + str(intlog)
		G.add_edge(edge[0],edge[1], popup=popstring, k=kval) #shows contact time and log(contact time) on edge popup
		G.add_edge_style(edge[0],edge[1], width=(intlog/2)) #makes edge thickness = 1/2 of log(contact time)
	graphspace.update_graph(G) #makes updates to the existing graph
	print("Graph Updated")

main()