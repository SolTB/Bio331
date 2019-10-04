'''
Lab 5: Random Walks
Author: Sol Taylor-Brill
10/02/19
'''
import matplotlib.pyplot as plt

def main():
	foutput = FileReader("test-edges.txt")
	nodes = foutput[0]
	edges = foutput[1]
	Gin, Gout = Neighbor_Dicts(nodes,edges)
	X = rw_probs(nodes,Gin,Gout,"A",25,1)
	Y = transpose(X)
	make_heatmap(Y,"HeatMap10.png", nodes)
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

# This function intakes a list of nodes, the Gin/Gout dictionaries, a starting node, s, timesteps, T, and probability value, q.
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
	print(table_list)
	return table_list

#X = Y from transpose(X)
#copied from Lab 5 instructions
def make_heatmap(X,figname,ordered_nodes):
	plt.figure(figsize=(5,3)) # make a 6.5" wide by 4" tall figure. 
	plt.imshow(X,interpolation='nearest') # plots the table 
	plt.colorbar(shrink=.5) ## Add a colorbar

	## Add x label and x ticks for every column.
	plt.xlabel('Time Step')
	plt.xticks(range(len(X[0])))
        
    ## Add y label and row names.
	plt.ylabel('Gene')
	plt.yticks(range(len(ordered_nodes)),ordered_nodes)
        
    ## Tighten the layout and save the figure.
	plt.tight_layout()
	plt.savefig(figname)
	print('Wrote to ',figname)
	return

#X = table_list
#copied from Lab 5 instructions
def transpose(X):
	Y = []
	for j in range(len(X[0])):  # loop through cols first
		Y += [[0]*len(X)]   # make row of 0's the right size
		for i in range(len(X)):  # loop through rows
			Y[j][i] = X[i][j]  # update entries
	return Y

main()