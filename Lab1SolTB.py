## Lab1 -- Warm Up
## September 4 2019
## Requirements: Python3

def main():
	"""
	Main function. There are no inputs 
	(nodes and matrix are specified within this function)
	"""

	## List of node names
	nodes = ['A','B','C','D','E','F']
	## Adjacency matrix (assume this represents an undirected graph & has same order as nodes)
	adj_mat = [[0,1,0,1,0,1],[1,0,0,0,1,1],[0,0,0,1,0,0],[1,0,1,0,0,0],[0,1,0,0,1,1],[1,1,0,0,1,0]] 
	#adj_mat_challenge = [[0,1,0,1,0,1],[1,0,0,0,1,1],[0,0,0,1,0,0],[1,0,1,0,0,0],[0,1,0,0,1,1],[1,1,0,0,1,0]]

	## function calls
	print('\nAdjacency Matrix:')
	print_mat(nodes,adj_mat)
	print("Number of Edges from Adj. Matrix:",num_edges_from_mat(adj_mat))

	print('\nMaking Adjacency List...')
	adj_list = mat_to_list(nodes,adj_mat)

	print('\nAdjacency List:')
	print_list(adj_list)
	print("Number of Edges from Adj. List:",num_edges_from_list(adj_list))

	output = list_to_mat(adj_list) #Optional Challenge Function
	print("\nNode List from Adj List:", output[0])
	print("\nAdj Mat from Adj List:", output[1])
	print('\nDone!')
	return # done with main function

def print_mat(nodes,adj_mat):
	"""
	Prints the adjacency matrix
	Inputs: nodes (list of strings) and adj_mat (list of lists)
	Returns: nothing
	"""
	LineA = " "
	for n in nodes:
		LineA = LineA + "  " + n
	print(LineA)
	for i in range(len(nodes)):
		print(nodes[i],adj_mat[i])	
	return

def num_edges_from_mat(adj_mat):
	"""
	Counts the number of edges from the adjacency matrix
	Inputs: adj_mat (list of lists)
	Returns: the number of edges (int)
	"""
	symedge = 0
	selfedge = 0
	for x in range(len(adj_mat)):
		for y in range(len(adj_mat[x])):
			if adj_mat[x][y] == 1:
				if x==y: #if it's on the diagonal (a selfedge)
					selfedge = selfedge + 1
				else: #if it's in the main part of the matrix - symmetrical edge bewtween two nodes
					symedge = symedge + 1
	realedge = selfedge + (symedge/2)

	return realedge

def mat_to_list(nodes,adj_mat):
	"""
	Converts the adjacency matrix to an adjacency list
	Inputs: nodes (list of strings) and adj_mat (list of lists)
	Returns: adjacency list (dictionary of (node,neighbor list) pairs).
	"""
	adj_dict = {}
	for x in range(len(nodes)):
		entry = []
		for y in range(len(adj_mat[x])):
			if adj_mat[x][y] == 1:
				entry = entry + [nodes[y]]
		adj_dict[nodes[x]] = entry
	return adj_dict

def print_list(adj_list):
	"""
	Prints the adjacency list
	Inputs: adj_list (dictionary of (node,list) pairs)
	Returns: nothing
	"""
	for x in adj_list:
		stringboy = "" #initializing a string to make the adj list look pretty
		for n in adj_list[x]:
			stringboy = stringboy + n + " " 
		call = x + ":"
		print(call, stringboy)

	return 

def num_edges_from_list(adj_list):
	"""
	Counts the number of edges from the adjacency list
	Inputs: adj_list (dictionary of (node,list) pairs)
	Returns: the number of edges (int)
	"""
	normedgedouble= 0
	selfedge = 0
	for x in adj_list:
		degree = 0
		for y in adj_list[x]:
			if x == y:
				selfedge = selfedge + 1
			else:
				normedgedouble = normedgedouble + 1
	normedge = normedgedouble/2
	numedges = normedge + selfedge
	return numedges

def list_to_mat(adj_list): #Optional Challenge Function
	nodes = []
	new_adj_mat = []
	for x in adj_list:
		nodes = nodes + [x] #I start by making an ordered list of nodes
	for x in adj_list:
		xlist = []
		for nodenum in range(len(nodes)):
			if nodes[nodenum] in adj_list[x]: #Then I go through the ordered list of nodes and for every node I ask if there is a corresponding value in the dictionary output
				xlist = xlist + [1] #if that node is there add a 1
			else:
				xlist = xlist + [0] # if not add a 0
		new_adj_mat = new_adj_mat + [xlist] #making a list of lists
	return nodes, new_adj_mat

"""
 Leave this is at the bottom of the file. Once all functions are loaded, then 
 main() is called UNLESS you are importing this file into another script. 
 See https://docs.python.org/3/library/__main__.html
"""
if __name__ == '__main__':
	main()
