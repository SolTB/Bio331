import Lab4_utils
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session

def main():
	node_list, edge_list, adj_list, adj_mat = Lab4_utils.get_graph('lab')
	D = shortest_paths(node_list,edge_list,adj_list, 'A')
	Graph_Maker(node_list, edge_list, D)
	return

def shortest_paths(node_list, edge_list, adj_list, node_s):
	D={}
	for node in node_list:
		D[node] = 1000 #sets the layer to 1000
	D[node_s] = 0 #sets the layer of node_s to 0 since it is the starting node
	Q = [node_s] #the queue starts off with node_s
	while len(Q)>0: #if there are items in Q
		w = Q[0] #the node to be explored is the first item on the list
		Q.remove(w) #first item on the list is removed
		for neighbor in adj_list[w]:
			if D[neighbor] == 1000:
				D[neighbor] = D[w] + 1 #layer set to whatever the previous layer was + 1
				Q.append(neighbor)#adds the neighbor to the end of Q
	print(D)
	return D

def Graph_Maker(node_list, edge_list, D):
	G = GSGraph()
	G.set_name('Sol_BFS_Graph')
	G.set_tags(['Lab 4'])
	Dlist = [] #I'm just making a list of layers to use to calculate node color
	for noden in D: #just making a list of layers (useful for nodecolor)
		if D[noden] not in Dlist:
			Dlist.append(D[noden])
	for n in node_list:
		G.add_node(n, label=n)
		if D[n] == 1000:
			nodecolor = '#FFFFFF' # if node not connected, color = white
		else:
			blue = D[n]/max(Dlist) #first node should result in 0 last layer should result in 1
			red = 1-blue #first node = 1; last layer = 0
			nodecolor = rgb_to_hex(red,0,blue)
		G.add_node_style(n, width = 20, height = 20, color=nodecolor)
	for e in edge_list:
		G.add_edge(e[0],e[1])
	graphspace.update_graph(G)
	print("Graph Updated")

'''
This function just outputs a hex color code from inputs of the amount of red and blue calculated from layer in Graph_Maker()
This function was just copied from the lab handout.
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