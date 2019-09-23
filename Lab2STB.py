from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session

def main():
	#post_test_graph()
	post_dolphin_network()
	return

def post_test_graph(): #Posts a test graph to GraphSpace
	nodes = ['A','B','C','D','E']
	edges = [['A','B'],['A','C'],['A','E'],['B','D'],['B','C'],['D','E'],['E','B']]
	G = GSGraph()
	G.set_name('Sol_Test_Graph')
	G.set_tags(['Lab 2'])
	for n in nodes:
		G.add_node(n, label=n)
		if n=='A' or n == 'E':
			G.add_node_style(n,color='#D2B4DE')
		else:
			G.add_node_style(n,shape='star', color='#F8C407')
	for e in edges:
		G.add_edge(e[0],e[1])
		G.add_edge_style(e[0],e[1], edge_style='dotted')
	graphspace.update_graph(G)
	print("Test Graph Updated")

def post_dolphin_network():
	G = GSGraph()
	G.set_name('Sol_Dolphin_Graph')
	G.set_tags(['Lab 2'])
	G.set_data(data={'description': 'females=pink; males=blue; unknown=purple; ULT=star; SF=rectangle'})
	FemaleList = open("females.txt", 'r')
	females = []
	for n in FemaleList:
		x = n.split()
		females = females + [x[0]]
	MaleList = open("males.txt",'r')
	males = []
	for n in MaleList:
		x = n.split()
		males = males + [x[0]]
	ULT = open("upside-down-lobtailers.txt",'r')
	SF = open("side-floppers.txt")
	ULTlist = []
	SFlist = []
	for n in ULT:
		x = n.split()
		ULTlist = ULTlist + [x[0]]
	for n in SF:
		x = n.split()
		SFlist = SFlist + [x[0]]
	DolphinEdgeList = open("dolphin_edgelist.txt",'r')
	nodes = []
	edges = []
	for l in DolphinEdgeList:
		line = str(l)
		nodesinline = line.split()
		edges = edges + [nodesinline]
		for n in nodesinline:
			if n not in nodes:
				nodes = nodes + [n]
	for dolphin in nodes:
		gender = ""
		behavior = ""
		if dolphin in males:
			gender = "male"
			if dolphin in ULTlist:
				behavior = "upside-down-lobtailer"
			elif dolphin in SFlist:
				behavior = "side-flopper"
			else:
				behavior = "none"
		elif dolphin in females:
			gender = "female"
			if dolphin in ULTlist:
				behavior = "upside-down-lobtailer"
			elif dolphin in SFlist:
				behavior = "side-flopper"
			else:
				behavior = "none"
		else:
			gender = "unknown"
			if dolphin in ULTlist:
				behavior = "upside-down-lobtailer"
			elif dolphin in SFlist:
				behavior = "side-flopper"
			else:
				behavior = "none"
		popupstring = 'name:' + dolphin +' sex:' + gender +' behavior:'+ behavior
		G.add_node(dolphin, label=dolphin, popup=popupstring)
		if dolphin in females:
			if dolphin in ULTlist:
				G.add_node_style(dolphin, color='#F48CDA', shape='star', )
			elif dolphin in SFlist:
				G.add_node_style(dolphin, color='#F48CDA', shape='rectangle')
			else:
				G.add_node_style(dolphin, color='#F48CDA')
				print("Female Added")
		elif dolphin in males:
			if dolphin in ULTlist:
				G.add_node_style(dolphin, color='#77D7F6', shape='star')
			elif dolphin in SFlist:
				G.add_node_style(dolphin, color='#77D7F6', shape='rectangle')
			else:
				G.add_node_style(dolphin, color='#77D7F6')
			print("Male Added")
		else:
			if dolphin in ULTlist:
				G.add_node_style(dolphin, color='#D692EE', shape='star')
			elif dolphin in SFlist:
				G.add_node_style(dolphin, color='#D692EE', shape='rectangle')
			else:
				G.add_node_style(dolphin, color='#D692EE')
			print("Other Added")
	for e in edges:
		G.add_edge(e[0],e[1])
	graphspace.update_graph(G)
	print("Graph Finished Updating")
	

		

main()

