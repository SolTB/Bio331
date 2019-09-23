## Lab2 GraphSpace Solution
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
import sys

def main():

	## connect to GraphSpace
	graphspace = GraphSpace('aritz@reed.edu', 'platypus')

	post_test_graph(graphspace)
	post_dolphin_graph(graphspace)

def post(G,gs):
	try:
		graph = gs.update_graph(G)
	except:
		graph = gs.post_graph(G)
	return graph

def post_test_graph(graphspace):
	nodes = ['A','B','C','D','E']
	edges = [['B','E'],['D','E'],['C','A'],['C','E'],['B','C']]

	## create a graph
	G = GSGraph()
	G.set_name('Test Graph')
	G.set_tags(['Lab2'])

	## add nodes:
	for node in nodes:
		G.add_node(node,label=node)
		G.add_node_style(node, shape='ellipse', color='#73C8F3', width=90, height=90)

	## add edges:
	for edge in edges:
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edges[0],edges[1],width=2.0)

	## post graph
	graph = post(G,graphspace)
	print('done posting toy graph')

def post_dolphin_graph(graphspace):
	## read dolphin edges
	edges = []
	nodes = set()
	with open('dolphin_edgelist.txt') as fin:
		for line in fin:
			row = line.strip().split()
			nodes.update(row)
			edges.append(row)
	print(edges)

	## read females
	females = set()
	with open('females.txt') as fin:
		for line in fin:
			females.add(line.strip())

	males = set()
	with open('males.txt') as fin:
		for line in fin:
			males.add(line.strip())

	unknown = set()
	with open('unknown-sex.txt') as fin:
		for line in fin:
			unknown.add(line.strip())

	## read side-floppers
	SF = set()
	with open('side-floppers.txt') as fin:
		for line in fin:
			SF.add(line.strip())

	## read ULTs
	ULT = set()
	with open('upside-down-lobtailers.txt') as fin:
		for line in fin:
			ULT.add(line.strip())


	## create a graph
	G = GSGraph()
	G.set_name('Dolphin Graph')
	G.set_tags(['Lab2'])
	G.set_data(data={'description': '<b>Node Shape</b>: circles: females; rectangles: males; triangles: unknown sex<br><b>Node Color</b> blue: upside-down-lobtailers; green: side-floppers '})


	## add nodes:
	for node in nodes:
		G.add_node(node,label=node)
		if node in unknown:
			node_shape = 'triangle'
		elif node in females:
			node_shape = 'ellipse'
		elif node in males:
			node_shape = 'rectangle'
		else:
			sys.exit('ERROR! '+node+' missing')

		node_color = '#BBBBBB' # gray
		if node in ULT:
			node_color = '#73C8F3'
		if node in SF:
			node_color = '#DAF7A6'
		G.add_node_style(node, shape=node_shape, color=node_color, width=80, height=50)

	## add edges:
	for edge in edges:
		G.add_edge(edge[0],edge[1])
		G.add_edge_style(edges[0],edges[1],width=2.0)

	## post graph
	graph = post(G,graphspace)
	print('done posting dolphin graph')


if __name__ == '__main__':
	main()