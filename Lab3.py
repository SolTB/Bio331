import random
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

def main():
	Output = BA(3,2,50)
	V = Output[0]
	E = Output[1]
	nk = Output[2]
	ek = Output[3]
	GraphMaker(V,E,nk,ek)
	return

def ER(n,m):
	V=[]
	S= set()
	for x in range(n):
		V.append(str(x))
	for i in range(n):
		for j in range(n):
			if i!=j and (V[j],V[i]) not in S:
				S.add((V[i],V[j]))
	S = list(S)
	random.shuffle(S)
	E = S[:m]
	return V,E

def BA(no,mo,t):
	V = []
	for x in range(no):
		V.append(str(x)) #make original node list
	EROutput = ER(no,(no*(no-1))) #making original graph
	E = EROutput[1] #original edges
	d = degree_calculator(V,E) #original degree dictionary
	ekdict = {}
	nkdict = {}
	for n in V:
		nkdict[n] = 0
	for e in E:
		ekdict[e] = 0
	for i in range(t): #for each time step
		n1 = no + i
		u = str(n1)
		pr_list = []#generating a new pr_list each time step
		newnodelist = [] #initiating weighted nodelist each time step
		nkdict[u] = i
		d = degree_calculator(V,E)
		for v in V: #generating a probability list
			numE = len(E)
			dictcall= d[str(v)]
			p=pr(dictcall,numE)
			pr_list = pr_list + [p]
		for x in range(len(pr_list)): #for every value in the probability list
			num = (pr_list[x]/min(pr_list))
			intnum = int(num)
			for j in range(intnum):
				newnodelist.append(V[x]) #adds it proportionally to the probability
		for x in range(mo):
			nodeA = random.choice(newnodelist)
			E.append((u,nodeA))
			ekdict[(u,nodeA)] = i
			newnodelist = remAll(newnodelist,nodeA)
		V = V + [str(u)]
	return V,E, nkdict, ekdict

def remAll(L, item):
    answer = []
    for i in L:
        if i!=item:
            answer.append(i)
    return answer

def degree_calculator(nodes,edges):
	DegreeDict = {}
	for n in nodes: #goes through all the nodes
		degree = 0
		for e in edges: #checks if that node is in an edge by going through each edge
			if n in e:
				degree = degree + 1 #if there is an edge add 1 to the degree
		DegreeDict[n] = degree #creates dictionary of node to degree
	return DegreeDict

def pr(d,E):
	p = d/(2*E)
	return p

def GraphMaker(V,E,nk,ek):
	graphspace = GraphSpace("soltb@reed.edu", "solTB") #Starting GraphSpace session
	G = GSGraph()
	G.set_name('Sol_BA_Graph')
	G.set_tags(['Lab 3'])
	for node in V:
		nkval = nk[node]
		G.add_node(node, label=node, k=nkval)
	for edge in E:
		ekval = ek[edge]
		G.add_edge(edge[0],edge[1], k=ekval)
	graphspace.update_graph(G)
	print("Graph Posted")
main()