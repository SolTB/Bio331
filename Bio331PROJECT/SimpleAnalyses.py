'''
Program: SimpleAnalyses.py
Author: Sol Taylor-Brill
Date of Last Edit: 12/09/19
This program just runs some very basic analyses on the candidates from various methods
'''
import math
import statistics

def main():
	candidatelist_list = [] #list containing lists of candidates from each file (weighted + unweighted)
	weightedcandidatelist_list = [] #list containing lists of candidates from each file (WEIGHTED ONLY)
	weightslist_list = [] #list containing lists of weights from each file
	weightD_list = [] #list containing dictionaries of gene -> weight from each file

	#Candidate list files
	CN_f = "Common_Neighbor_Candidates.txt"
	WN_f = "Weighted_Neighbor_Candidates.txt"
	Co_f = "Clustering_Candidates.txt"
	Maddy_f = "Maddy_candidates_10.txt"

	weighted_f_list = [CN_f, WN_f, Co_f]
	unweighted_f_list = [Maddy_f]

	#Runs through all the weighted files to compile master lists
	for f in weighted_f_list:
		candidates, weights, WeightD = candidate_reader(f)
		weightedcandidatelist_list.append(candidates) #candidates added to weighted only list
		candidatelist_list.append(candidates) #list added to weighted + unweighted list of lists
		weightslist_list.append(weights) #list added to list of lists
		weightD_list.append(WeightD)

	#Runs through all the unweighted files to compile masterlists
	for f in unweighted_f_list:
		candidates = unweighted_candidate_reader(f)
		candidatelist_list.append(candidates) #adds candidate list to masterlist

	#Analyses
	intersectionD, category_D = intersections(candidatelist_list) #dictionary of node-> how many lists it appears in + dict #methods -> list of nodes
	
	print(category_D)
	return

#Input: txt file containing unweighted ranked list of candidates
#Output: ordered list of candidates
def unweighted_candidate_reader(candidate_f):
	candidates = [] #list of candidates from first to last

	f= open(candidate_f,'r') #opens the graph file
	print("Opening " + candidate_f)

	for line in f: #goes through every gene in the file
		l = str(line)
		l_sublist = l.split()
		candidates.append(l_sublist[0]) #Add gene to list of candidates

	return candidates

#Input: candidate file which contains two rows: candidates (from highest ranked down) and weights of the candidates. 
#Output = candidates (list), weights(list), WeightD (candidate -> weight dictionary)
def candidate_reader(candidate_f):
	candidates = [] #list of candidates from first to last
	weights = [] #list of weights from first to last
	WeightD = {} #candidate to weight dictionary

	f= open(candidate_f,'r') #opens the graph file
	print("Opening " + candidate_f)

	for line in f: #goes through every gene in the file
		l = str(line)
		l_sublist = l.split()
		candidates.append(l_sublist[0]) #Add gene to list of candidates
		weights.append(l_sublist[1]) #Add weight to list of weights
		WeightD[l_sublist[0]] = l_sublist[1]

	return candidates, weights, WeightD

#Input: List of lists of candidates
#Output: intersectionD (node->#methods dictionary) and category_D (#methods -> list of corresponding nodes)
def intersections(candidatelist_list):
	intersectionD = {} #dictionary of node->how many methods suggest it
	intersection_nums = [] #list of number of intersections there are (ie categories)
	intersection_list = [] #list of (node, #lists) tuples
	category_D = {} # dictionary of number of methods-> list of nodes

	#Makes a dictionary of nodes > #methods that include them
	for l in candidatelist_list: #go through all the candidate lists
		for node in l:
			nodesum = 0 #number of methods that include the node
			for l2 in candidatelist_list: #check if the node is in any of the other candidate lists
				if node in l2:
					nodesum = nodesum + 1 #if it's in a list sum = +1
					intersectionD[node] = nodesum #node = nodesum
					intersection_list.append((nodesum, node))  # + (#methods, node)
					if nodesum not in intersection_nums: #if a new number of methods
						intersection_nums.append(nodesum)
					l2.remove(node) #after it's been counted its removed from the list so it won't be counted 2x (also faster)
	
	#sorts the nodes into categories based on how many methods include them
	for x in intersection_nums: #ie [4,3,2,1]
		sublist = [] #list for the category
		for node in intersection_list: #every node in the list
			if node[0] == x: #if the node is contained in x number of methods
				sublist.append(node[1])#add the node to the list
		category_D[x] = sublist

	return intersectionD, category_D






main()