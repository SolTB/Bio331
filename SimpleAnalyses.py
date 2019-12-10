'''
Program: SimpleAnalyses.py
Author: Sol Taylor-Brill
Date of Last Edit: 12/09/19
This program just runs some very basic analyses on the candidates from various methods
'''
import math
import statistics

def main():
	k = 10 #NUMBER OF TOP NODES TO OUTPUT
	candidatelist_list = [] #list containing lists of candidates from each file (weighted + unweighted)
	weightedcandidatelist_list = [] #list containing lists of candidates from each file (WEIGHTED ONLY)
	weightslist_list = [] #list containing lists of weights from each file
	weightD_list = [] #list containing dictionaries of gene -> weight from each file

	#Candidate list files
	CN_f = "Common_Neighbor_Candidates.txt" #weighted, 100
	WN_f = "Weighted_Neighbor_Candidates.txt" #weighted, 100
	Co_f = "Clustering_Candidates.txt" #weighted, 100
	Maddy_f = "Maddy_candidates_10.txt" #unweighted, 10

	weighted_f_list = [CN_f, WN_f, Co_f] #list of txt files with ranked nodes and corresponding weights
	unweighted_f_list = [Maddy_f] #list of txt files with ranked nodes and NO weights

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
	unweighted_intersectionD, unweighted_category_D, UW_allcandidates = intersections(candidatelist_list) #dictionary of node-> how many lists it appears in + dict #methods -> list of nodes
	weighted_intersectionD, weighted_category_D, W_allcandidates = intersections(weightedcandidatelist_list) #with weighted lists only for numeric analyses

	Weight_Means(weighted_intersectionD, W_allcandidates, weightD_list, k)
	
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
# and all_candidates_set (a set containing all unique nodes)
def intersections(candidatelist_list):
	all_candidates_set = set()
	intersectionD = {} #dictionary of node->how many methods suggest it
	intersection_nums = [] #list of number of intersections there are (ie categories)
	intersection_list = [] #list of (node, #lists) tuples
	category_D = {} # dictionary of number of methods-> list of nodes
	seennodes = set()

	#Makes a dictionary of nodes > #methods that include them
	for l in candidatelist_list: #go through all the candidate lists
		for node in l:
			if node not in seennodes:
				all_candidates_set.add(node) #adds nodes to a list of all unique nodes
				nodesum = 0 #number of methods that include the node
				for l2 in candidatelist_list: #run through the other lists
					if node in l2 and node not in seennodes:
						nodesum = nodesum + 1 #if it's in a list sum = +1
			
				seennodes.add(node) #node is added  after its checked to seennode set so it won't be run over a million times
				intersectionD[node] = nodesum #node = nodesum
				intersection_list.append((nodesum, node))  # + (#methods, node)

			if nodesum not in intersection_nums: #if a new number of methods
				intersection_nums.append(nodesum)

					

	#sorts the nodes into categories based on how many methods include them
	for x in intersection_nums: #ie [4,3,2,1]
		sublist = [] #list for the category
		for node in all_candidates_set: #every node in the set of all nodes
			if intersectionD[node] == x: #if the node is in the category of #methods
				sublist.append(node)#add the node to the list
		category_D[x] = sublist

	return intersectionD, category_D, all_candidates_set



##NOTE: A WAY TO COMBINE BOTH WEIGHT + INTERSECTION DATA IS TO DIVIDE BUT TOTAL NUMBER OF LISTS
## RATHER THAN JUST THE NUMBER OF LISTS THAT THE NODE IS IN TO NATURALLY DECREASE THE WEIGHT OF
## NODES THAT ARE IN FEWER METHODS. COME BACK TO THIS

#Input: weighted_intersectionD (node->#methods dictionary), W_allcandidates (all unique nodes in the weighted lists)
# weightD_list (list of node->weight dictionaries from each method)
#Outputs = MeanD (dictionary of mean weight of each node), TopMeanList (list of top nodes sorted by meanweight), TopIntMeans (list of top nodes sorted by sumweight/#ALL methods)
def Weight_Means(weighted_intersectionD, W_allcandidates, weightD_list, k):
	MeanD = {} #dictionary of mean weight of each node
	listofmeans = [] #list of (mean,node) tuples
	intersection_listofmeans = []
	
	#goes through every unique node and finds the mean weight across methods
	for node in W_allcandidates: 
		sumweight = 0.0 #combined weight from all methods
		for D in weightD_list: #goes through all the node ->weight dictionaries
			if node in D:
				sumweight = sumweight + float(D[node]) #add weight to sum weight
		meanweight = sumweight/weighted_intersectionD[node] # meanweight = Sum of weights/#methods
		listofmeans.append((meanweight,node)) #add (mean,node) to list
		MeanD[node] = meanweight # add node-> mean weight to dictionary
		intmean = sumweight/len(weightD_list) #sum weight/ divided by # lists
		intersection_listofmeans.append((intmean, node))

	listofmeans.sort(reverse=True)
	TopMeanList = listofmeans[0:k]

	intersection_listofmeans.sort(reverse = True)
	TopIntMeans = intersection_listofmeans[0:k]

	return MeanD, TopMeanList, TopIntMeans










main()