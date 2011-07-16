"""
Floyd-Warshall algorithm adapted to finding Longest Paths in Directed
Acyclic Graphs (DAGs). (A Dynamic Programming Approach)
I should investigate using Johnson's Algorithm that solves the problem
for sparse DAGs. 
"""

from util.matrix import print_matrix

def longest_path(graph):
	# Path does not exist for infinite values. 
	inf = float('Infinity')
	DNE = [inf, -inf]

	length = len(graph)

	# Holds the length of the path between i and j as c[i][j]
	# We build this Dynamic Programming matrix as we increment k. 
	cost = [[inf for x in range(length)] for xx in range(length)]

	# Path reconstruction matrix.
	# For the longest path between i and j, path[i][j] is the next
	# vertex in the path. Then path[vert][j], and so forth...
	path = [[inf for x in range(length)] for xx in range(length)]

	# Initialize cost^(0) level.
	for i in range(length):
		for j in range(length):
			cost[i][j] = graph[i][j]

	print "Init\n========="
	print_matrix(cost)
	print "\n"

	# Adapted from Floyd's Algorithm
	for k in range(length):
		for i in range(length):
			for j in range(length):
				#if cost[i][k] in DNE or cost[k][j] in DNE:
				#	continue
				if cost[i][j] > cost[i][k] + cost[k][j]:
				#if cost[i][k] + cost[k][j] < cost[i][j]:
					val = cost[i][k] + cost[k][j]
					x = cost[i][k]
					y = cost[k][j]
					print "%f, %f, %f" % (x, y, val)
					cost[i][j] = cost[i][k] + cost[k][j]
					path[i][j] = k

		print "Level %d\n=========" % k
		print_matrix(cost)
		print "\n"

	return (cost, path)

