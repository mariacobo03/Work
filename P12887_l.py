'''
Minimum spanning trees via hand-coded Kruskal

NOTE: The statement allows for parallel edges
'''

from pytokr import item, items
from heapq import heapify, heappop

# using Union-Find, you can quickly determine whether two vertices belong to the same connected component or not
# Union-Find can help detect cycles efficiently --> you don't want them


class MyUF(dict): # Union Find
	
	def new(self, x): # adds a new element x to the disjoint set data structure
		self[x] = x # sets x as its own parent
	
	def find(self, v):
		if v == self[v]:
			return v # returns the representative (root) of the set that v
		r = self.find(self[v]) # ecursively finds the representative of its parent until it reaches the root
		self[v] = r           # path compression
		return r # returns the representative of the set containing v
		
	def union(self, x, y): # merges the sets containing elements x and y
		# finds the representatives (roots) of the sets containing x and y
		xx = self.find(x)
		yy = self.find(y)
		# sets the representative of the set containing x to be the representative of the set containing y
		self[xx] = yy  # merges the two sets into one

class MyGraph(dict):
	
	def add_edge(self, u, v, w):
		assert u != v # edges should not connect a vertex to itself
		u, v = min(u, v), max(u, v) # store undirected edge in this direction
		if u not in self: # vertex u present in the graph?
			self[u] = dict()
		if v not in self:
			self[v] = dict()
		# checks if an edge between u and v already exists
		if v not in self[u]: # if not 
			self[u][v] = float("inf") # initializes the weight of the edge between u and v to positive infinity 
		self[u][v] = min(self[u][v], w) # in case of parallel edges: keep cheapest
		# only the minimum weight among parallel edges is retained

	def edges(self):
		for u in self:
			for v in self[u]:
				yield u, v, self[u][v]  # edge between vertex u and vertex v, along with the weight of the edge self[u][v]

	def nodes(self):
		for u in self:
			yield u # yields all vertices of the graph

def kruskal(g):
	uf = MyUF()
	# Initializes each vertex in the graph g as a separate set in the Union-Find data structure.
	for u in g.nodes():
		uf.new(u)
	h = list() #  store edges in the form of tuples 
	for u, v, w in g.edges():
		h.append((w, u, v)) # list of tuples containing edges along with their weights from the graph g
	# edge with the minimum weight is always at the top of the heap
	heapify(h) # Converts the list h into a min-heap based on the weights of the edges
	mst = MyGraph() # Initializes an empty graph mst to represent mst
	mst_edges = 0 # keep track of the number of edges added
	cost = 0 # the total cost of the MST
	# loop continues until either all vertices are connected or there are no more edges to process (h)
	while mst_edges < len(g) - 1 and h:
		w, u, v = heappop(h) # Pops the edge with the minimum weight from the min-heap h
		if uf.find(u) != uf.find(v): # check if adding the edge (u,v) would create a cycle, if not:
			mst.add_edge(u, v, w) # Adds the edge (u, v) to the MST with weight w
			mst_edges += 1 # Increments the count of edges
			cost += w # Increments the count of weight
			# to prevent creating a cycle
			uf.union(u, v) # Merges the sets containing vertices u and v, mantain connecivity 
	return mst, cost


for size in items(): # size --> number of items in the collection
	g = MyGraph()  # initialize empty graph
	size = int(size)  # convert to integer
	edges = int(item()) # reading number of edges 
	for _ in range(edges):
		u, v, w = int(item()), int(item()), int(item())
		g.add_edge(u, v, w)
	# ~ for u, v, w in g.edges():
		# ~ print(u, v, g[u][v], w)
	st, cost = kruskal(g) # using Kruskal's algorithm (kruskal(g)) to find the minimum spanning tree (st) of the graph
	print(cost)
	# ~ for u, v, w in st.edges():
		# ~ print(u, v, st[u][v], w)



