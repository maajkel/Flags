from queue import PriorityQueue
import random
import sys
from datetime import datetime
import gc

from graph import Vertex, Edge, Graph2D
from geometry import Point2D
from kdtree import KdTree
from gparser import GraphParser

from heapdict import heapdict

NR_OF_REGIONS = 128

# Arc-flag vector - stores a vector of arc-flags
# There is exactly one bit for each region
class ArcFlagVector:
    def __init__(self, vec=[]):
        if not vec:
            self.vect = 1 << (NR_OF_REGIONS + 4)
        else:
            self.vect = vec
        
    def get_flag_value(self, index):
        return (self.vect & (1 << index)) >> index


# Represents an edge with additional info (arg-flag vector)
class EdgeWithArcs(Edge):
    def __init__(self, id, vert1, vert2, weight, nr_of_regions=4):
        Edge.__init__(self, id, vert1, vert2, weight)
        self.vec = ArcFlagVector()
    
    def __str__(self):
        vecStr = ""
        for reg in range(NR_OF_REGIONS):
            if self.vec.get_flag_value(reg) == 1:
                vecStr += '1'
            else:
                vecStr += '0'
                
        return str(self.id) + " " + str(self.vert1.id) + " " + str(self.vert2.id) + " " + str(self.weight) + " " + vecStr #"<Edge id:{} vert1:{} vert2:{} weight:{} arc-flag vector=<{}>>".format(self.id, self.vert1, self.vert2, self.weight, vecStr)
    
    def get_arc_flag_value(self, index):
        return (self.vec.vect & (1 << index)) >> index
        
    def update_arc_flag_value(self, index, val):     
        self.vec.vect = self.vec.vect | (val << index)


# A graph subclass with the arc-flag vector included
class Graph2DWithArcs(Graph2D):
    def __init__(self):
        Graph2D.__init__(self)
        
    def __init__(self, graph : Graph2D):
        Graph2D.__init__(self)
        
        for vertex in graph.get_vertices():
            self.add_vertex(vertex.id, vertex.location)
        
        for edge in graph.get_edges():
            self.add_edge(edge.id, edge.vert1.id, edge.vert2.id, edge.weight)
    
    def add_edge(self, id, vert_id1, vert_id2, weight, nr_of_regions=4):
        vert1 = self.get_vertex(vert_id1)
        vert2 = self.get_vertex(vert_id2)
        
        Graph2D.add_edge(self, id, vert_id1, vert_id2, weight)
        edge = EdgeWithArcs(id, vert1, vert2, weight, nr_of_regions)
        self.edges[id] = edge
        self.get_vertex(vert_id1).add_edge(edge)
        self.get_vertex(vert_id2).add_edge(edge)

test_graph_with_arcs = Graph2DWithArcs(Graph2D())

# Reverses a given graph and returns the reversed graph
# Time complexity: O(n + m)
def reverse_graph(graph):
    result_graph = Graph2DWithArcs()
    
    for vertex in graph.get_vertices():
        result_graph.add_vertex(vertex.id, vertex.location)
       
    for edge in graph.get_edges():
        result_graph.add_edge(edge.id, edge.vert2, edge.vert1, edge.weight, 4)
        
    return result_graph

regions_info = dict()

# Returns a dictionary with (vertex_id, region_nr) pairs
def get_regions(graph, tree=None, max_accumulation=100000):
    if tree == None:
        tree = KdTree(graph.get_bounds(), max_accumulation)
        for vert in graph.get_vertices():
            tree.add(vert)
        tree.divide()
    
    region_set = dict()
    
    reg_nr = 0
    for region in tree.get_regions():
        for point in region.items:
            region_set[point.id] = reg_nr
        
        reg_nr += 1
        
    return region_set


# Returns the set of boundary nodes for each region in a given graph
def get_boundary_nodes(graph, tree=None):
    region_set = regions_info
    
    boundary_nodes = dict()
    
    for vertex in graph.get_vertices():
        for edge in vertex.get_edges():
            
            id_1 , id_2 = edge.vert1.id , edge.vert2.id
            reg_1 , reg_2 = region_set[id_1] , region_set[id_2]
            
            if reg_1 != reg_2:
                if reg_2 not in boundary_nodes:
                    boundary_nodes[reg_2] = set()
                if reg_1 not in boundary_nodes:
                    boundary_nodes[reg_1] = set()
                boundary_nodes[reg_2].add(id_2)
                boundary_nodes[reg_1].add(id_1)
    
    return boundary_nodes


# Runs Dijkstra's algorithm and returns the distance array or the distance to vertex t, if specified
# Sets arg-flag values for the start vertex region and if t is specified, searches shortest paths using
# precomputed arg-flag values 
def dijkstra(s, t=None, use_arcs=False):
    start_region = regions_info[s.id]
    target_region = 0
    
    vertices_size = len(test_graph_with_arcs.get_vertices())
    
    if t != None:
        target_region = regions_info[t.id]
    
    dist = [ -1 for i in range(vertices_size + 5) ]

    
    visited = [0 for i in range(vertices_size + 5) ]
    dist[s.id] = 0
    
    pq = heapdict()
    pq[s] = dist[s.id]
    
    if use_arcs == False:
        for vertex in test_graph_with_arcs.get_vertices():
            if regions_info[vertex.id] == start_region:
                dist[vertex.id] = 0
                pq[vertex] = dist[vertex.id]
    
    while len(pq) > 0:
        (vert, dist_1) = pq.popitem()
        
        if visited[vert.id] == 1:
            continue
        visited[vert.id] = 1
        
        for edge in vert.get_edges():
            vert_2 = edge.get_other_vert(vert)
                
            if (use_arcs == True and edge.get_arc_flag_value(target_region) == 0):
                continue
            
            new_dist = dist_1 + edge.weight
                
            if dist[vert_2.id] == -1 or new_dist < dist[vert_2.id]:
                if use_arcs == False:
                    edge.update_arc_flag_value(start_region, 1)
                dist[vert_2.id] = new_dist
                pq[vert_2] = new_dist


# Runs a DFS to find one path from s to t. The DFS is run only on edges with t region's arc-flag set to true
def dfs(graph, visited, region_data, s, t, distance):
    print("DFS in", s, "(dist = ", distance, ")")
    
    if s == t:
        return True
    visited.add(s)
    
    good_edges = list()
    
    for edge in graph.get_vertex(s.id).get_edges():
        if edge.get_arc_flag_value(region_data[t.id]) == 1:
            good_edges.append(edge)
    
    random.shuffle(good_edges)
    
    for i in range(len(good_edges)):
        chosen_edge = good_edges[i]
    
        new_s = chosen_edge.get_other_vert(s)
        
        print(i, "Probuje", chosen_edge.id, ", czyli", s.id, "<->", new_s.id)
    
        if new_s not in visited:
            if dfs(graph, visited, region_data, new_s, t, distance + chosen_edge.weight) == True:
                return True
    
    # Every possible edge with region's t bit set to true has been tested
    return False


#### DEBUG ####
    
vert_filename = sys.argv[1]
edge_filename = sys.argv[2]
max_accumulation = int(sys.argv[3])
save_filename = sys.argv[4]

def print_time(message):
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print(message, current_time, flush=True)

graph = Graph2D()

def simple_test():  
    for i in range(13):
        graph.add_vertex(i, Point2D(0, i))

    graph.add_edge(0, 1, 2, 10)
    graph.add_edge(1, 2, 1, 20)
    graph.add_edge(2, 5, 7, 30)
    graph.add_edge(3, 5, 7, 40)
    graph.add_edge(35, 0, 4, 100)
    graph.add_edge(4, 0, 1, 50)# 8, 50, 4)
    graph.add_edge(5, 4, 1, 60)
    graph.add_edge(6, 7, 8, 70)
    graph.add_edge(7, 8, 7, 80)
    graph.add_edge(8, 8, 5, 90)
    graph.add_edge(9, 7, 2, 100)
    graph.add_edge(10, 11, 12, 10)
    graph.add_edge(11, 11, 6, 10)
    
    ea = EdgeWithArcs(17, 18, 19, 20, 4)
    ea.update_arc_flag_value(3, 1)
    for i in range(10):
        print(ea.get_arc_flag_value(i))



#with open(vert_filename) as vert_filestream, open(edge_filename) as edge_filestream:
#    parser = GraphParser()
#    parser.parse_csv(vert_filestream, edge_filestream)

#    graph = parser.get_graph()
 
simple_test()
 
test_graph_with_arcs = Graph2DWithArcs(graph)

del graph
gc.collect()
print("Graph imported", flush=True)

print_time("Current Time 1")

regions_info = get_regions(test_graph_with_arcs, max_accumulation=6)

print(regions_info, flush=True)
print_time("Current Time 2 =")

boundary_nodes_info = get_boundary_nodes(test_graph_with_arcs)

#### PREPROCESSING ####


print_time("Current Time 3 =")
print("Starting preprocessing...", flush=True)

# Setting arc-flag value for edges with equal vertex-regions
for vertex in test_graph_with_arcs.get_vertices():
    for edge in vertex.get_edges():
        v1_id , v2_id = edge.vert1.id , edge.vert2.id
        
        if regions_info[v1_id] == regions_info[v2_id]:
            edge.update_arc_flag_value(regions_info[v1_id], 1)

print_time("Current Time 4 =")

def dijktras_cnt():
    cnt = 0
    
    visited = [False for i in range(len(test_graph_with_arcs.get_vertices()))]

    for region in boundary_nodes_info.keys():
        for boundary_node in boundary_nodes_info[region]:
            for edge in test_graph_with_arcs.get_vertex(boundary_node).get_edges():
                v1 , v2 = edge.vert1 , edge.vert2
                
                if regions_info[v1.id] != regions_info[v2.id]:
                    if not visited[v1.id]:
                        cnt += 1
                        visited[v1.id] = True
                    if not visited[v2.id]:
                        cnt += 1
                        visited[v2.id] = True
    return cnt

sum_cnt = dijktras_cnt()
global_cnt = 0
visited = [False for i in range(len(test_graph_with_arcs.get_vertices()))]

# Run one Dijkstra per region
def dijkstra_with_supervertices():
    for region in boundary_nodes_info.keys():
        for boundary_node in boundary_nodes_info[region]:
            print("Dijkstra from region = " + str(region) + ", node = " + str(boundary_node), flush=True)
            print_time("Current Time =")
            
            dijkstra(test_graph_with_arcs.get_vertex(boundary_node))
            break

# Setting arc-flag values for tree made by Dijkstra starting from boundary nodes
def old_dijkstra():
    for region in boundary_nodes_info.keys():
        for boundary_node in boundary_nodes_info[region]:
            for edge in test_graph_with_arcs.get_vertex(boundary_node).get_edges():
                v1 , v2 = edge.vert1 , edge.vert2
                
                if regions_info[v1.id] != regions_info[v2.id]:
                    if not visited[v1.id]:
                        dijkstra(test_graph_with_arcs.get_vertex(v1.id))
                        visited[v1.id] = True
                        global_cnt += 1
                    if not visited[v2.id]:
                        dijkstra(test_graph_with_arcs.get_vertex(v2.id))
                        visited[v2.id] = True
                        global_cnt += 1

dijkstra_with_supervertices()

##### TESTS #####

print_time("Current Time 5 =")

for edge in test_graph_with_arcs.get_edges():
    print(edge)

dijkstra(test_graph_with_arcs.get_vertex(8), test_graph_with_arcs.get_vertex(4), use_arcs=True)

#print(test_graph_with_arcs, flush=True)
dfs(test_graph_with_arcs, set(), regions_info, test_graph_with_arcs.get_vertex(1), test_graph_with_arcs.get_vertex(5), 0)

 