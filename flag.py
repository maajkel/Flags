from queue import PriorityQueue
import random

from graph import Vertex, Edge, Graph2D
from geometry import Point2D
from kdtree import KdTree

INFINITY_DIST = 1000 * 1000 * 1000

# Single arc-flag - stores a boolean value
class ArcFlag:
    def __init__(self, val=0):
        self.val = val


# Arc-flag vector - stores a vector of arc-flags
# There is exactly one flag for each region
class ArcFlagVector:
    def __init__(self, vec=[]):
        if not vec:
            self.vect = [ 0 for i in range(4) ]
        else:
            self.vect = vec
        
    def get_flag_value(self, index):
        return self.vect[index]


# Represents an edge with additional info (arg-flag vector)
class EdgeWithArcs(Edge):
    def __init__(self, id, vert1, vert2, weight, nr_of_regions=4):
        Edge.__init__(self, id, vert1, vert2, weight)
        self.vec = ArcFlagVector([ 0 for i in range(nr_of_regions) ])
    
    def __str__(self):
        vecStr = ""
        for reg in range(4):
            vecStr += str(self.vec.get_flag_value(reg)) + " "
        
        return "<Edge id:{} vert1:{} vert2:{} weight:{} arc-flag vector=<{}>>".format(self.id, self.vert1, self.vert2, self.weight, vecStr)
    
    def get_arc_flag_value(self, index):
        return self.vec.vect[index]
        
    def update_arc_flag_value(self, index, val):
        self.vec.vect[index] = val


# A graph subclass with the arc-flag vector included
class Graph2DWithArcs(Graph2D):
    def __init__(self):
        Graph2D.__init__(self)
        
    def add_edge(self, id, vert_id1, vert_id2, weight, nr_of_regions):
        Graph2D.add_edge(self, id, vert_id1, vert_id2, weight)
        edge = EdgeWithArcs(id, vert_id1, vert_id2, weight, nr_of_regions)
        self.edges[id] = edge
        self.get_vertex(vert_id1).add_edge(edge)
        self.get_vertex(vert_id2).add_edge(edge)


# Reverses a given graph and returns the reversed graph
# Time complexity: O(n + m)
def reverse_graph(graph):
    result_graph = Graph2DWithArcs()
    
    for vertex in graph.get_vertices():
        result_graph.add_vertex(vertex.id, vertex.location)
       
    for edge in graph.get_edges():
        result_graph.add_edge(edge.id, edge.vert2, edge.vert1, edge.weight, 4)
        
    return result_graph


# Returns a dictionary with (vertex_id, region_nr) pairs
def get_regions(graph, tree=None):
    if tree == None:
        tree = KdTree(graph.get_bounds(), 6)
        for vert in graph.get_vertices():
            tree.add(vert)
        tree.divide()
    
    region_set = dict()
    
    reg_nr = 0
    for region in tree.get_regions():
        print(reg_nr, ": ", end="")
        for point in region.points:
            print(point.id, end = ", ")
            region_set[point.id] = reg_nr
        print()
        reg_nr += 1
        
    return region_set


# Returns the set of boundary nodes for each region in a given graph
def get_boundary_nodes(graph, tree=None):
    region_set = get_regions(graph, tree)
    
    boundary_nodes = dict()
    
    for vertex in graph.get_vertices():
        for edge in vertex.get_edges():
            id_1 , id_2 = edge.vert1 , edge.vert2
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
def dijkstra(graph, s, t=None, use_arcs=False):
    start_region = get_regions(graph)[s.id]
    target_region = 0
    
    if t != None:
        target_region = get_regions(graph)[t.id]
    
    dist = [ INFINITY_DIST for i in range(len(graph.get_vertices())) ]
    dist[s.id] = 0
    
    print("Starting Dijkstra with node " + str(s.id) + " in region " + str(start_region))
    
    pq = PriorityQueue()
    pq.put((dist[s.id] , s.id))
    
    while not pq.empty():
        dist_1 , vert_id =  pq.get()
        
        print(vert_id, "(dist = ", dist_1, ")")
        
        for edge in graph.get_vertex(vert_id).get_edges():
            if use_arcs == True and edge.get_arc_flag_value(target_region - 1) == 0:
                continue
                
            vert_1_id , new_dist = edge.vert1 , dist_1 + edge.weight
            vert_2_id = edge.vert2
            
            if vert_id != vert_1_id:
                vert_2_id = vert_1_id
                
            if new_dist <= dist[vert_2_id]:
                print("Edge " + str(edge) + " is equal/better!")
                edge.update_arc_flag_value(start_region - 1, 1)
                dist[vert_2_id] = new_dist
                pq.put((new_dist , vert_2_id))
    
    #if t == None:
    return dist
    #return dist[t.id]


# Runs a DFS to find one path from s to t. The DFS is run only on edges with t region's arc-flag set to true
def dfs(graph, visited, region_data, s, t, distance):
    print("DFS in", s, "(dist = ", distance, ")")
    
    if s == t:
        return True
    visited.add(s)
    
    good_edges = dict()
    
    counter = 0
    
    for edge in graph.get_vertex(s).get_edges():
        if edge.get_arc_flag_value(region_data[t] - 1) == 1:
            counter += 1
            good_edges[counter] = edge
    
    print("Counter =", counter)
    
    for edge_nr in range(1, counter + 1):
        #chosen_edge_nr = random.randint(1, edge_nr)
        chosen_edge = good_edges[edge_nr]
        
        print("Probuje", edge_nr, ", czyli", chosen_edge.vert1, "<->", chosen_edge.vert2)
    
        new_s = chosen_edge.vert1
   
        if chosen_edge.vert2 != s:
            new_s = chosen_edge.vert2
    
        if new_s not in visited:
            if dfs(graph, visited, region_data, new_s, t, distance + chosen_edge.weight) == True:
                return True
                
                
    return False
    

#### DEBUG ####
    
test_graph_with_arcs = Graph2DWithArcs()

for i in range(13):
    test_graph_with_arcs.add_vertex(i, Point2D(0, i))

test_graph_with_arcs.add_edge(0, 1, 2, 10, 4)
test_graph_with_arcs.add_edge(1, 2, 1, 20, 4)
test_graph_with_arcs.add_edge(2, 5, 7, 30, 4)
test_graph_with_arcs.add_edge(3, 5, 7, 40, 4)
test_graph_with_arcs.add_edge(35, 0, 4, 100, 4)
test_graph_with_arcs.add_edge(4, 0, 1, 50, 4)# 8, 50, 4)
test_graph_with_arcs.add_edge(5, 4, 1, 60, 4)
test_graph_with_arcs.add_edge(6, 7, 8, 70, 4)
test_graph_with_arcs.add_edge(7, 8, 7, 80, 4)
test_graph_with_arcs.add_edge(8, 8, 5, 90, 4)
test_graph_with_arcs.add_edge(9, 7, 2, 100, 4)
test_graph_with_arcs.add_edge(10, 11, 12, 10, 4)
test_graph_with_arcs.add_edge(11, 11, 6, 10, 4)

regions_info = get_regions(test_graph_with_arcs)

print(test_graph_with_arcs)
boundary_nodes_info = get_boundary_nodes(test_graph_with_arcs)
print(boundary_nodes_info)

#for edge in test_graph_with_arcs.get_edges():
#    if dijkstra(test_graph_with_arcs, test_graph_with_arcs.get_vertex(node.id))


### Preprocessing ###

# Setting arc-flag value for edges with equal vertex-regions
for vertex in test_graph_with_arcs.get_vertices():
    for edge in vertex.get_edges():
        v1 , v2 = edge.vert1 , edge.vert2
        
        if regions_info[v1] == regions_info[v2]:
            edge.update_arc_flag_value(regions_info[v1] - 1, 1)

# Setting arc-flag values for tree made by Dijkstra starting from boundary nodes
for region in boundary_nodes_info.keys():
        
    for boundary_node in boundary_nodes_info[region]:
        print("X: ", boundary_node)
        #dijkstra(test_graph_with_arcs, test_graph_with_arcs.get_vertex(boundary_node))
        for edge in test_graph_with_arcs.get_vertex(boundary_node).get_edges():
            v1 , v2 = edge.vert1 , edge.vert2
            print(v1, v2)
            
            if regions_info[v1] != regions_info[v2]:
                dijkstra(test_graph_with_arcs, test_graph_with_arcs.get_vertex(v1))
                dijkstra(test_graph_with_arcs, test_graph_with_arcs.get_vertex(v2))

##### TESTS #####

print(test_graph_with_arcs)
dfs(test_graph_with_arcs, set(), regions_info, 4, 8, 0)

 