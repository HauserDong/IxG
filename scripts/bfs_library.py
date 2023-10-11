from collections import deque
import util


def bfs_with_paths(graph, start):
    visited = set()
    queue = deque([(start, [start])])
    paths = {vertex: [] for vertex in graph}

    while queue:
        vertex, path = queue.popleft()
        visited.add(vertex)
        paths[vertex] = path

        for neighbor in graph[vertex]:
            if neighbor not in visited:
                new_path = path + [neighbor]
                queue.append((neighbor, new_path))

    return paths


# Function to convert a list of edges to an adjacency dictionary
def edges_to_adjacency_dict(edges):
    adjacency_dict = {}
    for u, v in edges:
        adjacency_dict.setdefault(u, []).append(v)
        adjacency_dict.setdefault(v, []).append(u)  # If the graph is undirected

    return adjacency_dict


# Load edges
edges = util.DeserializeEdges('../examples/insatxgcs/resources/maze2d/maze_edges.csv')

# Convert edges to an adjacency dictionary
graph = edges_to_adjacency_dict(edges)

path_dict = {}
for st_node in graph:
    print(f"BFS starting from vertex '{st_node}':")
    node_path_dict = bfs_with_paths(graph, st_node)
    for go_node in node_path_dict:
        path_dict[[st_node, go_node]] = node_path_dict[go_node]

dump_list = []
for terms in path_dict:
    l = terms
    l.append(path_dict[terms])
    dump_list.append(l)

# Specify the file path where you want to write the data
file_path = '../examples/insatxgcs/resources/maze2d/paths.txt'

# Open the file in write mode
with open(file_path, 'w') as file:
    for inner_list in dump_list:
        line = ' '.join(map(str, inner_list))  # Convert inner list elements to strings and join with spaces
        file.write(line + '\n')  # Write the line to the file and add a newline character



