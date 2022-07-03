import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from scipy.linalg import eigvalsh


def plemelij_sokhotski_delta(x):
    return -1/np.pi*1/x.imag


def spectral_density_plemelij_sokhotski(n, x, s):
    somme = 0
    for i in range(n):
        somme += 1/(x-s[i])
    return -1/(n*np.pi)*somme


# Defining the Laplacian matrix of an Erdös-Renyi graph of n nodes with parameter p
def LER(n, p):
    X = np.random.random([n, n])
    X = (X+X.T)/2
    np.fill_diagonal(X, 0)
    A = (X <= p).astype(float)
    L = np.diag(A.sum(axis=0)) - A
    return L


def spectral_density_laplacian_er(n, p, x, ite, eps=1E-1):
    def resolvent_trace(x, lambdai):
        return np.sum([1.0/(xi+1j*eps-x) for xi in lambdai])

    def average_resolvent_trace(x):
        return np.mean([resolvent_trace(x, eigvalsh(LER(n, p))) for r in range(0, ite)])
    return [-1/(np.pi*n)*np.imag(average_resolvent_trace(z)) for z in x]


def rand_prob_node():
    nodes_probs = []
    for node in G.nodes():
        node_degr = G.degree(node)
        # print(node_degr)
        node_proba = node_degr / (2 * len(G.edges()))
        # print("Node proba is: {}".format(node_proba))
        nodes_probs.append(node_proba)
        # print("Nodes probablities: {}".format(nodes_probs))
    random_proba_node = np.random.choice(G.nodes(),p=nodes_probs)
    # print("Randomly selected node is: {}".format(random_proba_node))
    return random_proba_node


def add_edge():
        if len(G.edges()) == 0:
            random_proba_node = 0
        else:
            random_proba_node = rand_prob_node()
        new_edge = (random_proba_node, new_node)
        if new_edge in G.edges():
            print("!ככה לא בונים חומה")
            add_edge()
        else:
            print("!מזל טוב")
            G.add_edge(new_node, random_proba_node)
            print("Edge added: {} {}".format(new_node + 1, random_proba_node))


plt.figure(figsize=(8, 5))
plt.subplot(221)
# Drawing a 125 nodes Barabasi-Albert random graph, with parameter 12 and 85 initial nodes
COLOR = 'blue'
# Get parameters
init_nodes = 85
final_nodes = 125
m_parameter = 12
G = nx.complete_graph(init_nodes)
count = 0
new_node = init_nodes

for f in range(final_nodes - init_nodes):
    print("----------> Step {} <----------".format(count))
    G.add_node(init_nodes + count)
    print("Node added: {}".format(init_nodes + count + 1))
    count += 1
    for e in range(0, m_parameter):
        add_edge()
    new_node += 1
nx.draw(G, alpha=0.3, edge_color=COLOR, node_color=COLOR, node_size=50)


plt.subplot(222)
# Drawing a 125 nodes geometric graph with parameter 0.125
# Use seed when creating the graph for reproducibility
G = nx.random_geometric_graph(125, 0.125, seed=896803)
# position is stored as node attribute data for random_geometric_graph
pos = nx.get_node_attributes(G, "pos")

# find node near center (0.5,0.5)
dmin = 1
ncenter = 0
for n in pos:
    x, y = pos[n]
    d = (x - 0.5) ** 2 + (y - 0.5) ** 2
    if d < dmin:
        ncenter = n
        dmin = d

# color by path length from node near center
p = dict(nx.single_source_shortest_path_length(G, ncenter))
nx.draw_networkx_edges(G, pos, alpha=0.2)
nx.draw_networkx_nodes(
    G,
    pos,
    nodelist=list(p.keys()),
    node_size=80,
    node_color=list(p.values()),
    cmap='viridis',
)

plt.xlim(-0.05, 1.05)
plt.ylim(-0.05, 1.05)
plt.axis("off")

plt.subplot(223)
# Drawing a 125 nodes Watts-Strogatz graph with k=4 and p=0.65
p = 0.65
G = nx.watts_strogatz_graph(125, 4, p, seed=None)
nx.draw(G, alpha=0.3, edge_color=COLOR, node_color=COLOR, node_size=50)

plt.subplot(224)
# Drawing a 4-regular graph (with parameter 0.3)
d = 4
G = nx.random_regular_graph(d, 125, seed=None)
nx.draw(G, alpha=0.3, edge_color=COLOR, node_color=COLOR, node_size=50)

plt.show()
