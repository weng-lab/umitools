#!/usr/bin/env python3
import networkx as nx
import editdistance as ed


class UmiGraph:
    def __init__(self, uc, max_ed=1):
        '''Given a dictionary in which UMI is key and UMI count is the 
        value, this function constructs a directed graph. Whenever two
        nodes (a, b) have the edit distance <= max_ed AND a >= 2 * b -1, there
        is a edge from a to b.
'''
        self.UG = nx.DiGraph()
        for u in uc:
            self.UG.add_node(u, count=uc[u])
        for u in uc:
            for u2 in uc:
                if u != u2 and ed.eval(u, u2) <= max_ed:
                    if self.get_count(u) >= 2 * self.get_count(u2) - 1:
                        # The edge direction:
                        # from large node (u) to small node (u2)
                        self.UG.add_edge(u, u2)

    def number_true_umi(self):
        '''This function returns the number of UMIs after collapsing small nodes.
        In other words, it returns the number of weakly connected components
        in a directed graph
'''
        return nx.algorithms.components.number_weakly_connected_components(self.UG)

    def get_umi_clusters(self):
        """Prints the number of clusters and return the main UMI for each cluster
        """
        ret = []
        for comp in nx.algorithms.components.weakly_connected_components(self.UG):
            ret.append([(u, self.get_count(u)) for u in comp])
        return ret

    def get_repr_umi(self):
        """It returns the representative UMIs. In other words, it just report the
        true UMIs and does not report the UMIs caused by PCR duplicates
        """
        ret = []
        for comp in nx.weakly_connected_components(self.UG):
            max_count = 0
            max_u = 0
            for u in comp:
                cur_u = self.UG.nodes[u]
                if cur_u["count"] > max_count:
                    max_count = cur_u["count"]
                    max_u = u
            ret.append(max_u)
        return ret
            
    def print_graph(self):
        for i in self.UG.nodes:
            print("Edge(s) from {}".format(i))
            # for j in self.UG.successors(i):
            for j in self.UG.neighbors(i):            
                print("{}({}) -> {}({})".format(i, self.get_count(i), j, self.get_count(j)))
            print("-" * 80)

    def get_count(self, u):
        return self.UG.nodes[u]["count"]

#     def simplify_graph(self):
#         '''It simplfies the graph by collapsing small nodes into large nodes
#         when they are connected by edges. It finishes when there are no edges
#         left in the graph.
# '''
#         # This can be done with number_connected_components()
#         pass
#         while len(self.UG.edges) > 0:
#             # Strategy: find the smallest terminal node, collapse it into its neighbor
#             pass


def main():
    
    uc = {"AAAA": 1, "AAAC": 2, "AAAG": 1, "AAAT": 10,
          "GGGG": 100, "GGGA": 99, "GGGC": 50, "GGGT": 2,
          "CCCC": 101}
    G = UmiGraph(uc, max_ed=1)
    print("Number of true UMIs for this locus: {}".format(G.number_true_umi()))
    print("Subgraphs of weakly connected components")
    clusters = G.get_umi_clusters()
    for c in clusters:
        print(c)

    print("Representative UMIs for this UMI graph:")
    print(G.get_repr_umi())
    
if __name__ == "__main__":
    main()


    
