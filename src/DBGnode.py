class DBGnode:

    def __init__(self, k_mere_seq):
        '''Initialize a de Bruijn graph node with its k-mer sequence'''
        self.k_mere_seq = k_mere_seq
        self.edges_to = {}
        self.edges_from = {}
        
    def add_edge_to(self, eto):
        """
        Add an incoming edge from the DBGnode eto to this node.
        If the edge already exists, its weight is increased by 1.
        eto must be an instance of DBGnode (not just a sequence string).
        """
        self.edges_to[eto] = self.edges_to.get(eto, 0) + 1

    def add_edge_from(self, efrom):
        """
        Add an outgoing edge from this node to the DBGnode efrom.
        If the edge already exists, its weight is increased by 1.
        efrom must be an instance of DBGnode (not just a sequence string).
        """
        self.edges_from[efrom] = self.edges_from.get(efrom, 0) + 1

    def get_potential_from(self):
        """
        Return a list of all k-mer sequences that could theoretically
        have an edge leading to this node.
        The list contains 4 strings, created by prepending A, G, T, or C
        and removing the last base.
        """
        potential_from = []
        for base in ["A", "G", "T", "C"]:
            potential_from.append(base + self.k_mere_seq[:-1])
        return potential_from

    def get_potential_to(self):
        """
        Return a list of all k-mer sequences that this node could
        theoretically have outgoing edges to.
        The list contains 4 strings, created by removing the first base
        and appending A, G, T, or C.
        """
        potential_to = []
        for base in ["A", "G", "T", "C"]:
            potential_to.append(self.k_mere_seq[1:] + base)
        return potential_to

    def get_edge_to_weight(self, other):
        """
        Return the weight of the incoming edge from DBGnode other
        to this node (0 if no such edge exists).
        other must be an instance of DBGnode.
        """
        return self.edges_to.get(other, 0)

    def get_edge_from_weight(self, other):
        """
        Return the weight of the outgoing edge from this node
        to DBGnode other (0 if no such edge exists).
        other must be an instance of DBGnode.
        """
        return self.edges_from.get(other, 0)
