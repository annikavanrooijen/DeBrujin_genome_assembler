import DBGnode

class DBGraph:
    def __init__(self):
        '''Initialize an empty de Bruijn graph and its internal data structures.'''
        
        self.nodes = {}
        self.kmersize = None

    def add_kmers(self, kmers):
        '''Insert k-mers into the graph and connect them according
        to their (k-1)-overlaps.

        Each unique k-mer is represented as a node. Directed edges
        are added between nodes if the suffix of one k-mer overlaps
        with the prefix of another k-mer.

        The method ensures that all k-mers added to the graph have
        the same length. Adding a k-mer of a different length
        raises a ValueError.'''
        
        self.kmersize = len(next(iter(kmers)))

        for kmer in kmers:
            if len(kmer) != self.kmersize:
                 raise ValueError(f"A k-mer of another length can not be added")
            
            if kmer in self.nodes:
                node = self.nodes[kmer]
            else:
                node = DBGnode(kmer) 
                self.nodes[kmer] = node

            potential_from = node.get_potential_from()
            for x in potential_from:
                if x in self.nodes:
                        self.nodes[kmer].add_edge_from(self.nodes[x])
                        self.nodes[x].add_edge_to(self.nodes[kmer])

            potential_to = node.get_potential_to()
            for y in potential_to:
                if y in self.nodes:
                        self.nodes[kmer].add_edge_to(self.nodes[y])
                        self.nodes[y].add_edge_from(self.nodes[kmer])


    def count_edges(self):
        '''
        Count the total number of edges in the graph.
        Each edge is stored twice (incoming and outgoing),
        but is counted only once in this method.'''
        
        counter = 0
        for node in self.nodes.values():
             counter += len(node.eto)
        return counter
    
    def count_nodes(self):
        '''Return the total number of nodes in the graph,
        corresponding to the number of unique k-mers.'''

        return len(self.nodes)

    def __str__(self):
        '''Return a concise textual summary of the graph,
        including k-mer size, number of nodes, and number of edges.'''
        
        return f" x = {self.kmersize}, y = {self.count_nodes()}, z = {self.count_edges()}"
