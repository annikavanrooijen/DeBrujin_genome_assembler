from DBGnode import DBGnode
import itertools

class DBGraph:

    '''
    Represents a de Bruijn graph constructed from sequencing reads.

    Nodes correspond to k-mers (or merged contigs after simplification),
    and directed edges represent (k-1)-overlaps between k-mers.
    '''
    
    def __init__(self):
        '''
        Initialize an empty de Bruijn graph.

        Attributes:
            nodes (dict): Maps k-mer/contig sequences to DBGnode objects.
            kmerlen (int or None): Length of k-mers used in the graph.
            simplified (bool): Indicates whether the graph has been simplified.
        '''
        self.nodes = {}
        self.kmerlen = None
        self.simplified = False

    def add_kmers(self, kmers):
        '''Add k-mers (and their multiplicities) to the de Bruijn graph.

        The input `kmers` is expected to be a dict mapping k-mer strings to counts
        (as returned by Read.get_kmers()). For each k-mer, this method:
        - validates that its length matches the graph's k-mer length,
        - creates a new DBGnode if the k-mer is not present yet,
        - adds directed edges to/from already existing compatible k-mers
            (based on (k-1)-overlaps via get_potential_to / get_potential_from).

        Notes:
        - If the graph was already simplified, no further k-mers are added.
        - If `kmers` is empty, the method returns immediately.

        Raises:
            ValueError: if k-mers with incompatible lengths are added.
        '''
        if self.simplified:
            return
        if len(kmers) == 0:
            return
        if self.kmerlen is None:
            self.kmerlen = len(next(iter(kmers.keys())))
        for kmer_s in kmers.keys():
            if len(kmer_s) != self.kmerlen:
                raise ValueError("Incompatible k-mer lengths: " + str(self.kmersize) + " and " + str(len(kmer_s)))
            if kmer_s not in self.nodes.keys():
                self.nodes[kmer_s] = DBGnode(kmer_s)
            kmer = self.nodes[kmer_s]
            for pto in kmer.get_potential_to():
                if pto in self.nodes.keys():
                    self.nodes[pto].add_edge_from(kmer)
                    kmer.add_edge_to(self.nodes[pto])
            for pfrom in kmer.get_potential_from():
                if pfrom in self.nodes.keys():
                    self.nodes[pfrom].add_edge_to(kmer)
                    kmer.add_edge_from(self.nodes[pfrom])

    def count_edges(self):
        '''
        Count the total number of edges in the graph.
        Each edge is stored twice (incoming and outgoing),
        but is counted only once in this method.
        '''
        counter = 0
        for node in self.nodes.values():
             counter += len(node.edges_to)
        return counter
    
    def count_nodes(self):
        '''Return the total number of nodes in the graph,
        corresponding to the number of unique k-mers.
        '''
        return len(self.nodes)

    def __str__(self):
        '''Return a concise textual summary of the graph,
        including k-mer size, number of nodes, and number of edges.
        '''
        return "DBG(" + str(self.kmerlen) + ") with " + str(self.count_nodes()) + " nodes and " + str(
            self.count_edges()) + " edges"

    def simplify(self):
        '''Simplify the de Bruijn graph by merging non-branching paths into contigs.

        The method repeatedly extends each node to the right (using extend_next())
        as long as it can be uniquely extended, thereby collapsing linear chains
        into single nodes with longer sequences (contigs).

        After merging:
        - removed nodes are deleted from `self.nodes`,
        - the node dictionary is rebuilt so that keys match the updated node.seq,
        - `self.simplified` is set to True to prevent further k-mer insertions.
        '''
        for node in list(self.nodes.values()):
            while True:
                next_node = node.extend_next()
                if next_node is None:
                    break
                self.nodes.pop(next_node.seq, None)

        new_nodes = {}
        for node in self.nodes.values():
            new_nodes[node.seq] = node
        self.nodes = new_nodes

        self.simplified = True

    def get_FASTA(self):
        '''
        Export the current graph nodes as a FASTA-formatted string.

        Each node sequence is written as a contig entry:
            >contig1
            ACTG...
            >contig2
            ...

        Returns:
            str: FASTA representation of all node sequences (contigs).
        '''
        return "\n".join(f">contig{i+1}\n{node.seq}" for i, node in enumerate(self.nodes.values()))
    
    def get_numContigs(self):
        '''
        Returns:
        int or None: Number of contigs if the graph is simplified, otherwise None.
        '''
        if self.simplified:
            return len(self.nodes)
        return None
    
    def get_N50(self):
        '''
        Compute the N50 statistic of the contig set (after simplification).

        N50 is defined as the contig length L such that at least 50% of the total
        assembled length is contained in contigs of length >= L. This implementation:
        - computes all contig lengths,
        - sorts them in ascending order,
        - builds a cumulative sum until it exceeds half of the total length,
        - returns the contig length at that point.

        Returns:
            int: The N50 contig length.
        '''
        total_length = sum(len(k) for k in self.nodes)
        lengths = sorted(len(k) for 
                         k in self.nodes)
        cumSum = list(itertools.accumulate(lengths))
        half = total_length/2
        for cum, L in zip(cumSum, lengths):
            if cum > half:
                return L