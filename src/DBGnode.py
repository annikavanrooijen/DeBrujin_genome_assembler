class DBGnode:
    '''
    Represents a node in a de Bruijn graph.

    Each node corresponds to a k-mer sequence and stores weighted directed
    edges to predecessor and successor nodes. During graph simplification,
    nodes may be merged to form longer contig sequences.
    '''

    def __init__(self, seq):
        '''
        Initialize a de Bruijn graph node with its k-mer sequence.

        Args:
            seq (str): The k-mer sequence represented by this node.
        '''
        self.seq = seq
        # Dictionaries mapping DBGnode -> edge weight
        # edges_to   : incoming edges
        # edges_from : outgoing edges
        self.edges_to = {}
        self.edges_from = {}
        
    def add_edge_to(self, eto):
        '''
        Add an incoming edge from node `eto` to this node.

        If the edge already exists, its weight is increased by 1.

        Args:
            eto (DBGnode): The source node of the edge.
        '''
        self.edges_to[eto] = self.edges_to.get(eto, 0) + 1

    def add_edge_from(self, efrom):
        '''
        Add an outgoing edge from this node to node `efrom`.

        If the edge already exists, its weight is increased by 1.

        Args:
            efrom (DBGnode): The destination node of the edge.
        '''
        self.edges_from[efrom] = self.edges_from.get(efrom, 0) + 1

    def get_potential_from(self):
        '''
        Compute all possible predecessor k-mers for this node.

        These are the k-mers that share a (k-1)-suffix overlap with this node
        and could therefore have an incoming edge.

        Returns:
            list[str]: Four possible predecessor k-mer strings.
        '''
        potential_from = []
        for base in ["A", "G", "T", "C"]:
            potential_from.append(base + self.seq[:-1])
        return potential_from

    def get_potential_to(self):
        '''
        Compute all possible successor k-mers for this node.

        These are the k-mers that share a (k-1)-prefix overlap with this node
        and could therefore have an outgoing edge.

        Returns:
            list[str]: Four possible successor k-mer strings.
        '''
        potential_to = []
        for base in ["A", "G", "T", "C"]:
            potential_to.append(self.seq[1:] + base)
        return potential_to

    def get_edge_to_weight(self, other):
        '''
        Get the weight of the incoming edge from `other` to this node.

        Args:
            other (DBGnode): Source node.

        Returns:
            int: Edge weight (0 if no edge exists).
        '''
        return self.edges_to.get(other, 0)

    def get_edge_from_weight(self, other):
        '''
        Get the weight of the outgoing edge from this node to `other`.

        Args:
            other (DBGnode): Destination node.

        Returns:
            int: Edge weight (0 if no edge exists).
        '''
        return self.edges_from.get(other, 0)
    
    def extend_next(self):
        '''
        Extend (merge) this node to the right by combining it with its unique successor.

        If `can_extend_next()` is True, this method merges the successor node into
        the current node by computing the maximal overlap between:
        - the suffix of `self.seq` and
        - the prefix of `next_node.seq`.

        After merging:
        - `self.seq` is extended with the non-overlapping suffix from the successor,
        - the edge from this node to the successor is removed,
        - all outgoing edges of the successor are transferred to this node,
        - and all connections of the successor are cleared.

        Returns:
            DBGnode or None:
            The successor node that was merged (so it can later be removed from
            the graph), or None if no extension is possible.'''
        if not self.can_extend_next():
            return None
        
        (next_node, _ ), = self.edges_to.items()
        max_overlap = min(len(self.seq), len(next_node.seq))
        overlap = 0
        for i in range(max_overlap, 0, -1):
            if self.seq.endswith(next_node.seq[:i]):
                overlap = i
                break

        self.seq = self.seq + next_node.seq[overlap:]

        self.edges_to.pop(next_node)

        for eto_next_node, weight in next_node.edges_to.items():
            self.edges_to[eto_next_node] = self.edges_to.get(eto_next_node,0) + weight
            eto_next_node.edges_from.pop(next_node, 0)
            eto_next_node.edges_from[self] = eto_next_node.edges_from.get(self , 0) + weight

        next_node.edges_to.clear()
        next_node.edges_from.clear()
        return  next_node

    def extend_prev(self):
        '''
        Extend (merge) this node to the left by combining it with its unique predecessor.

        If `can_extend_prev()` is True, this method merges the predecessor node into
        the current node by computing the maximal overlap between:
        - the prefix of `self.seq` and
        - the suffix of `previous_node.seq`.

        After merging:
        - `self.seq` is extended with the missing prefix from the predecessor,
        - the edge from the predecessor to this node is removed,
        - all incoming edges of the predecessor are transferred to this node,
        - and all connections of the predecessor are cleared.

        Returns:
            DBGnode or None:
            The predecessor node that was merged (so it can later be removed from
            the graph), or None if no extension is possible.
        '''
        if not self.can_extend_prev():
            return None
        
        (previous_node, _ ), = self.edges_from.items()
        max_overlap = min(len(self.seq), len(previous_node.seq))
        overlap = 0
        for i in range(0, max_overlap):
            if self.seq.startswith(previous_node.seq[i:]):
                overlap = i
                break
        
        self.seq = previous_node.seq[:overlap] + self.seq 
        self.efroedges_fromm.pop(previous_node)

        for efrom_previous_node, weight in previous_node.edges_from.items():
            self.efredges_fromom[efrom_previous_node] = self.efedges_fromrom.get(efrom_previous_node,0) + weight
            efrom_previous_node.edges_to.pop(previous_node, 0)
            efrom_previous_node.edges_to[self] = efrom_previous_node.edges_to.get(self , 0) + weight

        previous_node.edges_to.clear()
        previous_node.edges_from.clear()
        return  previous_node

    def can_extend_next(self):
        '''
        Check whether this node can be safely extended to its successor.

        Extension is possible only if:
        - the node has exactly one outgoing edge,
        - the successor node is not the node itself,
        - the successor has exactly one incoming edge,
        - and this node is the only predecessor of the successor.

        Returns:
            bool:
            True if the node can be extended to the next node, otherwise False.'''
        if len(self.edges_to) != 1:
            return False
        
        (next_node, _ ), = self.edges_to.items()
        if next_node is self:
            return False
        
        if len(next_node.edges_from) != 1:
            return False
        
        (efrom_next_node, _ ), = next_node.edges_from.items()
        if efrom_next_node is not self:
            return False 
        
        return True 

    def can_extend_prev(self):
        '''
        Check whether this node can be safely extended to its predecessor.

        Extension is possible only if:
        - the node has exactly one incoming edge,
        - the predecessor node is not the node itself,
        - the predecessor has exactly one outgoing edge,
        - and this node is the only successor of the predecessor.

        Returns:
            bool:
            True if the node can be extended to the previous node, otherwise False.'''
        if len(self.edges_from) != 1:
            return False 
        
        (previous_node, _ ), = self.edges_from.items()
        if previous_node is self:
            return False
        
        if len(previous_node.edges_to) != 1:
            return False
        
        (eto_previous_node,_ ), = previous_node.edges_to.items()
        if eto_previous_node is not self:
            return False
            
        return True