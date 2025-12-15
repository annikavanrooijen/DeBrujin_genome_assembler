import Read, DBGraph

def read_fasta(readfile: str):
    '''Read a FASTA file and parse it into a list of Read objects'''
    reads = []
    with open(readfile, "r") as infile:
        header = None
        sequence = ""
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    reads.append(Read([header, sequence]))
                header = line
                sequence = ""
            else:
                sequence += line
        if header is not None:
            reads.append(Read([header, sequence]))
    return reads

def build_graph(filename, kmersize):
    """ Nun können Sie endlich den Graph aus einer FASTA-Datei erstellen: 
    Implementieren Sie eine top-level-Methode build_graph(filename, kmersize) 
    die aus den Reads in der FASTA-Datei im Pfad filename einen Graph mit der k-mer-Länge kmersize erstellt und zurückgibt."""
    graph = DBGraph()
    for read in read_fasta(filename):
        graph.add_kmers(read.get_kmers(kmersize))  
    return graph 