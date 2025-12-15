from DBGraph import DBGraph
from Read import Read

def read_fasta(readfile: str):
        '''
        Read a FASTA file and parse it into a list of Read objects.

        Args:
            readfile (str): Path to the FASTA file.

        Returns:
            list[Read]: List of parsed Read objects.
        '''
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
        '''Build and return a de Bruijn graph from a FASTA file.

        This top-level helper reads all sequences from `filename`, extracts k-mers
        of length `kmersize` for each read, and inserts them into a new DBGraph.

        Args:
            filename (str): Path to the FASTA file containing reads.
            kmersize (int): Length of k-mers to generate from each read.

        Returns:
            DBGraph: The constructed de Bruijn graph.
        '''
        graph = DBGraph()
        for read in read_fasta(filename):
            graph.add_kmers(read.get_kmers(kmersize))  
        return graph 

if __name__ == "__main__":
    kmersize= 15
    print(f"Kmer-size: {kmersize}")
    dbg = build_graph("data/virus_perfectreads.fasta", kmersize)
    print(dbg)
    dbg.simplify()
    print(dbg)
    print(f"N50: {dbg.get_N50()}")
    print(f"Number of contigs: {dbg.get_numContigs()}")
    print(dbg.get_FASTA())

