# Represents a single sequencing read from a FASTA file,  storing its identifier and nucleotide sequence.
class Read:

    def __init__(self, lines):
        '''Initialize the read name (without '>') and concatenate sequence lines'''
        self.name = lines[0].strip().replace(">", "")
        self.bases = "".join(lines[1:]).replace(" ", "").strip()

    def get_kmers(self, kmersize):
        '''Extract all k-mers of length kmersize from the sequence
        and count their occurrences'''
        kmers = {}
        for i in range(len(self.bases) - kmersize + 1):
            kmer = self.bases[i:i + kmersize]
            if kmer in kmers:
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
        return kmers

    def __str__(self):
        '''Return a human-readable string representation of the read
        (truncate long sequences for readability)'''
        if len(self.bases) > 20:
            return f"{self.name}: {self.bases[0:20]}..."
        return f"{self.name}: {self.bases}"

    def __repr__(self):
        '''Define how the object is represented in the console and during debugging'''
        return self.__str__()

    def __eq__(self, other): 
        '''Define equality: two Read objects are equal if both their names and sequences are identical'''
        if not isinstance(other, Read):
            return False
        return self.name == other.name and self.bases == other.bases


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
