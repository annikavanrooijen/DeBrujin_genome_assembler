# De Bruijn Graph–Based Genome Assembler

This project implements a simple genome assembler based on **de Bruijn graphs**.  
It reconstructs contigs from short sequencing reads by modeling k-mer overlaps as a directed graph and simplifying non-branching paths.

The implementation is intended for educational purposes and focuses on clarity and correctness rather than performance optimization.

---

## Project Context

This project was developed as part of the course  
**“Modelling and Visualization of Medical Data”**  
at **HTW Berlin** during the **Winter Semester 2025/26**.

**Instructor:** Ben Wulf  

---

## Overview

The assembler follows these main steps:

1. **Read parsing**  
   Sequencing reads are read from a FASTA file.

2. **k-mer extraction**  
   Each read is decomposed into overlapping k-mers of fixed length.

3. **de Bruijn graph construction**  
   - Nodes represent k-mers  
   - Directed edges represent (k-1)-overlaps between k-mers  

4. **Graph simplification**  
   Non-branching paths are merged into single nodes, producing contigs.

5. **Assembly statistics**  
   The final contigs can be exported in FASTA format, and assembly statistics
   such as the **N50** value are computed.

---

## Project Structure

src/

├── Read.py = FASTA parsing and k-mer extraction

├── DBGnode.py  = de Bruijn graph node implementation

├── DBGraph.py = de Bruijn graph construction and simplification

└── Main.py = Entry point and example usage

## Usage

### Build and simplify a graph

```bash
python3 Main.py
```

## Output
- **Graph summary:** Number of nodes (k-mers / contigs) and edges.
- **FASTA output:** The assembled contigs in FASTA format.
- **Assembly statistics:** Number of contigs and N50 value

## Assigning Sequences Using BLAST
The next important question is which organism the assembled sequence belongs to, potentially identifying a pathogen. To address this, the assembled genome sequence (or individual contigs) can be compared against databases of known sequences.

A widely used tool for this purpose is **BLAST (Basic Local Alignment Search Tool)**, provided by the***NCBI (National Center for Biotechnology Information)**.

BLAST can be accessed via:
https://blast.ncbi.nlm.nih.gov/Blast.cgi

## Notes
- All k-mers in the graph must have the same length.
- After simplification, the graph is considered final and cannot be extended with additional k-mers.
- The implementation does not perform error correction and assumes ideal input data.

## License
This project is released under the MIT License.

## Acknowledgements
Developed as part of coursework at HTW Berlin, 2025
Course: Modelling and Visualization of Medical Data
Instructor: Ben Wulf