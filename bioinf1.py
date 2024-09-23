"""
Project Goal:
~~~~~~~~~~~~~
To develop a Python program using libraries such as NumPy and Biopython to analyze DNA sequences, including:

Nucleotide composition calculation
Reverse complement calculation
GC content calculation
Identification of Open Reading Frames (ORFs)


Project Plan:
~~~~~~~~~~~~~
Input: The user provides a DNA sequence (as a string or FASTA file).
Nucleotide Composition: The program counts occurrences of each nucleotide (A, T, G, C).
Reverse Complement: Calculate the reverse complement of the DNA sequence.
GC Content: Calculate the percentage of G and C nucleotides in the sequence.
Find Open Reading Frames (ORFs): Identify all ORFs that start with a start codon (ATG) and end with a stop codon (TAA, TAG, TGA).
"""


import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO

# Function to calculate nucleotide composition
def nucleotide_composition(sequence):
    """
    Nucleotide Composition:
    Counts the number of each nucleotide (A, T, G, C) in the given sequence.
    """
    comp = {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'G': sequence.count('G'),
        'C': sequence.count('C')
    }
    return comp

# Function to calculate reverse complement of a sequence
def reverse_complement(sequence):
    seq = Seq(sequence)
    return str(seq.reverse_complement())

# Function to calculate GC content
def gc_content(sequence):
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    return (g_count + c_count) / len(sequence) * 100

# Function to find open reading frames (ORFs)
def find_orfs(sequence):
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    
    orfs = []
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        print ('', codon)
        if codon == start_codon:
            for j in range(i + 3, len(sequence) - 2, 3):
                stop_codon = sequence[j:j + 3]
                if stop_codon in stop_codons:
                    orfs.append(sequence[i:j + 3])
                    break
    return orfs

# Function to analyze DNA sequence
def analyze_dna(sequence):
    print("Nucleotide Composition:", nucleotide_composition(sequence))
    print("Reverse Complement:", reverse_complement(sequence))
    print("GC Content: {:.2f}%".format(gc_content(sequence)))
    orfs = find_orfs(sequence)
    print(f"Number of ORFs found: {len(orfs)}")
    for idx, orf in enumerate(orfs):
        print(f"ORF {idx + 1}: {orf}")

# Example sequence input (you can replace this with file reading if needed)
example_sequence = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

# Analyze the example sequence
analyze_dna(example_sequence)



"""
Reverse Complement:
A DNA sequence is complementary when A pairs with T and G pairs with C. The reverse complement is the sequence read in the opposite direction.
"""

"""
GC Content:
GC content refers to the percentage of guanine (G) and cytosine (C) bases. This metric is important as it relates to the stability of the DNA molecule.
"""

"""
Open Reading Frames (ORFs):
ORFs are sections of the DNA that can potentially code for proteins. They start with an ATG codon and end with one of the stop codons (TAA, TAG, or TGA).
"""



"""
4. Input and Output:
Input:
You can either pass the sequence as a string or read from a FASTA file using SeqIO from Biopython.

Output:
The program will output:
"""


"""
The nucleotide composition
The reverse complement of the sequence
The GC content (as a percentage)
The number and details of the ORFs in the sequence
"""


"""
5. Additional Enhancements:
I can further improve the project by:
Allowing the user to input sequences in FASTA format.
Identifying ORFs on both strands of the DNA.
Plotting GC content using libraries like Matplotlib for graphical representation.
Analyzing larger datasets or even incorporating machine learning to predict gene sequences.
"""