"""
Goal: Learn to Read & Write a Bio Sequence File
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Generate a Synthetic DNA Sequence

synthetic_seq0 = Seq("ATGATCATCGATTAACTGGGTATGGGTATGGCCCTATCTATCTATCTATCTATCAATATATATATATATATGA")
synthetic_seq1 = Seq("ATGTTTTTTTTTTTTTTTTTAATGGGTATGGCCCTATTATCTATCTATCAATATATATATATATATAA")
synthetic_seq2 = Seq("ATGATCATCGATTAACCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGATATATATAG")
synthetic_seq3 = Seq("ATGAAAAAAGCTAAAAAAGCTAAAAAAAAGCTAAAAAAAAAAAAAAAAAAGCTTTTTTTAAAAAAAAAAAACAATATATATATATATATGA")

record1 = SeqRecord (synthetic_seq0, id ="synthetic_gene0", description="Synthetic DNA Sequence 0")
record2 = SeqRecord (synthetic_seq1, id ="synthetic_gene1", description="Synthetic DNA Sequence 1")
record3 = SeqRecord (synthetic_seq2, id ="synthetic_gene2", description="Synthetic DNA Sequence 2") 
record4 = SeqRecord (synthetic_seq3, id ="synthetic_gene3", description="Synthetic DNA Sequence 3")


# Save this to FASTA File!
with open ("synthetic_gene.fasta", "w") as output_handle:
     SeqIO.write(record1, output_handle, "fasta")
     SeqIO.write(record2, output_handle, "fasta")
     SeqIO.write(record3, output_handle, "fasta")     
     SeqIO.write(record4, output_handle, "fasta")
# ####################################################################
# note: WITH statement handles exceptions AND closes file once DONE!
# ####################################################################

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Goodness of the with statement :) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The 'with' statement is popularly used with file streams, 
as shown above 
and with Locks, sockets, subprocesses and telnets etc.

'with' statement makes the code compact and much more readable. 
'with' statement helps avoiding bugs and leaks by ensuring 
 that a resource is properly released when the code using 
 the resource is completely executed.

"""


with open("synthetic_gene.fasta", "r") as input_handle:
    for record in SeqIO.parse (input_handle, "fasta"):
        print ("ID: ", (record.id))
        print ("Sequence Length:", len(record))
        print ("Description: ", record.description, "\n")
                             

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
How to Read a DNA Sequence :) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SeqIO.read()
Turn a sequence file into a single SeqRecord.

Arguments:
~~~~~~~~~~
handle - handle to the file, or the filename as a string
(note older versions of Biopython only took a handle).
format - string describing the file format.
alphabet - no longer used, should be None.
This function is for use parsing sequence files containing exactly one record. 
"""


"""
SeqIO.parse()
Turn a sequence file into an iterator returning SeqRecords.

Arguments:
~~~~~~~~~~
handle - handle to the file, or the filename as a string
(note older versions of Biopython only took a handle).
format - lower case string describing the file format.
alphabet - no longer used, should be None.
Typical usage, opening a file to read in, and looping over the record(s):

>>> from Bio import SeqIO
>>> filename = "Fasta/sweetpea.nu"
>>> for record in SeqIO.parse(filename, "fasta"):
...    print("ID %s" % record.id)
...    print("Sequence length %i" % len(record))
ID gi|3176602|gb|U78617.1|LOU78617
Sequence length 309
For lazy-loading file formats such as twobit, for which the file contents is read on demand only, ensure that the file remains open while extracting sequence data.

If you have a string 'data' containing the file contents, you must first turn this into a handle in order to parse it:

>>> data = ">Alpha\nACCGGATGTA\n>Beta\nAGGCTCGGTTA\n"
>>> from Bio import SeqIO
>>> from io import StringIO
>>> for record in SeqIO.parse(StringIO(data), "fasta"):
...     print("%s %s" % (record.id, record.seq))
Alpha ACCGGATGTA
Beta AGGCTCGGTTA
Use the Bio.SeqIO.read(...) function when you expect a single record only.
"""