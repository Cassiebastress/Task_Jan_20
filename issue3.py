# Demonstrating use of Seq in Biopython.
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess

# Exploring Biopython Seq
'''
my_seq = Seq("AGTACACTGGT")
print(my_seq)
complement = my_seq.complement()
print(complement)
reverse_complement = my_seq.reverse_complement()
print(reverse_complement)

'''

# Parsing sequence fasta file
'''
for seq_record in SeqIO.parse("FOXP2.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
'''

# Other functions
'''
print(my_seq.count("A"))
print(len(my_seq))
print(my_seq[4:7])
'''

# Creating a SeqRecord
'''
simple_seq = Seq("GATC")
simple_seq_r = SeqRecord(simple_seq)
simple_seq_r.id = "AC12345"
simple_seq_r.description = "Made up"
simple_seq_r.name = "something"
simple_seq_r.annotations["evidence"] = "None"
print(simple_seq_r)
'''

# Reading from a file
'''
record = SeqIO.read("FOXP2.fasta", "fasta") # need a file with only one record!
record_seq = record.seq
record_id = record.id
record_name = record.name
record_description = record.description
'''

# FUNCTIONS:

# Function 1:
# Take input file and output file names
# Use Mafft to align sequences
def mafft_align(input_fasta, output_fasta): 
    with open(output_fasta, "w") as outfile:
        subprocess.run(["mafft", "--auto", input_fasta], stdout=outfile, stderr=subprocess.DEVNULL)

mafft_align("input.fasta", "output.fasta")

# Function 2:
# Take input file name
# Return list of Seq objects

# Function 3:
# Take a list of Seq objects
# Return a list of Seq containing the allignments

# Function 4:
# Unsure of what this function should do