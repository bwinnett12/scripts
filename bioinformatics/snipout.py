# Author - Bill Winnett
# email - bwinnett12@gmail.com
# Simple command line to extract a fasta from an existing fasta but only in a certain segment

import sys
from Bio import SeqIO

file_in = sys.argv[1]
start_loc = int(sys.argv[2]) - 1
end_loc = int(sys.argv[3]) - 1

name = sys.argv[4] if len(sys.argv) > 3 else "new"

for record in SeqIO.parse(file_in, "fasta"):
    sequence = str(record.seq[start_loc:end_loc])

    outfile = open(record.name + "_" + name + ".fasta", "w")
    outfile.write("> " + record.id + "\n")

    for n in range(0, len(sequence), 75):
        outfile.write(sequence[n:n+75] + "\n")
    outfile.write("\n")


