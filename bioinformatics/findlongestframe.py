# Author - Bill Winnett
# email - bwinnett12@gmail.com
# A simple command line program to find lengths of reading frames. Default is using table 4
import sys
from Bio import SeqIO
from Bio.Seq import Seq

file_in = sys.argv[1]
# file_in = "cox1_segmented/AF538053.1_new.fasta"
for record in SeqIO.parse(file_in, "fasta"):

    for off in range(0, 3):
        for direc in ["to", "back"]:
            sequence = str(record.seq) if direc == "to" else str(record.seq)[::-1]

            translated = Seq(sequence[off:(len(sequence) - off) - ((len(sequence) - off) % 3)]).\
                translate(table=4).split("*")
            # translated = [x for x in translated if "*" not in x]
            longest = max(list(filter(None, translated)), key=len)

            print("Frame " + str(off) + " - " + direc + ": " + str(len(longest)))
            print(longest)
            r = 2



        # frame2_to = Seq(seq[off:(len(seq) - off) - ((len(seq) - off) % 3)]).translate(to_stop=True, table=4)

    # reverse = seq_to[::-1]
    #
    # frame1_rev = Seq(reverse[:(len(reverse)) - ((len(reverse)) % 3)]).translate(to_stop=True, table=4)
    # frame2_rev = Seq(reverse[1:(len(reverse) - 1) - ((len(reverse) - 1) % 3)]).translate(to_stop=True, table=4)
    # frame3_rev = Seq(reverse[2:(len(reverse) - 2) - ((len(reverse) - 2) % 3)]).translate(to_stop=True, table=4)

    # print("Frame 1 - to: " + str(len(frame1_to)))
    # print("Frame 2 - to: " + str(len(frame2_to)))
    # print("Frame 3 - to: " + str(len(frame3_to)))
    # print("Frame 1 - rev: " + str(len(frame1_rev)))
    # print("Frame 2 - rev: " + str(len(frame2_rev)))
    # print("Frame 3 - rev: " + str(len(frame3_rev)))
