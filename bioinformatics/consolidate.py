
__author__ = "Bill Winnett"
# For use on ctenophore evolution project
# https://github.com/bwinnett12/ctenophora_genome_analysis
import os, glob
from Bio import SeqIO


folder_a = "./fetch_data/faa/"
folder_b = "./mitos2_data/"
out_folder = "./consolidated/"

# adds all records from fetch data
def add_from_fetch(in_folder, out_folder):

    for file in glob.glob(out_folder + "*.faa"):
        open(file, "w")

    for file in glob.glob(in_folder + "*.faa"):  # Gets each file
        for record in SeqIO.parse(file, "fasta"):  # gets the records from each one

            description = record.description.split(":")
            species = description[-1]
            gene = description[1]
            seq = str(record.seq)

            if "ND" in gene:
                gene = gene.replace("ND", "NAD")

            # Makes a file if it doesn't exist
            if not os.path.isfile(out_folder + species + ".faa"):
                species_file = open(out_folder + species + ".faa", "x")
            species_file = open(out_folder + species + ".faa", "a")

            species_file.write(">" + gene + ":FETCH:" + species + "\n")
            for pos in range(0, len(seq), 75):
                species_file.write(seq[pos:pos+75] + "\n")
            species_file.write("\n")


# Adds information from mitos and appends to folders
def add_from_mitos(in_folder, out_folder):
    for file in glob.glob(in_folder + "/*/*.faa"):  # Gets each file in each mitos folder
        for record in SeqIO.parse(file, "fasta"):  # Gets each record

            description = record.description
            gene = description.split(";")[-1].upper().lstrip()
            species = file.split("/")[-2]
            seq = str(record.seq)

            # Makes file if it doesn't exist
            if not os.path.isfile(out_folder + species + ".faa"):
                species_file = open(out_folder + species + ".faa", "x")
            species_file = open(out_folder + species + ".faa", "a")

            # Generates each fasta entry
            species_file.write(">" + gene + ":MITOS:" + species + "\n")
            for pos in range(0, len(seq), 75):
                species_file.write(seq[pos:pos+75] + "\n")
            species_file.write("\n")


# Sorts each sequence alphabetically then Fetch, mitos
def organize(out_folder):
    def SortTuple(tup):

        # Getting the length of list
        # of tuples
        n = len(tup)

        for i in range(n):
            for j in range(n - i - 1):

                if tup[j][0] > tup[j + 1][0]:
                    tup[j], tup[j + 1] = tup[j + 1], tup[j]

        return sorted(tup, key = lambda x: x[0])

    for file in glob.glob(out_folder + "*.faa"):

        records = []
        for record in SeqIO.parse(file, "fasta"):
            records.append((record.description, record))

        # for record in records:
        records = SortTuple(records)
        file = open(file, "w")

        for record in records:
            file.write(">" + record[1].description + "\n")
            for pos in range(0, len(str(record[1].seq)), 75):
                file.write(str(record[1].seq)[pos:pos+75] + "\n")
            file.write("\n")









list = ["a", "b", "c", "d"]
print(list[list.index("b")::])
os.chdir("../")
add_from_fetch(folder_a, out_folder)
add_from_mitos(folder_b, out_folder)
organize(out_folder)
