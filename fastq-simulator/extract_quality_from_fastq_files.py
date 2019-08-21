# AUTHOR: Attayeb Mohsen
# DESCRIPTION: This script extracts phred quality from fastq files.

from Bio import SeqIO
import os
folder = "fastq_files_for_quality_scores/"
files = os.listdir(folder)

R1 = []
R1files = [x for x in files if "_R1_" in x]
for file in R1files:
    print(file)
    print(len(R1))
    f = open(folder+"{}".format(file))
    lines = f.readlines()
    x = 0
    while True:
        try:
            line = lines[x]
            x = x+1
            if x % 4 == 0:
                if line in R1:
                    print(".", end="")
                    pass
                else:
                    R1.add(line)
        except:
            break

with open('R1_quality.txt', 'w') as f:
    for item in R1:
        f.write(item)


R2 = []
R2files = [x for x in files if "_R2_" in x]
for file in R2files:
    print(file)
    print(len(R2))
    f = open(folder + "{}".format(file))
    lines = f.readlines()
    x = 0
    while True:
        try:
            line = lines[x]
            x = x + 1
            if x % 4 == 0:
                if line in R2:
                    print(".", end="")
                    pass
                else:
                    R2.add(line)
        except:
            break


with open("R2_quality.txt", 'w') as f:
    for item in R2:
        f.write(item)