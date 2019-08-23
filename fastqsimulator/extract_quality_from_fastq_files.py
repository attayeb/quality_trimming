# AUTHOR: Attayeb Mohsen
# DESCRIPTION: This script extracts phred quality from fastq files.

from Bio import SeqIO
import os


def extract_quality_scores(folder, R1q_filename="R1Quality.txt", R2q_filename="R2Quality.txt"):
    """
    extracts quality scores from fastq files and save them in two text files

    Parameters
    ----------
    folder : str
        The folder which contains the fastq files
    R1q_filename : str, optional
        R1 quality scores file name, by default "R1Quality.txt"
    R2q_filename : str, optional
        R2 quality scores file name, by default "R2Quality.txt"
    
    Examples
    --------
    >>> extract_quality_scores("samples/", "R1.qual", "R2.qual")
    """
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

    with open(R1q_filename, 'w') as f:
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

    with open(R2q_filename, 'w') as f:
        for item in R2:
            f.write(item)


if __name__ == "__main__":
    extract_quality_scores("samples/", "r1_q.txt", "r2q.txt")
