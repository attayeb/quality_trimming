#

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from Bio.SeqRecord import SeqRecord
import random
import click
from os.path import join


def trim_fasta(reference_fasta, output_file_name, outputfolder, reference_name, p5="CCTACGGGNGGCWGCAG", p3 == "CCTACGGGNGGCWGCAG", readlength=310):
    """
    trim the fasta reference file using a primer and produces 2 fasta files R1, and R2.

    Parameters
    ----------
    reference_fasta : str
        The reference database fasta file name
    output_file_name : str
        The output file name followed by _R1 or _R2 accordingly
    outputfolder: str
        The output folder
    reference_name : str
        the reference read id
    p5 : str
        p5 primer
    p3 : str
        p3 primer
    """

    ref = SeqIO.parse(reference_fasta, "fasta")

    fasta_R1 = []
    fasta_R2 = []
    x = 0
    while True:
        try:
            x += 1
            seq = next(ref)
            p5i = nt_search(str(seq.seq), p5)
            p3i = nt_search(str(seq.reverse_complement().seq), p3)
            if (len(p5i) == 2) and (len(p3i) == 2):
                Merged = seq.seq[p5i[1]:len(seq)-p3i[1]]
                R1 = Merged[:readlength]
                R2 = Merged[-readlength:]

                fasta_R1.append(SeqRecord(R1, id="{}:R1".format(seq.id),
                                          description="{}:{}".format(reference_name, seq.id)))
                fasta_R2.append(SeqRecord(R2.reverse_complement(), id="{}:R2".format(
                    seq.id), description="{}:{}".format(reference_name, seq.id)))

        except:
            print("finish")
            break

    SeqIO.write(fasta_R1, join(
        outputfolder, "{}_R1.fasta".format(output_file_name)), "fasta")
    SeqIO.write(fasta_R2, join(
        outputfolder, "{}_R2.fasta".format(output_file_name)), "fasta")


@click.command()
@click.option('-r', '--reference-fasta', required=True,
              type=click.Path(exists=True))
@click.option('-o', '--output-file-name', required=True)
@click.option('--p5', default="CCTACGGGNGGCWGCAG")
@click.option('--p3', default="GACTACHVGGGTATCTAATCC")
@click.option('--reference-name', required=True)
@click.option('--readlength', default=310)
@click.option('--outputfolder', default=".")
def trim_fasta_click():
    trim_fasta(reference_fasta, output_file_name, outputfolder, reference_name, p5, p3, readlength)


if __name__ == "__main__":
    trim_fasta_click()
