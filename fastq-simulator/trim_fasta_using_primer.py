# 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from Bio.SeqRecord import SeqRecord
import random
import click


@click.command()
@click.option('-r', '--reference-fasta', required=True,
              type=click.Path(exists=True))
@click.option('-o', '--output-file-name', required=True)
@click.option('--p5', default="CCTACGGGNGGCWGCAG")
@click.option('--p3', default="GACTACHVGGGTATCTAATCC")
@click.option('--reference-name', required=True)
def main(reference_fasta, output_file_name, reference_name, p5, p3):
    """
    This script extracts reads from reference fasta file depending on two primers.

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
                R1 = Merged[:310]
                R2 = Merged[-310:]

                fasta_R1.append(SeqRecord(R1, id="{}:R1".format(seq.id),
                                          description="{}:{}".format(reference_name, seq.id)))
                fasta_R2.append(SeqRecord(R2.reverse_complement(), id="{}:R2".format(seq.id),
                                          description="{}:{}".format(reference_name, seq.id)))

        except:
            print("finish")
            break

    SeqIO.write(fasta_R1, "{}_R1.fasta".format(output_file_name), "fasta")
    SeqIO.write(fasta_R2, "{}_R2.fasta".format(output_file_name), "fasta")


if __name__ == "__main__":
    main()
