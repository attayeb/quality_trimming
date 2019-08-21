# pylint: disable=no-value-for-parameter
"""
"""
# coding: utf-8

# AUTHOR: Attayeb Mohsen
# DESCRIPTION:

import random
from multiprocessing.dummy import Pool
import click
import numpy as np
from numpy.random import choice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# example:
# python simulator.py --reference gg_R1.fasta gg_R2.fasta --quality R1_quality.txt R2_quality.txt
# --sample-names-pattern stage2- --number-of-samples 15 --sample-size-range 10000 50000

def calcualte_p(quality_score):
    """
    calculate the p value for phred quality score

    Parameters
    ----------
    quality_score : int
        The quality score

    Returns
    -------
    float
        Probability of the error
    """
    return 10 ** (-quality_score/10)


def suggest_nucleotide(correct, quality_score):
    """
    Suggest a nucleotide depending on the quality score

    Parameters
    ----------
    correct : str
        The correct nucleotide ['A', 'C', 'T', 'G']
    quality_score : int
        Phred quality

    Returns
    -------
    str ['A', 'C', 'T', 'G']
        The suggested nucleotide depending on the quality score.
    """
    p = 1 - calcualte_p(quality_score)
    other_p = (1 - p) / 3
    p_list = [p, other_p, other_p, other_p]
    if correct in "ACTG":
        nucs = [n for n in "ACTG".replace(correct, "")]
        choices = [correct]
        choices.extend(nucs)
    else:
        p_list = [0.25, 0.25, 0.25, 0.25]
        choices = ['A', 'C', 'T', 'G']

    return choice(choices, p=p_list)


def quality_based_modification(original_seq, qual):
    """
    Modify a sequence string using a list of quality scores

    Parameters
    ----------
    original_seq : str
        the oriqinal sequence
    qual : list of ints
        the list of quality scores

    Returns
    -------
    str
        modified sequence
    """
    res_seq = ""
    for n in range(len(original_seq)):
        res_seq += suggest_nucleotide(original_seq[n], qual[n])
    return res_seq


def phred_to_values(string):
    """
    convert phred quality character to phred value

    Parameters
    ----------
    string : str
        text presentation of phred quality

    Returns
    -------
    list of int
        list of phred quality values
    """
    return [ord(x)-33 for x in string.strip()]


def get_quality(R1, R2, number):
    """
    read quality scores from the quality pool

    Parameters
    ----------
    R1 : str
        R1 quality scores file name
    R2 : str
        R2 quality scores file name
    number : int
        number of scores extracted

    Returns
    -------
    two numpy arrays
        quality scores of R1 and R2
    """
    r1f = open(R1)
    seqsR1 = np.array(r1f.readlines())
    r2f = open(R2)
    seqsR2 = np.array(r2f.readlines())
    length = min(len(seqsR1), len(seqsR2))
    items = random.sample(list(range(length)), number)
    ret_R1 = np.array([phred_to_values(record) for record in seqsR1[items]])
    ret_R2 = np.array([phred_to_values(record) for record in seqsR2[items]])
    return ret_R1, ret_R2


def create_fastq_parallel(sample, seqs1, quals1, seqs2, quals2, njobs=20,
                          create_template=False):
    """
    combine the quality and the sequence in fastq file

    Parameters
    ----------
    sample : str
        sample file name
    seqs1 : [str]
        list of R1 sequences
    quals1 : [int]
        list of R1 quality scores
    seqs2 : [str]
        list of R2 sequences
    quals2 : [int]
        list of R2 quality scores
    njobs : int, optional
        number of cpus used in parallel, by default 20
    create_template : bool, optional
        create a perfect set of files, by default False

    Returns
    -------
    None
        creates fastq files
    """
    def process(i, seqs1=seqs1, quals1=quals1,
                seqs2=seqs2, quals2=quals2, create_template=create_template):
        seq1 = seqs1[i]
        qual1 = quals1[i]
        if create_template:
            tem1 = seq1

        seq1 = quality_based_modification(seq1, qual1)
        record1 = SeqRecord(Seq(seq1), id="seq:R1:{}".format(i),
                            description="simulated 16S seq")
        record1.letter_annotations["phred_quality"] = qual1
        if create_template:
            perfect_template1 = SeqRecord(Seq(tem1), id="seqt:R1:{}".format(i),
                                          description="simulated 16S tempseq")
            perfect_template1.letter_annotations["phred_quality"] = [
                40 for x in qual1]

        seq2 = seqs2[i]
        qual2 = quals2[i]
        if create_template:
            tem2 = seq2
        seq2 = quality_based_modification(seq2, qual2)

        record2 = SeqRecord(Seq(seq2), id="seq:R2:{}".format(
            i), description="simulated 16S seq")
        record2.letter_annotations["phred_quality"] = qual2
        if create_template:
            perfect_template2 = SeqRecord(Seq(tem2), id="seqt:R2:{}".format(i),
                                          description="simulated 16S tempseq")
            perfect_template2.letter_annotations["phred_quality"] = [
                40 for x in qual2]

        if create_template:
            return (record1, record2, perfect_template1, perfect_template2)
        else:
            return (record1, record2)
    p = Pool(njobs)
    records = p.map(process, range(len(seqs1)))
    p.close()
    records1 = [x[0] for x in records]
    records2 = [x[1] for x in records]
    SeqIO.write(records1, "{}_R1.fastq".format(sample), format='fastq')
    SeqIO.write(records2, "{}_R2.fastq".format(sample), format='fastq')
    if create_template:
        perfect_templates1 = [x[2] for x in records]
        perfect_templates2 = [x[3] for x in records]
        SeqIO.write(
            perfect_templates1, "{}-perfecttemplate_R1.fastq".format(sample), format='fastq')
        SeqIO.write(
            perfect_templates2, "{}-perfecttemplate_R2.fastq".format(sample), format='fastq')


class SimulateFastq():
    """
    [summary]
    """

    def __init__(self, seq_template_R1, seq_template_R2):
        """
        [summary]

        Parameters
        ----------
        seq_template_R1 : [type]
            [description]
        seq_template_R2 : [type]
            [description]
        """
        self.seq_template_R1 = seq_template_R1
        self.seq_template_R2 = seq_template_R2
        self.t_R1 = np.array([str(seq.seq)
                              for seq in SeqIO.parse(seq_template_R1, 'fasta')])
        self.t_R2 = np.array([str(seq.seq)
                              for seq in SeqIO.parse(seq_template_R2, 'fasta')])

    def add_quality_template(self, qual_template_R1, qual_template_R2, length):
        """
        [summary]

        Parameters
        ----------
        qual_template_R1 : [type]
            [description]
        qual_template_R2 : [type]
            [description]
        length : [type]
            [description]
        """
        self.quality_template_length = length
        self.qR1, self.qR2 = get_quality(
            qual_template_R1, qual_template_R2, length)

    def generate_fastq(self, sample, number_of_otus, number_of_reads, create_template=False):
        """
        [summary]

        Parameters
        ----------
        sample : [type]
            [description]
        number_of_otus : [type]
            [description]
        number_of_reads : [type]
            [description]
        create_template : bool, optional
            [description], by default False
        """
        rnd_otus = np.random.choice(len(self.t_R1), number_of_otus)
        ps = np.random.choice(100000, number_of_otus)
        sumps = ps.sum()
        ps = [x/sumps for x in ps]
        seq_index = list(np.random.choice(
            rnd_otus, size=number_of_reads, replace=True, p=ps))
        quality_index = list(np.random.choice(
            self.quality_template_length, number_of_reads, replace=True))
        # print(quality_index)
        q1 = self.qR1[quality_index]
        q2 = self.qR2[quality_index]

        seqs1 = self.t_R1[seq_index]
        seqs2 = self.t_R2[seq_index]

        seqs1 = [seqs1[x][:len(q1[x])] for x in range(len(seqs1))]
        seqs2 = [seqs2[x][:len(q2[x])] for x in range(len(seqs2))]

        create_fastq_parallel(sample, seqs1, q1, seqs2,
                              q2, njobs=2, create_template=create_template)


@click.command()
@click.option("--reference", nargs=2, required=True)
@click.option("--quality", nargs=2, required=True)
@click.option("--sample-names-pattern", default="simSample")
@click.option("--number-of-samples", default=1, type=int)
@click.option("--number-otus-range", nargs=2, default=[10, 500], type=int)
@click.option("--sample-size-range", nargs=2, default=[50000, 100000], type=int)
def main(reference, quality, sample_names_pattern,
         number_of_samples, number_otus_range, sample_size_range):
    """
    [summary]

    Parameters
    ----------
    reference : [type]
        [description]
    quality : [type]
        [description]
    sample_names_pattern : [type]
        [description]
    number_of_samples : [type]
        [description]
    number_otus_range : [type]
        [description]
    sample_size_range : [type]
        [description]
    """
    simulator = SimulateFastq(reference[0], reference[1])
    print("...")
    simulator.add_quality_template(quality[0], quality[1], 500000)
    print("......")

    samples = [sample_names_pattern +
               "{:02d}".format(x) for x in range(1, number_of_samples)]
    otus = random.sample(
        range(number_otus_range[0], number_otus_range[1]), number_of_samples)
    sizes = random.sample(
        range(sample_size_range[0], sample_size_range[1]), number_of_samples)
    for sample, otu, size in zip(samples, otus, sizes):
        print(sample)
        simulator.generate_fastq(sample, otu, size, create_template=True)


if __name__ == "__main__":
    main()
