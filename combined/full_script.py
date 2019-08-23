# input
# 1. reference database.
# 2. sample fastq-files to extract quality values from
# 3. primers
# 4. number of samples
# 5. Trimming thresholds
# 6. Percentage of mismatches for qiime 1
import sys
from yaml import Loader, load
from os.path import join
pm = load(open("parameters.yaml", "r"), Loader)
sys.path.append("../fastqsimulator")
of = pm['out_put_folder']

from fastqsimulator.extract_quality_from_fastq_files import extract_quality_scores

# read the reference database
from fastqsimulator.trim_fasta_using_primer import trim_fasta
trim_fasta(reference_fasta=pm['database'], output_file_name="ref",
           outputfolder=of, reference_name=pm['reference_database_name'], )
# extract_quality_scores:
extract_quality_scores(pm['template_samples_folder'],
                       join(of, "R1q.txt"), join(of, "R2q.txt"))

#
