#!/home/nibiohnproj9/attayeb/miniconda3/envs/qiime2-2019.1/bin/python

import subprocess
import click
import qiime2
import biom
import os, zipfile
import pandas as pd
from qiime2 import Visualization, Artifact, Metadata
from os.path import join
from pandas import read_csv, unique
from os import listdir, mkdir
import shutil

# region functions


def extract_one_file(qzafile, filename, filepath):

    with zipfile.ZipFile(qzafile) as z:
        file = [x.filename for x in z.filelist if filename in x.filename][0]
        with open(filepath, 'wb') as f:
            f.write(z.read(file))


def extract_from_artifact(art, filename, dist):
    source = os.path.join(art._archiver.data_dir, filename)
    shutil.copy(source, dist)


def create_phylogeny_tree(seq):
    from qiime2.plugins.alignment.methods import mafft, mask
    from qiime2.plugins.phylogeny.methods import fasttree, midpoint_root
    return midpoint_root(fasttree(mask(mafft(seq, n_threads=35).alignment).masked_alignment).tree).rooted_tree


def create_manifest(folder, filename="std", filt=lambda x: "_R1_" in x):
    if filename == "std":
        printf = False
    else:
        printf = True
    if printf:
        f = open(filename, "w")
    files = listdir(folder)
    files = [f for f in files if "fastq" in f]
    files.sort()
    sample_names = [x.split("_")[0] for x in files]
    # print(sample_names)
    filenames = [folder+x for x in files]
    fr = ["forward" if filt(x) else 'reverse' for x in files]
    # print(fr)
    if printf:
        f.write("sample-id,absolute-filepath,direction\n")
    else:
        print("sample-id,absolute-filepath,direction")
    for i in range(len(files)):
        if printf:
            f.write("{},{},{}\n".format(sample_names[i], filenames[i], fr[i]))
        else:
            print("{},{},{}".format(sample_names[i], filenames[i], fr[i]))
    try:
        f.close()
    except:
        pass


def filenames(art):
    return os.listdir(art._archiver.data_dir)


def folder(art):
    return art._archiver.data_dir


def create_biom_table(otu_table, taxonomy_table, filename="", message="QIIME2"):
    """
    create_biom_table(otu_table, taxonomy_table, filename="", message = "QIIME2")
    """

    ret = otu_table.view(biom.Table)
    ret.type = "OTU table"
    md = taxonomy_table.view(pd.DataFrame)
    md.index.name = '#OTU ID'
    md.columns = ['taxonomy', 'confidence']
    ret.add_metadata(md.to_dict("index"), axis="observation")
    if filename != "":
        with open(filename, "w") as f:
            f.write(ret.to_json(message))
    else:
        return ret


def run_command(cmd, verbose=True):
    print("Running external command line application. This may print "
          "messages to stdout and/or stderr.")
    print("The command being run is below. This command cannot "
          "be manually re-run as it will depend on temporary files that "
          "no longer exist.")
    print("\nCommand:", end=' ')
    print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def bbduk(artifact, trimming_threshold):
    from q2_types.per_sample_sequences import (
        SingleLanePerSamplePairedEndFastqDirFmt,
        FastqManifestFormat,
        YamlFormat)
    art_ = artifact.view(SingleLanePerSamplePairedEndFastqDirFmt)
    manifest_o = pd.read_csv(os.path.join(
        str(art_), art_.manifest.pathspec), header=0, comment='#')
    manifest = manifest_o.copy()
    manifest.filename = manifest.filename.apply(
        lambda x: os.path.join(str(art_), x))
    id_to_fps = manifest.pivot(
        index='sample-id', columns='direction', values='filename')
    result = SingleLanePerSamplePairedEndFastqDirFmt()
    for _, (__, (fwd_fp, rev_fp)) in enumerate(id_to_fps.iterrows()):
        path1 = os.path.split(fwd_fp)[1]
        path2 = os.path.split(rev_fp)[1]

        p1 = str(os.path.join(result.path, path1))
        p2 = str(os.path.join(result.path, path2))

        cmd = ['bbduk.sh', '-in1=' + fwd_fp, '-in2=' + rev_fp, '-out1=' + p1, '-out2=' + p2, '-trimq='+str(trimming_threshold),
               '-k=18', '-ktrim=f', '-qtrim=r']
        run_command(cmd)

    result.manifest.write_data(art_.manifest.view(
        FastqManifestFormat), FastqManifestFormat)
    result.metadata.write_data(art_.metadata.view(YamlFormat), YamlFormat)
    return Artifact.import_data('SampleData[PairedEndSequencesWithQuality]', result)

# endregion


@click.command(options_metavar='<options>')
@click.argument("infolder", type=click.Path(exists=True), metavar='<infolder>')
@click.argument("outfolder", type=click.Path(exists=False), metavar='<outfolder>')
@click.option("--name", metavar='<name>')
@click.option("--cores", metavar='<cores>')
def main(infolder, outfolder, name, cores):
    """This script automates QIIME2 data analysis using DADA2 and SILVA"""
    cores = int(cores)
    mkdir(outfolder)

    manifest_file = join(outfolder, "{}_manifest.txt".format(name))

    print("create manifest")
    create_manifest(folder=infolder, filename=manifest_file)

    print("Import data")
    art = Artifact.import_data('SampleData[PairedEndSequencesWithQuality]', manifest_file,
                               'PairedEndFastqManifestPhred33')
    art.save(join(outfolder, "{}_artifact.qza".format(name)))
    print("Dada2")
    from qiime2.plugins.dada2.methods import denoise_paired

    denoised = denoise_paired(art, trunc_len_f=0, trunc_len_r=0,
                              min_fold_parent_over_abundance=8.0, n_threads=cores)
    denoised.table.save(join(outfolder, "{}_table.qza".format(name)))
    denoised.representative_sequences.save(
        join(outfolder, "{}_representative_sequences.qza".format(name)))
    denoised.denoising_stats.save(
        join(outfolder, "{}_denoising_stats.qza".format(name)))
    # endregion

    print("Load classifier")
    classifier = Artifact.load(
        "/home/nibiohnproj9/attayeb/databases/qiime2/classifiers/silva-132-99-nb-classifier.qza")
    # endregion
    print("Taxonomy assignment")
    from qiime2.plugins.feature_classifier.methods import classify_sklearn

    tax_sk = classify_sklearn(
        denoised.representative_sequences, classifier, n_jobs=cores)
    tax_sk.classification.save(join(outfolder, "{}_taxonomy.qza".format(name)))

    print("Create Phylogeny tree")
    from qiime2.plugins.alignment.methods import mafft, mask

    from qiime2.plugins.phylogeny.methods import fasttree, midpoint_root

    phytree = create_phylogeny_tree(denoised.representative_sequences)
    phytree.save(join(outfolder, "{}_tree.qza".format(name)))
    extract_from_artifact(phytree, "tree.nwk", join(
        outfolder, "{}.tre".format(name)))
    print("Create Biom table")
    create_biom_table(denoised.table, tax_sk.classification,
                      join(outfolder, "{}.biom".format(name)))
    print("analysis finished")


if __name__ == "__main__":
    main()
