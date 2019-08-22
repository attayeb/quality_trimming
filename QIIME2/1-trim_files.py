#!/usr/bin/env python
from subprocess import call
from subprocess import Popen
import os
from shutil import copy, rmtree
from os.path import join
from sys import exit
input_folder = "/mnt/data/"
data_folder = "/mnt/data/fastq/"
copy_folder = "/mnt/data/copy/"
output_parent = "/mnt/data/result/"

trim_list = ["0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20",
             "22", "24", "26", "28", "30"]
for trim in trim_list:
    os.mkdir(copy_folder)

    filenames = os.listdir(data_folder)

    target_names = [
        x.replace("stage2", "T" + trim) for x in filenames]
    target_names = [x.replace("_R1", "_L001_R1_001") for x in target_names]
    target_names = [x.replace("_R2", "_L001_R2_001") for x in target_names]

    for i in range(len(filenames)):
        copy(data_folder + filenames[i],
             copy_folder + target_names[i])

    output_folder = output_parent + "T" + trim + "/"
    command = "python /mnt/scripts/auto-q.py -i " + copy_folder + " -t "\
        + trim + \
        " -o " + output_folder + \
        " -c /qiime.cfg -s trimming -n 20"
    print(command)
    print(os.listdir(data_folder))
    call(command, shell=True)
    rmtree(copy_folder)


print("finish")
