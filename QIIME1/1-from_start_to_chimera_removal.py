import os
from subprocess import call
from shutil import copy, rmtree
from os.path import join
input_folder = "" # The base folder
data_folder = join(input_folder, "fastq")
copy_folder = join(input_folder, "copy")
output_parent = join(input_folder, "result")
# list of trimming thresholds values
trim_list = ["0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20"]
# list of percentage of maximum difference values
p_list = ["0", "2", "4", "6", "8", "10", "12", "14", "16", "18", "20"]
for trim in trim_list:
    for p in p_list:
        os.mkdir(copy_folder)

        filenames = os.listdir(data_folder)

        target_names = [
            x.replace("sample", "T" + trim + "P" + p + "-") for x in filenames]
        target_names = [x.replace("_R1", "_L001_R1_001") for x in target_names]
        target_names = [x.replace("_R2", "_L001_R2_001") for x in target_names]

        for i in range(len(filenames)):
            copy(data_folder + filenames[i],
                 copy_folder + target_names[i])

        output_folder = output_parent + "T" + trim + "P" + p + "/"
        command = "python /mnt/scripts/auto-q.py -i " + copy_folder + " -o / -t "\
            + trim + " -p " + p + " -o " + output_folder + \
            " -c /qiime.cfg -s chimera_removal -n 20"

        call(command, shell=True)
        rmtree(copy_folder)


print("finish")
