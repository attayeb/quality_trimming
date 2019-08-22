
import os
from subprocess import call
from shutil import copy2
from os.path import join
original_folder = "result/"
dest_folder = "chimera_removed/"
folders = os.listdir(original_folder)
print(folders)
for folder in folders:
    if (folder.endswith("P0")) or (folder.endswith("P4")) or (folder.endswith("P8")) or (folder.endswith("P12")):
        current_folder = join(original_folder, folder, 'chi')
        files = os.listdir(current_folder)
        for file in files:
            copy2(join(current_folder, file), join(dest_folder, file))

command = "python auto-q.py -i "
