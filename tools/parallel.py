import os
import random
from shutil import copyfile

base_dir = "../parallel"
params_file = open(base_dir + "/params.txt", "r")
params = params_file.readlines()

number_of_indexes = int(params[0])
abs_path = params[1] # absolute path to the repo, i.e. the parent directory of "parallel"

index_list_file = open(base_dir + "/index-list.txt", "r")
index_file_paths = index_list_file.readlines()
random.shuffle(index_file_paths)

files_per_index = (len(index_file_paths) // number_of_indexes) + 1
current_index = 1
current_file_paths = []

os.mkdir(base_dir + "/indexes")

for i, index_file_path in enumerate(index_file_paths):
    current_file_paths.append(index_file_path)
    if (len(current_file_paths) == files_per_index) or (i == len(index_file_paths) - 1):
        new_path = base_dir + "/indexes/index-" + str(current_index)
        os.mkdir(new_path)
        outfile = open(new_path + "/index-list.txt", "w")
        outfile.writelines(current_file_paths)
        copyfile("../build/src/tools/load", new_path + "/load")
        bash_file = open(new_path + "/load-" + str(current_index) + ".sh", "w")
        bash_file.write("cd " + abs_path.strip() + "/parallel/indexes/index-" + str(current_index) + "\n")
        bash_file.write("./load index-list.txt -r" + "\n")
        current_index += 1
        current_file_paths = []

