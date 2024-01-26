import os
import json
import shutil

base_dir = "./"
index_list_file = open(base_dir + "/index-list.txt", "r")
index_file_paths = index_list_file.readlines()

os.mkdir(base_dir + "/library")
os.mkdir(base_dir + "/library/partitions")
total_spectra_ref = {}
total_num_partitions = 0
total_num_spectra = 0
total_num_peaks = 0

for i, index_file_path in enumerate(index_file_paths):
    lib_file_path = index_file_path.strip() + "/library"
    meta_file = open(lib_file_path + "/meta.txt", "r")
    meta_data = meta_file.readlines()[0].split(" ")
    total_num_spectra = total_num_spectra + int(meta_data[0])
    total_num_peaks = total_num_peaks + int(meta_data[1])
    num_partitions = int(meta_data[2])
    spectra_ref_file = open(lib_file_path + "/spectra_ref.txt", "r")
    json_blob = spectra_ref_file.read()
    spectra_ref = json.loads(json_blob)
    if (i == 0):
        total_spectra_ref = spectra_ref
        shutil.copyfile(lib_file_path + "/config.txt", base_dir + "/library/config.txt")
        shutil.copyfile(lib_file_path + "/config.txt", base_dir + "/library/config.txt")
    else:
        total_spectra_ref["value0"] = total_spectra_ref["value0"] + spectra_ref["value0"]
    partition_dirs = os.listdir(lib_file_path + "/partitions")
    if (len(partition_dirs) != num_partitions):
        print(len(partition_dirs))
        print(num_partitions)
        raise Exception("Number of partitions does not match metadata")
    for p, partition_dir in enumerate(partition_dirs):
        new_dir = base_dir + "/library/partitions/partition-" + str(total_num_partitions)
        shutil.move(lib_file_path + "/partitions/" + partition_dir, new_dir)
        total_num_partitions = total_num_partitions + 1

with open(base_dir + "/library/spectra_ref.txt", 'w') as outfile:
    json.dump(total_spectra_ref, outfile)

new_meta_data = str(total_num_spectra) + " " + str(total_num_peaks) + " " + str(total_num_partitions)
new_meta_data_file = open(base_dir + "/library/meta.txt", "w")
new_meta_data_file.write(new_meta_data)
