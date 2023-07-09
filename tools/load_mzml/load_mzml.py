import os
import argparse
from subprocess import Popen, PIPE
import ming_spectrum_library
import time

# process_size = 1000 # number of spectra to process in memory before writing to disk

parser = argparse.ArgumentParser(description='Load mzML/mzXML into index')
parser.add_argument('input_file_list', help='List of input mzML/mzXML files to load into index')
parser.add_argument('binary_path', help='binary_path')
args = parser.parse_args()

p = Popen([args.binary_path + " " + args.input_file_list + " -s"], shell=True, stdout=PIPE, stdin=PIPE)
input_file_list = open(args.input_file_list, "r").readlines()

def commit_db(p):
    line = bytes("COMMIT\n", 'UTF-8')  # Needed in Python 3.
    p.stdin.write(line)
    p.stdin.flush()
    print("Committed changes to disk\n")


i = 0
while i < len(input_file_list):
    input_file_path = input_file_list[i]
    try:
        input_spectra = ming_spectrum_library.SpectrumCollection(input_file_path.strip())
        input_spectra.load_from_file()
        input_lines = input_spectra.to_mgf_lines()
        file_path_bytes = bytes(input_file_path.strip().replace(" ", "\ ") + "\n", 'UTF-8')  # Needed in Python 3.
        p.stdin.write(file_path_bytes)
        p.stdin.flush()
        for line in input_lines:
            line = bytes(line + "\n", 'UTF-8')  # Needed in Python 3.
            p.stdin.write(line)
            p.stdin.flush()
        line = bytes("!\n", 'UTF-8')  # Needed in Python 3.
        p.stdin.write(line)
        p.stdin.flush()
        print("Loaded spectra from file: " + input_file_path)
        # if ((i % process_size) == (process_size - 1)):
        #     commit_db(p)
    except:
        print("Unable to process file: " + input_file_path)
    i = i + 1

line = bytes("END_STREAM\n", 'UTF-8')  # Needed in Python 3.
p.stdin.write(line)
p.stdin.flush()
print("Committing index...")
line = bytes("DONE\n", 'UTF-8')  # Needed in Python 3.
p.stdin.write(line)
p.stdin.flush()
print("Done")
time.sleep(86400)
