import os

input_mgf = open("./my-file.mgf", "r")
mgf_lines = input_mgf.readlines()

num_parts = 8
parts = []
current_part = 0

for i in range(num_parts):
    parts.append([])

for i, line in enumerate(mgf_lines):
    parts[current_part].append(line)
    if (line[0:8] == "END IONS"):
        parts[current_part].append("")
        current_part = (current_part + 1) % num_parts

for i in range(num_parts):
    outfile = open("./part-" + str(i) + ".mgf", "w")
    outfile.writelines(parts[i])
