import os

num_queries = 20

def get_results(mode):
    results = []
    results.append("Load time\tSearch time")
    for i in range(num_queries):
        log_file = open(base_dir + "/" + mode + "/test-" + str(i) + "/log-search.txt", "r")
        log_lines = log_file.readlines()
        for j, line in enumerate(log_lines):
            load_time = 99999
            search_time = 99999
            if j == 2:
                k = 25
                while line[k] != "m":
                    k++
                load_time = int(line[25:(k+1)])
            if j == 5:
                k = 30
                while line[k] != " ":
                    k++
                search_time = int(line[30:(k+1)])
        results.append(load_time + "\t" + search_time)
    results_file = open(base_dir + "/results-" + mode + ".tsv", "w")
    results_file.writelines(results)

get_results("exact")
get_results("analog")
