# PAIRING+


Pairing+ is a part of networking+ , an improvement on the GNPS spectrometry networking tool (http://proteomics.ucsd.edu/research-areas/spectral-networks/). Pairing+ calculates pairwise dot-products between metabolomics mass spectrometry data and generates connected component graph for the molecular network. Pairing+ reduces pairwise dot product calculation time of molecular networking by over two orders of magnitude. It's capable of performing all versus all dot product calculation and generate connected component graph on databases of over 10 million spectra n a single thread, which is not feasible with spectral-networks. 


#Using PAIRING+ Source Code
The pairing_plus binary integrates the tools needed for performing networking on clustered mass spectrum database. The binary is only for x84-64 linux systems. You can always consult `pairing_plus -help` for a detailed description of the available commands and flags

### Compile source code with O3 flag:
```
    g++ -O3 pairing_plus.cpp -o pairing_plus
```

### Usage
 To perform all-versus-all pairwise product calculation from a single spectra file, run the following
```sh
  ./pairing_plus [OPTIONS] -i <PATH> -o <PATH>...
```
 To perform all-versus-all pairwise product calculation from a file containing paths to a list of spectra files, run the following

 ```sh
  ./pairing_plus [OPTIONS] -l <PATH> -o <PATH>...
```
The `-o <PATH>` should take in the path of the output file containing all pairwise dot products larger than a certain threshold. The threshold can be modified through the argument `-t <THRESHOLD>` , which is set to be 0.7 by default.

When unsure about parameters, it is sufficient to provide only the  `-i <PATH>` or `-l <PATH>` and `-o <PATHS>` arguments. PAIRING+ will set `[OPTIONS]` to default.

Note that you should only provide either `-i <PATH>` and `-l <PATH>`, when using `-l <PATH>` flag, the argument should be a file in which each line corresponds to the path of a spectra file. 


<details>
<summary>Available options</summary>

```
OPTIONS:
    -i <PATH>,
        This field is a path to a single input spectra file.
    -l <PATH>
        This field is a path to a list file containing the paths to input spectra files
    -o <output_path> 
        output path for dot product results 
    
    -h, --help
            Print help information
    -t <THRESHOLD>
        The minimum threshold of preserving dot product between two spectra
        While performing networking we only preserve the edge between two spectra  if their similarity score is larger than this threshold
        [default: 0.7]
    
    --peak_resolution 
        Absolute peak m/z tolerance 
        
        The error tolerance used to match peak m/zs from different spectras in the database.
        Peaks from different spectras in the spectral database that have a difference in m/z bigger than this are
        considered to have distinct m/zs, otherwise they are considered to be (approximately)
        equal. The indexing strategy relies on this - a smaller tolerance means we divide 
        spectras into more bins and this implies a more fine-grained comparison.
        [default: 0.01]
        
    --min_peakmatches <MIN_MATCHES>
        The minimum number of matched peaks between two spectras
        
        Any pair of spectras need to have at least this number of matched peaks in the same location for its edge to be preserved in the network
        
        [default: 6]
    
    --edge_limit <EDGE_LIMIT>: 
        The largest number of refscans each query is connected to in the network
        
        For a given query inside the database, we only consider the TOP N spectra in the databse with the largest similarity score when building the connected component graph
        
        [default: 10]
    
    --connected_components_out <CONNECTED_COMPONENTS_OUT>
        The output file containing the connected components information of spectral networking
        
        [default:  cluster_info_with_cc.tsv]
```
</details>

### Examples 
```sh
./pairing_plus -i path/to/input.mgf  -o path/to/out.tsv 
```
This reads in the spectras from input.mgf and outputs the networking dot product results to out.tsv with the default options

```sh
./pairing_plus -l path/to/list.txt  -o path/to/out.tsv 
```
This reads in the spectras from spectra file path contained in list.txt and outputs the networking dot product results to out.tsv with the default options

```sh
./pairing_plus -i path/to/list.txt  -o path/to/out.tsv -t 0.8 --edge_limit 20
```
This reads in the spectras from input.mgf and outputs the networking dot product results to out.tsv. Each pair of spectra need to have a similarity score of at least 0.8 for them to be connected. If a spectra is connected to more than 20 other spectras, we only preserve the top 20 connections with the highest similarity score when forming the connected component graph.

```sh
./pairing_plus -l path/to/list.txt  -o path/to/out.tsv -t 0.8 --peak_resolution 0.02 --min_peakmatches 4
```
This reads in the spectras from spectra file path contained in list.txt and outputs the networking dot product results to out.tsv. When calculating the similarity scores between two spectras, their peaks needs to be within 0.02 Dalton range to be considered in the same location. Each pair of spectra need to have a similarity score of at least 0.8 and at least 4 different matched peaks for them to be connected. 


###Result Format
The output will be in `tsv` format. Each row of the `tsv` output represents a match above the threshold. The columns of the output represent:
- `scan_1` is a unique ID assigned to the first spectra in the pair
- `mz_1` is the precursor mass of the first spectra in the pair
- `scan_2` is a unique ID assigned to the second spectra in the pair
- `mz_2` is the precursor mass of the second spectra in the pair
- `dot_product` is the similarity score between the two spectras
- `dot_product_shared` is the shared peak product component in the similarity score
- `dot_product_shifted` is shifted peak product component in the similarity score

An example output is as follows:

|scan_1 |  mz_1 |   scan_2 | mz_2 | dot_product | dot_product_shared | dot_product_shifted
| --- | --- | --- | --- | --- | --- | --- |
| 8453819 | 2223.57 | 8453821 | 2224.57 | 0.948327 | 0.948327 | 0 |
| 8453818 | 2222.58 | 8453820 | 2223.57 | 0.970785 | 0.970785 | 0 |
| 8453817 | 2221.58 | 8453818 | 2222.58 | 0.981046 | 0.981046 | 0 |
| 8453817 | 2221.58 | 8453820 | 2223.57 | 0.967376 | 0.967376 | 0 |
| 8453815 | 2211.6  | 8453819 | 2223.57 | 0.949034 | 0.949034 | 0 |
| 8453815 | 2211.6  | 8453821 | 2224.57 | 0.945403 | 0.945403 | 0 |
