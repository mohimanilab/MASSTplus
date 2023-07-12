# README Outline

Below we will detail how to compile and run the source code for MASST+ and Networking+ (which itself is the combination of Clustering+ and Pairing+).   In the second to last section of this README, we further detail how to use a MASST+ server (accesible with a url link)  that allows searching any given spectrum against GNPS. In the last section, we show how to access the results of running molecular networking on the entirety of GNPS.

## How to Build 
Requirements: gcc-8 or above, cmake, and x84-64 linux systems.

Both MASST+ and Networking+ can be built simultaneously with the following:
```sh
 cd <root_directory>
 mkdir build
 cd build
 cmake ..
 make
```


# MASST+

### Location of the compiled binary:
```
   build/masst_plus/tools
```

### Creating Index

In order to run a MASST+ search, first the database spectra need to be indexed. This can be done using the `load` binary.  

```
./load <library_file(list)_name>
              [--reference-list]
       --reference-list, -r: indicate library is a file list
```

#### Examples:

- `./load my_library.mgf` (for single mgf file)
- `./load my_library_list.txt -r` (see example library list below)
- `./load my_library_list.txt -r -l /my/library/path` (specify path to library, absolute or relative to current directory)

Library list file should look like:
```
file_1.mgf
file_2.mgf
file_3.mgf
file_xyz.mgf
```

By default, the first time this is run, a `library` directory will be created automatically in the current directory. On subsequent runs, spectra will be added to the existing library. (Command must be run from the same directory for this to work.) Alternatively, just pass in the `-l` flag to the `load` and `search` commands to specify a library directory. This directory will be created if it doesn't exist when the `load` command is run.

Searches on this library need to be executed from the same directory, or with the same `-l` argument.

### Search
MASST+ can conduct both an analog (error tolerant) and exact search.

```
./search <query_file_name>
                [--analog]
                [--peaktol <peak-tolerance>]
                [--thresh <threshold>]
       --analog, -a: Run analog search (without this, exact search is run)
       --peaktol, -p <peak-tolerance>: Specify the peak tolerance (peak masses +- this amount will be considered a match; default 0.02)
       --precursortol, -q <precursor-tolerance>: Specify the precursor tolerance (precursor masses +- this amount will be considered a match; default 0.025)
       --thresh, -t <threshold>: specify matching score threshold for search, default 0.7
```

#### Examples:

- `./search my_query.mgf` (exact search with 0.02 peak tolerance, 0.025 precursor tolerance, 0.7 score threshold)
- `./search my_query.mgf -a` (analog search with 0.02 peak tolerance, 0.025 precursor tolerance, 0.7 score threshold)
- `./search my_query.mgf -a -p 0.01` (analog search with 0.01 peak tolerance, 0.025 precursor tolerance, 0.7 score threshold)
- `./search my_query.mgf -a -q 0.03` (analog search with 0.02 peak tolerance, 0.03 precursor tolerance, 0.7 score threshold)
- `./search my_query.mgf -a -p 0.01 -t 0.8` (analog search with 0.01 peak tolerance, 0.025 precursor tolerance, 0.8 score threshold)
- `./search my_query.mgf -a -l my/library/path` (analog search with library path specified)
- `./search my_query.mgf -a -o my/output/file.tsv` (analog search with matches output file specified)

Spectra with a match score above the threshold will be listed in the output file `matches-all.tsv` created in the directory where the search is run, unless specified otherwise with the `-o` flag.


<!---               ### Setting up on linux server 

This will set up the clustered and unclustered index from pre-built indexes.

1. After cloning this repo, go to the [setup](/setup) folder
1. `bash setup.sh`

Note that the last step will be downloading upwards of a terabyte of data and then extracting it. This will take a long time.

### Test cases

After setting up on a server as described above, you can test your setup with the test cases in the [test](test) directory. It is divided into analog and exact search for both the clustered and unclustered index. For each combination, the query (.mgf file) along with the expected output produced using that query (.tsv file) are given.

You can run the given bash scripts in this directory to run the tests (the queries and script will need to be moved to the location of the index, or you will need to specify the library location with the `-l` argument). --> 




# NETWORKING+

NETWORKING+ is a combination of CLUSTERING+ and PAIRING+. CLUSTERING+ reduces the redundancy in a spectral dataset by clustering similar spectra together and PAIRING+ computes a shifted score between all pairs of clusters. NETWORKING+ is almost three orders of magnitude faster than the previous state of the art GNPS Molecular Networking software (https://ccms-ucsd.github.io/GNPSDocumentation/networking/).


### Location of the compiled binary:
 CLUSTERING+ can be found at
```
build/networking_plus/clustering_plus
```


PAIRING+ can be found at
```
build/networking_plus/pairing_plus
```

### CLUSTERING+ Usage
 To perform clustering on a single spectra file, run the following
```sh
  ./clustering_plus [OPTIONS] -i <PATH> -o <PATH>...
```
 To perform clustering on a file containing paths to a list of spectra files, run the following

 ```sh
  ./clustering_plus [OPTIONS] -l <PATH> -o <PATH>...
```
The `-o <PATH>` should take in the path of the output file containing the cluster information. The threshold can be modified through the argument `-t <THRESHOLD>` , which is set to be 0.7 by default.

When unsure about parameters, it is sufficient to provide only the  `-i <PATH>` or `-l <PATH>` and `-o <PATHS>` arguments. CLUSTERING+ will set `[OPTIONS]` to default.

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
        output path for cluster results 
    
    -h, --help
            Print help information
    -t <THRESHOLD>
        The minimum similarity threshold with a cluster center for joining the cluster 
        [default: 0.7]
    
    -c <CLUSTER_CENTER>
        Output path containing the representative spectra of each cluster
        
        [default: centers.mgf]
    
    -s <CLUSTER_INFORMATION>
        Output CSV file containing the content information of each cluster
        
        [default: cluster_info.csv]
    
    -f <PEAK_FILTER_WINDOW> <TOP_PEAKS>
        Window size for peak filtering (Da) and number of preserved peaks in each window
        This value is used in the filtering pre-processing step of the spectrum. This removes
        low intensity peaks by keeping only the <TOP_PEAKS> most intense peaks in non-overlapping m/z windows of size <PEAK_FILTER_WINDOW>. 
        
        [default: 50.0, 5]
    
    --pepmass_resolution
        Absolute pepmass tolerance (in Daltons)
        
        This error-tolerance is used during the indexing process to divide the spectra into bins
        of approximately similar pepmass. While scoring a query spectrum we only look at the
        spectra in the same bin (and a couple neighboring bins) as the query spectrum's pepmass.
        The indexing strategy relies on this - a smaller tolerance means we divide spectra into
        more bins and this implies a more fine-grained comparison.
        
        [default: 1.0]
        
    --peak_resolution  <PEAK_RESOLUTION>
        Absolute peak m/z tolerance 
        
        The error tolerance used to match query spectra peak m/zs to cluster center spectra m/zs.
        Peaks from cluster center spectras that have a difference in m/z bigger than this are
        considered to have distinct m/zs, otherwise they are considered to be (approximately)
        equal. The indexing strategy relies on this - a smaller tolerance means we divide 
        spectras into more bins and this implies a more fine-grained comparison.
        
        [default: 0.01]
    
    --topfour_threshold <TOPFOUR_THRESHOLD>
        The minimum number of spectra within within the same --pepmass_resolution range to apply top-four peaks filtering
        
        Top-four peaks filtering only considers the possibility of a query belonging to a cluster 
        
        if there are matched peaks within --topfour_resolution between the top four peaks of the query and top four peaks of the cluster center.
        
        can be skipped with --no-topfour.
        
        [default: 5000]
        
    --topfour_resolution <TOPFOUR_RESOLUTION>
        Absolute m/z tolerance of topfour peaks.
        
        The error tolerance used to match query spectra peak m/zs to cluster center spectra m/zs when conducting topfour filtering.
        
        can be skipped with --no-topfour.
        
        [default: 0.1]
         
    --no-topfour
        Disable topfour filtering
    
    --cluster_minsize <CLUSTER_MINSIZE>: 
        the minimum cluster size for writing the cluster center into the output file
        
        Usually for each peptide we observe multiple spectras in the database. Small clusters or singletons are often resulted from noise when generating mass spectrometry file.
        
        Clusters smaller than this amount will be considered as noise data and filtered out before writing the output file

        [default: 2]
    
    --connected_components_out <CONNECTED_COMPONENTS_OUT>
        The output file containing the connected components information of spectral networking
        
        [default:  cluster_info_with_cc.tsv]
        
    --bruteforce 
        Use bruteforce clustering technique rather than indexing 
        
        The clustering will be much slower if applied
```
</details>

### CLUSTERING+ Examples 
```sh
./clustering_plus -i path/to/input.mgf  -o path/to/out.tsv 
```
This reads in the spectras from input.mgf and outputs the clustering results to out.tsv with the default options

```sh
./pairing_plus -l path/to/list.txt  -o path/to/out.tsv 
```
This reads in the spectras from spectra file path contained in list.txt and outputs the clustering results to out.tsv with the default options

```sh
./clustering_plus -i path/to/list.txt  -o path/to/out.tsv -t 0.9 --no_topfour
```
This reads in the spectras from input.mgf and outputs the clustering results to out.tsv. The exact similarity score has to be at least 0.9 between a query spectra and cluster center for joining the cluster. The clustering performs without top four filtering.

```sh
./clustering_plus -l path/to/list.txt  -o path/to/out.tsv -c path/to/centers.mgf -s cluster_info.tsv -t 0.9 --peak_resolution 0.02 --pepmass_resolution 2.0 
```
This reads in the spectras from spectra file path contained in list.txt and outputs the clustering  results to out.tsv. The query spectra cluster center need to have a pepmass difference within 2.0 Dalton range and a similarity score of at least 0.9 for the spectra to be added to the cluster. When calculating the similarity scores between query spectra and clsuter center, their peaks needs to be within 0.02 Dalton range to be considered in the same location. 


### CLUSTERING+ Output Format

CLUSTERING+ produce three output files. 

The first output file shows a mapping between each input spectra and the cluster it associated with. 

An example output is as follows:

|cluster_idx |  scan | mz | RTINSECONDS | index in file | source filename | 
| --- | --- | --- | --- | --- | --- |
| 1 | 25542893 | 225.04 | 343.079 | 230 | GA204_M3_FT_IT.mgf |
| 1 | 25542935 | 225.04 | 332.666 | 226 | GL0422_M3_FT_IT.mgf |
| 8 | 25547132 | 225.06 | 13.9954 | 106 | 140709_PMA_NM_2_L2_ddMS2_pos.mgf |
| 8 | 25547655 | 225.06 | 4.38721 | 31 | 140709_PMA_NM_2_L2_ddMS2_pos.mgf |
| 8 | 25547668 | 225.05 | 12.482 | 94 | 140709_PMA_NM_2_L2_ddMS2_pos.mgf |

Here
- `cluster_idx` is a unique ID assigned to the cluster containing this spectra
- `scan` is a unique ID assigned to the spectra
- `mz` is the precursor mass of the spectra
- `RTINSECONDS` is the 
- `index in file` is the index of spectra in the original source file
- `source filename` is the path to the original source file 

The second output file contains statistics about each cluster. 

An example output cluster information is as follows:

| cluster_idx | average mz | average RT | num spectra | 
| --- | --- | --- | --- | 
| 1 | 225.04 | 337.873 | 2 | 
| 8 | 225.09 | 367.5 | 357 |
| 9 | 225.175 | 332.244 | 1579 | 
| 18 | 225.422 | 350.623 | 2 |
| 20 | 225.215 | 666.844 | 15 |

The generated clusters will be in `tsv` format. Each row of the `tsv` output represents the cluster information of a spectra:
- `cluster_idx` is a unique ID assigned to the cluster containing this spectra
- `average mz` is the average precursor mass of all spectra inside the cluster
- `average RT` is the average RTINSECONDS of all spectra inside the cluster
- `num spectra` is the number of spectra in the cluster

The third output file will be an mgf file which contains the spectrum for each cluster_idx.

```
BEGIN IONS   
CLUSTERINDEX=45   
CLUSTERSIZE=1107   
FILENAME=KP_243_Positive.mgf     
LOCAL_SCAN=3574    
PEPMASS=125.096   
RTINSECONDS=742.131   
105.067 30   
105.075 60  
105.322 40   
107.083 106  
108.086 44  
109.069 50  
109.198 18  
110.071 23  
125.096 183  
145.258 52  
1437.3 57  
1640.82 43  
END IONS  
```


### PAIRING+ Usage
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

### PAIRING+ Examples 
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


### PAIRING+ Result Format
The output will be in `tsv` format. Each row of the `tsv` output represents a match above the threshold. The columns of the output represent:
- `scan_1` is a unique ID assigned to the first spectra in the pair
- `mz_1` is the precursor mass of the first spectra in the pair
- `scan_2` is a unique ID assigned to the second spectra in the pair
- `mz_2` is the precursor mass of the second spectra in the pair
- `dot_product` is the similarity score between the two spectras
- `dot_product_shared` is the shared peak product component in the similarity score
- `dot_product_shifted` is shifted peak product component in the similarity score

An example output is as follows:

|scan_1 | mz_1 | scan_2 | mz_2 | dot_product | dot_product_shared | dot_product_shifted |
| --- | --- | --- | --- | --- | --- | --- |
| 8453371 | 1667.18 | 8453423 | 1669.19 | 0.800478 | 0.433597 | 0.366881 |  
| 8453819 | 2223.57 | 8453821 | 2224.57 | 0.948327 | 0.948327 | 0 |  
| 8453818 | 2222.58 | 8453820 | 2223.57 | 0.970785 | 0.970785 | 0 |  
| 8453467 | 1671.02 | 8453526 | 1674.04 | 0.889376 | 0.829941 | 0.0594347 |  
| 8453817 | 2221.58 | 8453820 | 2223.57 | 0.967376 | 0.967376 | 0 |  
| 8453815 | 2211.6  | 8453819 | 2223.57 | 0.949034 | 0.949034 | 0 |  
| 8453453 | 1669.7  | 8453474 | 1670.7  | 0.842288 | 0.823904 | 0.0183844 |  


# MASST+ Server

https://masst.ucsd.edu/masstplus/

MASST+ is an improvement on GNPS Mass Spectrometry Search Tool ([MASST](https://masst.ucsd.edu/)). MASST+ provides fast and error tolerant search of metabolomics mass spectrometry data while reducing the search time by two orders of magnitude. It is capable of querying against databases of billions of mass spectra, which was not feasible with MASST. Like MASST, MASST+ is publicly available as a web service on GNPS.

## Using MASST+ Server

### With a Spectrum USI

If you know the [spectrum USI](https://github.com/mwang87/MetabolomicsSpectrumResolver) of a spectrum you want to search with MASST+, you can enter it directly at https://masst.ucsd.edu/masstplus/.

### Searching a spectrum in the GNPS library

![spectrum](spectrum.png)

(a) First, navigate to the spectrum of interest on the GNPS library. Here, a [Malyngamide C spectrum](https://proteomics3.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID=CCMSLIB00000001549#%7B%7D) is viewed. Next, click the "MASST+" link. (c) This opens the MASST+ tab which runs a mass spectral search and presents the results.

### Integration with Molecular Networking

![molecular networking](molecular-networking.png)

(a) Start by submitting a [new molecular networking](https://proteomics3.ucsd.edu/ProteoSAFe/index.jsp?task=a9c32880f76b4786a5a89682ed101d8f) job on GNPS (this will require you to be logged in to a GNPS account). (b) When the job has completed, click "View All Clusters With IDs". (c) This will open a new tab, where you can click "Advanced MASST" and then "MASST+ Search" (or "MASST+ Analog Search") in order to start a new MASST+ search. (d) This will open a new tab for MASST+, where the search results will display after a few seconds.



# GNPS Molecular Network

We performed molecular networking (both clustering and spectral networking) using NETWORKING+ on the entirety of GNPS. We stored the results of CLUSTERING+ and PAIRING+ in tsv format. 

## Clustering+ Results for GNPS
We split the GNPS library into 9 divisions according to different precursor mass ranges and performed clustering+ algorithm on each of them. We provide the cluster information of each spectra and the centers for all clusters. 

### Cluster information for each spectra

The output is in `tsv` format. Each row of the `tsv` output represents a spectra from GNPS library. The columns of the output represent:
- `cluster_idx` is a unique ID assigned to each cluster in the division
- `scan` is a unique ID assigned to the each spectra in the division
- `mz` is the precursor mass of the spectra
- `RTINSECONDS` is the retention time of the spectra
- `MSV_source` is the MSV library it belongs to
- `Filename` is the GNPS source file of this spectra inside MSV library
- `Local_scan` is the spectra's scan number inside its GNPS source file


The clustering+ output files for all 9 divisions can be downloaded via the following links:


[CLUSTERING+ output for division 0](https://drive.google.com/file/d/1JAIOsRXYlNVMPkUxhwZO4vogCujFPEYZ/view?usp=sharing)

[CLUSTERING+ output for division 1](https://drive.google.com/file/d/1t9ofD5fX0MPvAKxK8sUJwSb_LY9Cyb4_/view?usp=sharing)

[CLUSTERING+ output for division 2](https://drive.google.com/file/d/1B-hRy_PuJ1kUc1Hcb9VIXTlNn9uZmxiu/view?usp=sharing)

[CLUSTERING+ output for division 3](https://drive.google.com/file/d/1Tjj-8CZoMWjk-DLAYcAovOleeqVpBSZk/view?usp=sharing)

[CLUSTERING+ output for division 4](https://drive.google.com/file/d/14VufKemZJKpwcxW5FjWJ9TZIN4S9FxhN/view?usp=sharing)

[CLUSTERING+ output for division 5](https://drive.google.com/file/d/1CYI0r0TNhYt1JBY-4Dn8tKTwSWVUQKst/view?usp=sharing)

[CLUSTERING+ output for division 6](https://drive.google.com/file/d/1jgh-V6zWWzXq36hOVB88TfiHY_NjyUQb/view?usp=sharing)

[CLUSTERING+ output for division 7](https://drive.google.com/file/d/15PCqvvegkoPv7BN0XpMr42ggPUjLisf4/view?usp=sharing)

[CLUSTERING+ output for division 8](https://drive.google.com/file/d/130ZwPV7aBrHp2Y09Vuj5Xv8jSbrvKC6A/view?usp=sharing)




### Cluster centers

We write the representative spectrum of each cluster_idx into a mgf file for each division.

The representative spectra for each cluster contains: 
- `CLUSTERINDEX` is the cluster index in the division
- `CLUSTERSIZE` is the number of spectra in the cluster
- `MSV_LIB` is the source MSV library of the representative spectra
- `FILENAME` is the source GNPS file of the representative spectra
- `LOCAL_SCAN` is the scan number of the representative spectra inside the GNPS source file
- `PEPMASS` is the percursor mass of the representative spectra
- `RTINSECONDS` is the retention time of the representative spectra

```
BEGIN IONS
CLUSTERINDEX=9
CLUSTERSIZE=282
MSV_LIB=MSV000083789
FILENAME=pos_Cd10MYY_33.mgf
LOCAL_SCAN=1069
PEPMASS=53.0051
RTINSECONDS=336.875
31.991 36
38.0024 36
38.0076 36
49.9917 111
51.9917 36
52.8466 75
53.0038 1338
53.0203 40
67.9882 72
END IONS
```


The spectrum file for each division can be downloaded via the following links:

[cluster centers for division 0](https://drive.google.com/file/d/15rYEmjh_w0NWViPAyhdAyuJTLcKoKjYr/view?usp=sharing)

[cluster centers for division 1](https://drive.google.com/file/d/1QvgNqgwssMD1xwK9E6eUwyqPrDvs_4mz/view?usp=sharing)

[cluster centers for division 2](https://drive.google.com/file/d/11yIs8dTsehk6F_Lm88Y7EdX5P1_z7UUX/view?usp=sharing)

[cluster centers for division 3](https://drive.google.com/file/d/1pmFx5eEvIddgK5cevZeYuEfC-igX0FrK/view?usp=sharing)

[cluster centers for division 4](https://drive.google.com/file/d/1Bkuz0kKzQCxaf_Qh8qdpt1WC0pJABzJ6/view?usp=sharing)

[cluster centers for division 5](https://drive.google.com/file/d/14k9btagtz-a2RAyY5hzmKmHrbjZ8_o2Z/view?usp=sharing)

[cluster centers for division 6](https://drive.google.com/file/d/1AouC-HZ-qPc6bMv94lSkIzDaVM3K48z_/view?usp=sharing)

[cluster centers for division 7](https://drive.google.com/file/d/1Hjkdvz__fSdI5Uvec_XRLN0GJxhKczsp/view?usp=sharing)

[cluster centers for division 8](https://drive.google.com/file/d/18KcAnNMdQ6R7_ypVy2my2fYr-CovLG40/view?usp=sharing)

## PAIRING+ Results for GNPS
We apply PAIRING+ to the clusters resulting from CLUSTERING+ to compute the molecular network. The network is stored in two files

### Network Nodes
 The first output file stores general information for the [nodes of the GNPS molecular network](https://drive.google.com/file/d/1_jaBgPQJRxO2Gpr8hxD1Na2QZEX_Kyju/view?usp=sharing) in `tsv` format. The network contains over 8M nodes (total number of non-singleton clusters resulting from CLUSTERING+) in total.  Each row of the `tsv` output represents a node in the network. The columns of the output represent:
- `scan_number_among_centers` is a unique ID assigned to each cluster in the network
- `component_index` is a unique ID assigned to the each connected component in the network
- `source_division` is the division this cluster came from (ranges from division0 to division8)
- `cluster_index_in_division` is the index of this cluster in its source division
- `cluster_size` is the size of the cluster
- `center_MSV` is the source MSV library of the representative spectra
- `center_source_file` is the source GNPS file of the representative spectra
- `center_scan_in_source_file` is the scan number of the representative spectra inside the GNPS source file
- `center_pepmass` is the percursor mass of the representative spectra
- `center_RT` is the retention time of the representative spectra

### Network Edges

 The second output file stores general information for the [edges of the GNPS molecular network](https://drive.google.com/file/d/1hSDVBihMiM-GxZ6wHGdBh2wJdggvF1E-/view?usp=sharing) in `tsv` format. Each row of the `tsv` output represents an edge in the network. The columns of the output represent:
 
 - `connected_component_index ` is a unique ID assigned to each connected component in the network
- `first_center_scan_number` is a unique ID assigned to the each connected node in the network
- `second_center_scan_number` is a unique ID assigned to the each connected node in the network
- `product` is the similarity dot-product between the two nodes 
- `product_shared` is the contribution of shared peak matches  in the similarity score
- `product_shifted` is the contribution of shifted peak matches in the similarity score








