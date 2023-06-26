# MASST+

https://masst.ucsd.edu/masstplus/

MASST+ is an improvement on GNPS Mass Spectrometry Search Tool ([MASST](https://masst.ucsd.edu/)). MASST+ provides fast and error tolerant search of metabolomics mass spectrometry data while reducing the search time by two orders of magnitude. It is capable of querying against databases of billions of mass spectra, which was not feasible with MASST. Like MASST, MASST+ is publicly available as a web service on GNPS.

## Using MASST+

### With a Spectrum USI

If you know the [spectrum USI](https://github.com/mwang87/MetabolomicsSpectrumResolver) of a spectrum you want to search with MASST+, you can enter it directly at https://masst.ucsd.edu/masstplus/.

### Searching a spectrum in the GNPS library

![spectrum](spectrum.png)

(a) First, navigate to the spectrum of interest on the GNPS library. Here, a [Malyngamide C spectrum](https://proteomics3.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID=CCMSLIB00000001549#%7B%7D) is viewed. Next, click the "MASST+" link. (c) This opens the MASST+ tab which runs a mass spectral search and presents the results.

### Integration with Molecular Networking

![molecular networking](molecular-networking.png)

(a) Start by submitting a [new molecular networking](https://proteomics3.ucsd.edu/ProteoSAFe/index.jsp?task=a9c32880f76b4786a5a89682ed101d8f) job on GNPS (this will require you to be logged in to a GNPS account). (b) When the job has completed, click "View All Clusters With IDs". (c) This will open a new tab, where you can click "Advanced MASST" and then "MASST+ Search" (or "MASST+ Analog Search") in order to start a new MASST+ search. (d) This will open a new tab for MASST+, where the search results will display after a few seconds.

# Using the MASST+ binary

The `masst_plus` binary integrates all tools needed to index and query mass spectrum database. It is designed to be fast and memory-efficient. The provided binary is only for x84-64 linux systems.

The tool comes with two types of search: exact and error-tolerant (analog). Separate indices are needed for each search type. This tool is thus split into two subcommands: exact and analog, each of which has its own set of subcommands for creating, adding to, and querying indices.

You can always consult `masst_plus --help` for a detailed description of the available commands and flags. 

## The search types

**Exact search** matches query spectra against database spectra based on shared peaks only, while **analog search** (also known as error-tolerant search) matches query spectra against database spectra based on shared and shifted peaks. The index subcommand ingests files into the disk-based database indexes. 

## Creating Index

To get started, create a spectrum index for the desired search type using the `index` subcommand. This will ingests files into disk-based database indexes. The parameters available for each search type differ slightly, but both require a directory for the index and some spectra files to add to the index.

### Usage

To create an index for **exact search**, run the following:

```sh
./masst_plus exact index [OPTIONS] --dir <DIR> --spectra-files <PATHS>...
```
or
```sh
./masst_plus exact index [OPTIONS] -d <DIR> -s <PATHS>...
```

When unsure about parameters, it is sufficient to provide only the  `--dir <DIR>` and `--spectra-files <PATHS>` arguments. MASST+ will set `[OPTIONS]` to default. 

For advanced use, `[OPTIONS]` could include the following. Note that `--peak-tol` and `--pepmass-tol` will be fixed for the index. Flags under `Spectra` will be applied to preprocess the provided spectra files before they are added to the index, so they cannot be undone. 

Note that `--spectra-files` can take one spectrum file or a list of multiple spectrum files separated by space. MASST+ currently supports files in `mgf`, `mzXML`, and `mzML` files.

Also note that `<DIR>` should not exist as it will be created by this program. Setting `<DIR>` to an existing directory will result in an error.

<details>
<summary>Available options</summary>

```
OPTIONS:
    -d, --dir <DIR>
            This field is a path to directory to be used by index.

            This directory cannot already exist and will be created by this program. The spectra
            read from the files will be preprocessed and normalized before ingestion.

    -h, --help
            Print help information

        --peak-tol <PEAK_TOL>
            Absolute peak m/z tolerance

            This error tolerance is used to match database spectra peak m/zs to query spectra peaks.
            Peaks in the spectral database that have a difference in m/z bigger than this are
            considered to have distinct m/zs, otherwise they are considered to be (approximately)
            equal. The indexing strategy relies on this - a smaller tolerance means we divide
            spectra into more bins and this implies a more fine-grained comparison.

            [default: 0.02]

        --pepmass-tol <PEPMASS_TOL>
            Absolute pepmass tolerance (in Daltons).

            This error-tolerance is used during the indexing process to divide the spectra into bins
            of approximately similar pepmass. While scoring a query spectrum we only look at the
            spectra in the same bin (and a couple neighboring bins) as the query spectrum's pepmass.
            The indexing strategy relies on this - a smaller tolerance means we divide spectra into
            more bins and this implies a more fine-grained comparison.

            [default: 2.0]

Spectra:
        --min-peak-count <MIN_PEAK_COUNT>
            Minimum number of spectral peaks

            Any spectra with fewer peaks than this will be filtered out. This also serves as the
            minimum number of peaks that can be returned after the peak merging operation.

            [default: 1]

        --no-filter
            Disable peak filtering

        --no-merge
            Disable peak merging

        --normalize
            Enable peak intensity normalization

        --peak-filter-window <PEAK_FILTER_WINDOW>
            Window size for peak filtering (Da)

            This value is used in the filtering pre-processing step of the spectrum. This removes
            low intensity peaks by keeping only the --peaks-per-window most intense peaks in
            non-overlapping m/z windows of size --peak-filter-window. Filtering pre-preprocessing
            can be skipped with --no-filter.

            [default: 50.0]

        --peak-merge-thresh <PEAK_MERGE_THRESH>
            Threshold used to merge peaks (Da)

            Peaks with m/zs within this value of each other will be merged into a single peak with
            intensity aggregated additively. The minimum number of peaks that can remain after the
            merge pre-processing step is defined by --min-peak-count. Merging pre-preprocessing can
            be skipped with --no-merge.

            [default: 0.05]

        --peaks-per-window <PEAKS_PER_WINDOW>
            Number of peaks to keep per window

            This value is used in the filtering pre-processing step of the spectrum. This removes
            low intensity peaks by keeping only the --peaks-per-window most intense peaks in
            non-overlapping m/z windows of size --peak-filter-window. Filtering pre-preprocessing
            can be skipped with --no-filter.

            [default: 5]

    -s, --spectra-files <PATHS>...
            Paths to files containing the spectra of interest

            To use more than one spectrum add multiple, space-separated files.
```

</details>
    
To create an index for **analog (error-tolerant) search**, run the following:

```sh
./masst_plus analog index [OPTIONS] --dir <DIR> --spectra-files <PATHS>...
```
or
```sh
./masst_plus analog index [OPTIONS] -d <DIR> -s <PATHS>...
```

Similarly, feel free to leave out `[OPTIONS]` and MASST+ will use default parameters. 

For advanced use,  `[OPTIONS]` are as follows. Similar to exact indexes, `--peak-tol` will be fixed for the index, though `--pepmass-tol` is not an option here. This is because `--pepmass-tol` will be set to a very large number (`9999.99`) to enable shifted matches. Flags under `Spectra` will be applied to preprocess the provided spectra files before they are added to the index, so they cannot be undone.

<details>
<summary>Available options</summary>

```
OPTIONS:
    -d, --dir <DIR>
            This is a path to the directory to be used by the error-tolerant index.

            This directory cannot already exist and will be created by this program.

    -h, --help
            Print help information

        --peak-tol <PEAK_TOL>
            Absolute peak m/z tolerance

            The error tolerance used to match database spectra peak m/zs to query spectra peaks.
            Peaks in the spectral database that have a difference in m/z bigger than this are
            considered to have distinct m/zs, otherwise they are considered to be (approximately)
            equal. The indexing strategy relies on this - a smaller tolerance means we divide
            spectra into more bins and this implies a more fine-grained comparison.

            [default: 0.02]

Spectra:
        --min-peak-count <MIN_PEAK_COUNT>
            Minimum number of spectral peaks

            Any spectra with fewer peaks than this will be filtered out. This also serves as the
            minimum number of peaks that can be returned after the peak merging operation.

            [default: 1]

        --no-filter
            Disable peak filtering

        --no-merge
            Disable peak merging

        --normalize
            Enable peak intensity normalization

        --peak-filter-window <PEAK_FILTER_WINDOW>
            Window size for peak filtering (Da)

            This value is used in the filtering pre-processing step of the spectrum. This removes
            low intensity peaks by keeping only the --peaks-per-window most intense peaks in
            non-overlapping m/z windows of size --peak-filter-window. Filtering pre-preprocessing
            can be skipped with --no-filter.

            [default: 50.0]

        --peak-merge-thresh <PEAK_MERGE_THRESH>
            Threshold used to merge peaks (Da)

            Peaks with m/zs within this value of each other will be merged into a single peak with
            intensity aggregated additively. The minimum number of peaks that can remain after the
            merge pre-processing step is defined by --min-peak-count. Merging pre-preprocessing can
            be skipped with --no-merge.

            [default: 0.05]

        --peaks-per-window <PEAKS_PER_WINDOW>
            Number of peaks to keep per window

            This value is used in the filtering pre-processing step of the spectrum. This removes
            low intensity peaks by keeping only the --peaks-per-window most intense peaks in
            non-overlapping m/z windows of size --peak-filter-window. Filtering pre-preprocessing
            can be skipped with --no-filter.

            [default: 5]

    -s, --spectra-files <PATHS>...
            Paths to files containing the spectra of interest

            To use more than one spectrum add multiple, space-separated files.
```

</details>

### Examples

```sh
./masst_plus exact index -d path/to/exact_index1 -s path/to/one_million_50.mgf
```
This creates an exact index for the spectra in `one_million_50.mgf` and adds it to the index at `path/to/exact_index1`, with the default options for spectrum preprocessing and filtering.

```sh
./masst_plus analog index -d path/to/analog_index1 -s path/to/one_million_49.mgf path/to/one_million_50.mgf
```
This creates an error-tolerant index for the spectra in `one_million_49.mgf` and `one_million_50.mgf` and adds it to the index at `path/to/analog_index1`, with the default options for spectrum preprocessing and filtering.

```sh
./masst_plus exact index --peak-tol 0.01 --pepmass-tol 2.0 --no-filter --no-merge --normalize -d path/to/exact_index2 -s path/to/one_million_50.mgf
```
This creates an exact index for the spectra in `one_million_50.mgf` and adds it to the index at `path/to/exact_index2`, with peak tolerance of 0.01 Daltons, peptide mass tolerance of 2.0 Daltons, no peak filtering, no peak merging, and peak intensity normalization.

## Adding Spectrum to Index

The `add` command enables adding more spectra to an existing index. It is similar to that of `index` except there is no option to specify `--peak-tol` nor `--pepmass-tol` as they are already set when the target was initiated. The same options for spectrum preprocessing or filtering as in `index` are still available. Run `./masst_plus <exact/analog> add --help` to see all available options. 

### Usage

```sh
masst_plus <exact/analog> add [OPTIONS] --index-path <INDEX_PATH> --spectra-files <PATHS>...
```

or with short flags:

```sh
masst_plus <exact/analog> add [OPTIONS] -i <INDEX_PATH> -s <PATHS>...
```

### Examples

```sh
./masst_plus exact add -i path/to/exact_index1 -s path/to/one_million_0.mgf
```
This adds the spectra in `one_million_0.mgf` to the exact index at `path/to/exact_index1`, with the default options for spectrum preprocessing and filtering.

```sh
./masst_plus analog add -i path/to/analog_index1 -s path/to/one_million_0.mgf path/to/one_million_1.mgf
```
This adds the spectra in `one_million_0.mgf` and `one_million_1.mgf` to the analog index at `path/to/analog_index1`, with the default options for spectrum preprocessing and filtering.

```sh
./masst_plus exact add --no-filter --no-merge --normalize -i path/to/exact_index2 -s path/to/one_million_0.mgf
```
This adds the spectra in `one_million_0.mgf` to the exact index at `path/to/exact_index2`, with no peak filtering, no peak merging, and peak intensity normalization.

## Scoring

The `score` command enables querying a specified index for matches to a set of spectra. Two spectra are considered a match if their similarity score is above a specified threshold. 

### Usage

The general usage of the `query` command is as follows:

```sh
./masst_plus exact score [OPTIONS] --index-path <INDEX_PATH> --spectra-files <PATHS>... -t <THRESHOLD> -o <OUT_CSV>
```
```sh
./masst_plus ert score [OPTIONS] --index-path <INDEX_PATH> --spectra-files <PATHS>... -t <THRESHOLD> -o <OUT_CSV>
```

Or the shortened version:
```sh
./masst_plus exact score [OPTIONS] -i <INDEX_PATH> -i <PATHS>... -t <THRESHOLD> -o <OUT_CSV>
```
```sh
./masst_plus ert score [OPTIONS] -i <INDEX_PATH> -i <PATHS>... -t <THRESHOLD> -o <OUT_CSV>
```

where `threshold` is the the minimum similarity score which we report as a match, and `out_csv` is the path to the output file. The same options for spectrum preprocessing or filtering as in `index` are still available. Run `./masst_plus <exact/analog> query --help` to see all available options.

### Examples

```sh
./masst_plus exact score -i path/to/exact_index1 -s path/to/query1.mgf -t 0.7 -o path/to/output1.csv
```
This queries the exact index at `path/to/exact_index1` with the spectrum in `query1.mgf` and reports all matches above the threshold of 0.7 in `path/to/output1.csv`.

```sh
./masst_plus analog score -i path/to/analog_index1 -s path/to/query1.mgf -t 0.7 -o path/to/output2.csv
```
This queries the analog index at `path/to/analog_index1` with the spectrum in `query1.mgf` and reports all matches above the threshold of 0.7 in `path/to/output2.csv`.

```sh
./masst_plus exact score --no-filter --no-merge --normalize -i path/to/exact_index2 -s path/to/query1.mgf -t 0.8 -o path/to/output3.csv
```
This queries the exact index at `path/to/exact_index2` with the spectrum in `query1.mgf` and reports all matches above the threshold of 0.8 in `path/to/output3.csv`, with no peak filtering, no peak merging, and peak intensity normalization.

### Result Format

The output will be in `csv` format. Each row of the `csv` output represents a match above the threshold. The columns of the output represent:
- `query_file` is the spectrum file where a query spectrum comes from
- `query_scan` is a unique ID assigned to the query spectrum
- `db_file` is the spectrum file where where the matched spectrum comes from
- `db_scan` is a unique ID assigned to the matched spectrum
- `score` is the similarity score between the query and matched spectrum

An example output (with whitespaces added for readability) is as follows:

```csv
query_file,                query_scan,    db_file,                      db_scan,    score
.../query_spectrum.mgf,    68,            .../database_spectrum.mgf,    313,        0.8580977244141167
.../query_spectrum.mgf,    68,            .../database_spectrum.mgf,    518,        0.8611318986541258
.../query_spectrum.mgf,    72,            .../database_spectrum.mgf,    1492,       0.884158366782821
.../query_spectrum.mgf,    80,            .../database_spectrum.mgf,    1446,       0.7219262309070071
.../query_spectrum.mgf,    82,            .../database_spectrum.mgf,    1251,       0.8573599725776676
.../query_spectrum.mgf,    87,            .../database_spectrum.mgf,    1259,       0.846181240320242
.../query_spectrum.mgf,    87,            .../database_spectrum.mgf,    1266,       0.819756932669402
.../query_spectrum.mgf,    87,            .../database_spectrum.mgf,    1272,       0.7644365254005531
.../query_spectrum.mgf,    162,           .../database_spectrum.mgf,    313,        0.8654356609871068
.../query_spectrum.mgf,    162,           .../database_spectrum.mgf,    593,        0.8409077876547436
.../query_spectrum.mgf,    162,           .../database_spectrum.mgf,    626,        0.8345826387897034
```

which corresponds to the following table:

| query_file | query_scan | db_file | db_scan | score |
| --- | --- | --- | --- | --- |
| .../query_spectrum.mgf | 68 |         .../database_spectrum.mgf | 313 |     0.8580977244141167 | 
| .../query_spectrum.mgf | 68 |         .../database_spectrum.mgf | 518 |     0.8611318986541258 | 
| .../query_spectrum.mgf | 72 |         .../database_spectrum.mgf | 1492 |    0.884158366782821 | 
| .../query_spectrum.mgf | 80 |         .../database_spectrum.mgf | 1446 |    0.7219262309070071 | 
| .../query_spectrum.mgf | 82 |         .../database_spectrum.mgf | 1251 |    0.8573599725776676 | 
| .../query_spectrum.mgf | 87 |         .../database_spectrum.mgf | 1259 |    0.846181240320242 | 
| .../query_spectrum.mgf | 87 |         .../database_spectrum.mgf | 1266 |    0.819756932669402 | 
| .../query_spectrum.mgf | 87 |         .../database_spectrum.mgf | 1272 |    0.7644365254005531 | 
| .../query_spectrum.mgf | 162 |        .../database_spectrum.mgf | 313 |     0.8654356609871068 | 
| .../query_spectrum.mgf | 162 |        .../database_spectrum.mgf | 593 |     0.8409077876547436 | 
| .../query_spectrum.mgf | 162 |        .../database_spectrum.mgf | 626 |     0.8345826387897034 | 


# PAIRING+


Pairing+ is a part of networking+ , an improvement on the GNPS spectrometry networking tool (http://proteomics.ucsd.edu/research-areas/spectral-networks/). Pairing+ calculates pairwise dot-products between metabolomics mass spectrometry data and generates connected component graph for the molecular network. Pairing+ reduces pairwise dot product calculation time of molecular networking by over two orders of magnitude. It's capable of performing all versus all dot product calculation and generate connected component graph on databases of over 10 million spectra n a single thread, which is not feasible with spectral-networks. 


#Using PAIRING+ Source Code
The pairing_plus binary integrates the tools needed for performing networking on clustered mass spectrum database. The binary is only for x84-64 linux systems. You can always consult `pairing_plus --help` for a detailed description of the available commands and flags

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
