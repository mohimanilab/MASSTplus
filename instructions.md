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
