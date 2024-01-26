# MASST+

## Indexing

#### Usage:

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

#### Notes:
- By default, the first time this is run, a `library` directory will be created automatically in the current directory. On subsequent runs, spectra will be added to the existing library. (Command must be run from the same directory for this to work.) Alternatively, just pass in the `-l` flag to the `load` and `search` commands to specify a library directory. This directory will be created if it doesn't exist when the `load` command is run.
- Searches on this library need to be executed from the same directory, or with the same `-l` argument.

## Searching

#### Usage:

```
./search <query_file_name>
                [--analog]
                [--peaktol <peak-tolerance>]
                [--thresh <threshold>]
       --analog, -a: Run analog search (without this, exact search is run)
       --peaktol, -p <peak-tolerance>: Specify the peak tolerance (peak masses +- this amount will be considered a match; default 0.02)
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

#### Notes:
- Spectra with a match score above the threshold will be listed in the output file `matches-all.tsv` created in the directory where the search is run, unless specified otherwise with the `-o` flag.
