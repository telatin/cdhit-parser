# cdhit-parser

CD-HIT file reader.

## Examples

Basic usage

```python
input = "cluster.fa.clstr"
for cluster in read_cdhit(input):
    print(f"{cluster.name} refSequence={cluster.refname} size={len(cluster)}")

    for member in cluster.sequences:
        print(f" {member.name} ({member.length}) identity={member.identity}% {'(Reference sequence)' if member.is_ref else ''}")
```

Load all clusters in to a list:
```python

# Load all clusters to a list
clusters = read_cdhit(input).read_items()
```

## Install

```bash
pip install cdhit-reader
```

## Demo applications

### Cluster stats

The module ships a demo program called `cdhit-reader.py`.

![`cdhit-parser -h`](docs/chdit.svg)

### Compare two fasta files

:warning: This requires cd-hit installed and available in the system path.

`cdhit-compare` allows to compare two fasta files and print the sequences that are in common, those which are only
present in one of the files or those which are redundant.

![`cdhit-compare --help`](docs/compare.svg)

Example:

```bash
cdhit-compare data/input1.faa data/input2.faa
```

will produce:

```text
input1_ IBJJOHBJ_00007
input2_ IBJJOHBJ_00007
input2_ IBJJOHBJ_00002
both    IBJJOHBJ_00003:_IBJJOHBJ_00003
both    IBJJOHBJ_00005:_IBJJOHBJ_00005
both    IBJJOHBJ_00004:_IBJJOHBJ_00004
dupl    IBJJOHBJ_00001:IBJJOHBJ_000F1
```

where records starting with _file1_ or _file2_ are only present in one of the files,
records starting with _both_ are present in both files (one per file),
records starting with _dupl_ are duplicates (two in one of the files),
and records starting with _multi_ are present multiple times in at least one of the datasets

## Author

* [Andrea Telatin](https://github.com/telatin)

## License

This project is licensed under the MIT License.

## Acknowledgments

This module was based on [fasta_reader](https://github.com/EBI-Metagenomics/fasta-reader-py)
by [Danilo Horta](https://github.com/horta)

