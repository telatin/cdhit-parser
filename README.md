# cdhit-parser

[![Python package](https://github.com/telatin/cdhit-parser/actions/workflows/python-package.yml/badge.svg)](https://github.com/telatin/cdhit-parser/actions/workflows/python-package.yml)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/cdhit-reader)](https://anaconda.org/bioconda/cdhit-reader)
[![pipy](https://img.shields.io/pypi/dm/cdhit-reader?color=blue&label=pipy%20downloads)](https://pypi.org/project/cdhit-reader/)
![PyPI](https://img.shields.io/pypi/v/cdhit-reader)

CD-HIT file reader.

## Read CD-HIT .clstr file

Basic usage

```python
from cdhit_reader import  read_cdhit
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

## Read FASTA file

```python
if os.path.exists(fileName):
    for seq in cdhit_reader.read_fasta(fileName, line_len=60):
        print(seq) # will be wrapped at 60 chars per line, use 0 to disable wrapping
        
        # to access individual attributes:
        # print(">" + seq.name + " " + seq.comment + "\n" + seq.sequence)
```

## Install

```bash
pip install cdhit-reader
```

or via [Miniconda](https://telatin.github.io/microbiome-bioinformatics/Install-Miniconda/), which will also install cd-hit

```bash
conda install -c bioconda -c conda-forge cdhit-reader
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
cdhit-compare data/input1.faa data/input2.faa  --id 0.99
```

will produce:

```text
input1  BJJOHBJ_00007
input2  BJJOHBJ_00007
input2  BJJOHBJ_00002
both    BJJOHBJ_00003:BJJOHBJ_00003
both    BJJOHBJ_00005:BJJOHBJ_00005
both    BJJOHBJ_00004:BJJOHBJ_00004
multi   input1#IBJJOHBJ_00006,input1#BBJJOHBJ_000B6,input1#CBJJOHBJ_000C6,input2#IBJJOHBJ_00006,input2#BBJJOHBJ_000B6,input2#CBJJOHBJ_000C6
dupl_input1     BJJOHBJ_00001:BJJOHBJ_000F
```

where records starting with _file1_ or _file2_ are only present in one of the files,
records starting with _both_ are present in both files (one per file),
records starting with _dupl_ are duplicates (two in one of the files),
and records starting with _multi_ are present multiple times in at least one of the datasets. 

## Author

* [Andrea Telatin](https://github.com/telatin)

## License

This project is licensed under the MIT License.

## Acknowledgments

This module was based on [fasta_reader](https://github.com/EBI-Metagenomics/fasta-reader-py)
by [Danilo Horta](https://github.com/horta)
