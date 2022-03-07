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

## Author

* [Andrea Telatin](https://github.com/telatin)

## License

This project is licensed under the MIT License.

## Acknowledgments

This module was based on [fasta_reader](https://github.com/EBI-Metagenomics/fasta-reader-py) by [Danilo Horta](https://github.com/horta)

