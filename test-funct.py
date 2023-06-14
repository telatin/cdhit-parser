#!/usr/bin/env python3
import os
from cdhit_reader import read_cdhit

input = "./cdhit_reader/test/aa.clstr"


if not os.path.isfile(input):
    print("File not found:", input)

# Load all clusters to a list
clusters = read_cdhit(input).read_items()


# Parser cluster file
clust_count = 0
for cluster in read_cdhit(input):
    clust_count += 1
    print(
        f"\n[{clust_count}] {cluster.name} refSequence={cluster.refname} size={len(cluster)}"
    )
    m_count = 0
    for member in cluster.sequences:
        m_count += 1
        print(
            f"{'🌖' if member.is_ref else '🌑'}\t{m_count}: {member.name} ({member.length}) %ID={member.identity}%"
        )


# Print cdhit_reader version
import cdhit_reader

print("\nCD-HIT_reader version:", cdhit_reader.__version__)
