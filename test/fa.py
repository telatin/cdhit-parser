#!/usr/bin/env python
import cdhit_reader

import os,sys

if os.path.exists(sys.argv[1]):
    for i in cdhit_reader.read_fasta(sys.argv[1], line_len=60):
        print(i)