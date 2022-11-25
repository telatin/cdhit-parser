#!/usr/bin/env python
import cdhit_reader

import os,sys
script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
input_fasta = os.path.join(parent_dir, "data", "test.fa.gz")

if os.path.exists(input_fasta):
    tot = 0
    comments = 0
    bp = 0
    for i in cdhit_reader.read_fasta(input_fasta, line_len=60):
        print("[seqname=", i.name, "|comment=", i.comment, "]", sep="",file=sys.stderr)
        tot += 1
        bp  += len(i)
        if i.comment:
            comments += 1
    if tot == 3 and comments == 2 and bp == 29:
        print("OK: {} sequences from {} comments, total {} bp".format(tot, comments, bp), file=sys.stderr)
    else:
        print("FAIL: {} sequences from {} comments,  total {} bp (expected 3, 2, 29)".format(tot, comments, bp), file=sys.stderr)
        sys.exit(1)