#!/usr/bin/env python
import cdhit_reader

import os,sys

script_dir = os.path.dirname(os.path.realpath(__file__))
parent_dir = os.path.dirname(script_dir)
input_clst = os.path.join(parent_dir, "data", "aa.clstr")


if os.path.exists(input_clst):

    total = 0
    members = 0
    for i in cdhit_reader.read_cdhit(input_clst):
        total += 1
        members += len(i)
        print(i, file=sys.stderr)

    
    if total == 7 and members == 10:
        print("OK: {} clusters from {} sequences".format(total, members), file=sys.stderr)
    else:
        print("FAIL: {} clusters from {} sequences (expected 7, 10)".format(total, members), file=sys.stderr)
        sys.exit(1)