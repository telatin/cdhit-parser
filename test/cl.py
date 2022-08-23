#!/usr/bin/env python
import cdhit_reader

import os,sys

if os.path.exists(sys.argv[1]):
    for i in cdhit_reader.read_cdhit(sys.argv[1]):
        print(i)