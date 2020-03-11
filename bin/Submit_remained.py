import sys
import os
import subprocess
import readline
import string

directory = '/user/rgoldouz/NewAnalysis2020/Analysis/bin'
for subdir, dirs, files in os.walk(directory):
    for f in files:
        if '.sh' not in f:
            continue
        os.system("qsub -q localgrid  -o STDOUT/" + f.split('.')[0] + ".stdout -e STDERR/" + f.split('.')[0] + ".stderr " + f)
