import sys
import os
import subprocess
import readline
import string
import csv, subprocess

directory = '/user/rgoldouz/NewAnalysis2020/Analysis/bin/Jobs'
notSubmitted = []
for subdir, dirs, files in os.walk(directory):
    for f in files:
        if '.sh' not in f:
            continue
        qsub = "qsub -q localgrid  -o STDOUT/" + f.split('.')[0] + ".stdout -e STDERR/" + f.split('.')[0] + ".stderr Jobs/" + f
        exit_status = subprocess.call(qsub, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            notSubmitted.append(qsub)

print 'submit those jobs that are unsubmitted'
for q in notSubmitted:
    os.system(qsub)
