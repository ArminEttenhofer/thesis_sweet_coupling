#! /usr/bin/env python3

from mule.JobGeneration import JobGeneration
from mule.JobData import JobData
from mule.InfoError import InfoError

import glob
import os
import re
import sys



if __name__ == "__main__":

    verbosity = 0

    jobdirs = []

    if len(sys.argv) > 1:
        jobdirs = sys.argv[1:]
    else:
        print("")
        print("Usage:")
        print("    "+sys.argv[0]+" [jobdir 1] [jobdir 2] ...")
        print("")
        sys.exit(1)

    for j in jobdirs:

        if verbosity > 5:
            print("*"*80)
            print("Job directory: "+j)
            print("*"*80)

        j = JobData(jobdir = j, verbosity=verbosity)
        d = j.get_flattened_data()

        for key, value in d.items():
            if isinstance(value, str):
                print(f"{key} => '{value}'")
            else:
                print(f"{key} => {value}")
