######################
# Libraries
######################

import pandas as pd
import numpy as np
import os
import itertools

######################
# Config file and sample sheets
######################

configfile: "config/config.yaml"

# Read in samples file
samples_df = pd.read_csv(config["samples_file"], comment = '#')

# Set variables

SAMPLES = samples_df["sample"]
ASSAYS = ["open_field", "novel_object"]
QUADRANTS = ["q1", "q2", "q3", "q4"]

# Remove faulty videos from sample list

# Read in videos to be excluded
excl_df = pd.read_csv(config["excluded_videos"], comment = "#")

## Create list of variable lists
full_list = [SAMPLES, ASSAYS, QUADRANTS]
## Create list of tuple combinations
combos = list(itertools.product(*full_list))

# Remove unavailable combinations
excl_df = excl_df.reset_index()
for index, row in excl_df.iterrows():
    combos.remove((row['sample'], row['assay'], row['quadrant']))

# Create new lists of variables
SAMPLES = [i[0] for i in combos]
ASSAYS = [i[1] for i in combos]
QUADRANTS = [i[2] for i in combos]