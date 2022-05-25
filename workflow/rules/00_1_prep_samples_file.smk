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


# Convert raw .xlsx metadata into .csv samples file
rule convert_metadata:
    input:
        config["raw_metadata"]
    output:
        "config/samples_converted.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/convert_metadata/all.log"
        ),
    resources:
        mem_mb = 1000
    container:
        config["R_4.1.2"]
    script:
        "../scripts/convert_metadata.R"

# Read in samples file
if os.path.exists(rules.convert_metadata.output[0]):
    samples_df = pd.read_csv(rules.convert_metadata.output[0], comment = '#')
    videos = sorted(list(set(samples_df['movie_name'])))
    VIDEOS = []
    for i in videos:
        out = i.strip('.avi')
        VIDEOS.append(out)
    SAMPLES = samples_df['sample']
    QUADRANTS = samples_df['quadrant']

    ASSAYS = ["open_field", "novel_object"]

    # Create full list for zipping
    SAMPLES_ZIP = SAMPLES * len(ASSAYS)
    QUADRANTS_ZIP = QUADRANTS * len(ASSAYS)
    ASSAYS_ZIP = np.repeat(ASSAYS, len(SAMPLES))


# Copy videos from FTP to Codon (for ease downstream)
rule copy_videos:
    input:
        os.path.join(
            config["raw_data_dir"],
            "{video}.avi"
        ),
    output:
        os.path.join(
            config["working_dir"],
            "raw_videos/{video}.avi"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/copy_videos/{video}.log"
        ),
    resources:
        mem_mb = 100
    shell:
        """
        cp {input} {output} \
            2> {log}
        """

# Create a list of frames per second for each video to add
# manually to `config/samples.csv`
rule add_fps:
    input:
        videos = expand(rules.copy_videos.output,
            video = VIDEOS
        ),
        samples_file = rules.convert_metadata.output,
    output:
        csv = config["samples_file"],
    log:
        os.path.join(
            config["working_dir"],
            "logs/add_fps.log"
        ),
    resources:
        mem_mb = 200
    container:
        config["opencv"]
    script:
        "../scripts/add_fps.py"



# Set variables

#SAMPLES = samples_df["sample"]
#ASSAYS = ["open_field", "novel_object"]
#QUADRANTS = ["q1", "q2", "q3", "q4"]
#
## Remove faulty videos from sample list
#
## Read in videos to be excluded
#excl_df = pd.read_csv(config["excluded_videos"], comment = "#")
#
### Create list of variable lists
#full_list = [SAMPLES, ASSAYS, QUADRANTS]
### Create list of tuple combinations
#combos = list(itertools.product(*full_list))
#
## Remove unavailable combinations
#excl_df = excl_df.reset_index()
#for index, row in excl_df.iterrows():
#    combos.remove((row['sample'], row['assay'], row['quadrant']))
#
## Create new lists of variables
#SAMPLES = [i[0] for i in combos]
#ASSAYS = [i[1] for i in combos]
#QUADRANTS = [i[2] for i in combos]