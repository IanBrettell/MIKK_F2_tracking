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

##########################
# Set variables
##########################

# Read in samples file
samples_df = pd.read_csv(config["samples_long"])

VIDEOS = []
for i in sorted(list(set(samples_df['movie_name']))):
    VIDEOS.append(i.strip('.avi'))

ASSAYS = ['open_field', 'novel_object']

VIDEOS_ZIP = []
for i in samples_df['movie_name']:
    VIDEOS_ZIP.append(i.strip('.avi'))
SAMPLES_ZIP = samples_df['sample']
ASSAYS_ZIP = samples_df['assay']

# Exclude samples that can't be tracked
## Read in excluded samples
excl_df = pd.read_csv("config/excluded_samples.csv")
## Create 
zip_df = pd.DataFrame(
    {
        'video': VIDEOS_ZIP,
        'sample': SAMPLES_ZIP,
        'assay': ASSAYS_ZIP
    }
)
## Remove rows matching excluded samples
for i, row in excl_df.iterrows():
    target_assay = row['assay']
    target_sample = row['sample']
    zip_df.drop(
        zip_df[
            (zip_df['sample'] == target_sample) & \
            (zip_df['assay'] == target_assay)
            ].index,
            inplace = True
    )

## Extract samples and assays
VIDEOS_ZIP_EX = zip_df['video']
SAMPLES_ZIP_EX = zip_df['sample']
ASSAYS_ZIP_EX = zip_df['assay']

##########################

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