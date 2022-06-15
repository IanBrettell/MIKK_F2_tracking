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

# Get samples with over 85% tracking success

## Read in .csv
ts_df = pd.read_csv('config/tracking_success.csv')
## Filter for samples with < 85% tracking success
ts_df = ts_df.loc[ts_df['prop_success'] < 0.85]
## Remove samples with 
for i, row in ts_df.iterrows():
    target_assay = row['assay']
    target_sample = row['sample']
    target_video = row['video']
    zip_df.drop(
        zip_df[
            (zip_df['sample'] == target_sample) & \
            (zip_df['assay'] == target_assay) & \
            (zip_df['video'] == target_video)
            ].index,
            inplace = True
    )
## Get filtered samples and assays
VIDEOS_ZIP_EX_TRK = zip_df['video']
SAMPLES_ZIP_EX_TRK = zip_df['sample']
ASSAYS_ZIP_EX_TRK = zip_df['assay']
## Multiply lists for combinations with `seconds_interval`
n_intervals = len(config["seconds_interval"])
VIDEOS_ZIP_EX_TRK_INT = list(VIDEOS_ZIP_EX_TRK.values) * n_intervals
SAMPLES_ZIP_EX_TRK_INT = list(SAMPLES_ZIP_EX_TRK.values) * n_intervals
ASSAYS_ZIP_EX_TRK_INT = list(ASSAYS_ZIP_EX_TRK.values) * n_intervals
## Multiply each element of `seconds_intervals` by length of original SAMPLES list
INTERVALS_ZIP_EX_TRK_INT = np.repeat(config["seconds_interval"], len(SAMPLES_ZIP_EX_TRK_INT))



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


