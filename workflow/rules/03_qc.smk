# Choose input trajectories .csv file
def get_final_csvs(wildcards):
    # Get path of csv files
    ## Trajectories without gaps
    traj_wo_gaps_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        wildcards.video,
        "session_" + wildcards.sample,
        "trajectories_wo_gaps",
        "trajectories_wo_gaps.trajectories.csv"
        )
    ## Trajectories (with gaps)
    traj_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        wildcards.video,
        "session_" + wildcards.sample,
        "trajectories",
        "trajectories.trajectories.csv"
        )
    # If there is no "without gaps" file, return the "trajectories" file
    if os.path.exists(traj_wo_gaps_file):
        return(traj_wo_gaps_file)
    else:
        return(traj_file)

# Get frames-per-second
def get_fps(wildcards):
    fps = samples_df.loc[
        (samples_df['sample'] == wildcards.sample) & \
        (samples_df['assay'] == wildcards.assay),
        "fps"].values[0]
    return(str(fps))

#samples_df.loc[
#    (samples_df['sample'] == SAMPLE) & \
#    (samples_df['assay'] == ASSAY),
#    "fps"].values[0]

# Get relative location of reference iCab fish
def get_ref_loc(wildcards):
    cab_loc = samples_df.loc[(samples_df['sample'] == wildcards.sample) & \
                             (samples_df['assay'] == wildcards.assay), \
                             'cab_coords'].values[0]
    # `ref_loc` is nan if there are two iCabs, so convert to string
    if pd.isna(cab_loc):
        cab_loc = "NA"
    return(cab_loc)

#samples_df.loc[(samples_df['sample'] == SAMPLE) & \
#                             (samples_df['assay'] == ASSAY), \
#                             'cab_coords'].values[0]

# Assign reference and test fish IDs, and filter for frames up to 10 minutes
rule assign_ref_test:
    input:
        get_final_csvs,
    output:
        os.path.join(
            config["working_dir"],
            "final_tracks/{assay}/{video}/{sample}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/assign_ref_test/{assay}/{video}/{sample}.log"
        ),
    params:
        fps = get_fps,
        ref_loc = get_ref_loc,
    resources:
        mem_mb = 200
    script:
        "../scripts/assign_ref_test.py"
    
# Create .csv with proportion of frames tracked
rule tracking_success:
    input:
        expand(rules.assign_ref_test.output,
            zip,
            assay = ASSAYS_ZIP_EX,
            video = VIDEOS_ZIP_EX,
            sample = SAMPLES_ZIP_EX                   
        ),
    output:
        "config/tracking_success.csv"
    log:
        os.path.join(
            config["working_dir"],
            "logs/tracking_success/tracking_success.log"
        ),
    container:
        config["rocker_tidyverse"]
    resources:
        mem_mb = 1000
    script:
        "../scripts/tracking_success.R"

# Function to pull trajectories.npy if trajectories_wo_gaps.npy does not exist
def get_trajectories_file(wildcards):
    # Get path of csv files
    ## Trajectories without gaps
    traj_wo_gaps_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        wildcards.video,
        "session_" + wildcards.sample,
        "trajectories_wo_gaps",
        "trajectories_wo_gaps.npy"
        )
    ## Trajectories (with gaps)
    traj_file = os.path.join(
        config["working_dir"],
        "split",
        wildcards.assay,
        wildcards.video,
        "session_" + wildcards.sample,
        "trajectories",
        "trajectories.npy"
        )
    # If there is no "without gaps" file, return the "trajectories" file
    if os.path.exists(traj_wo_gaps_file):
        return(traj_wo_gaps_file)
    elif os.path.exists(traj_file):
        return(traj_file)

