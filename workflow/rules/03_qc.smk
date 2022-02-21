# Choose input trajectories .csv file
def get_final_csvs(wildcards):
    # Get path of csv files
    traj_wo_gaps_file = os.path.join(
        config["working_dir"],
        "split/{assay}/session_{sample}_{quadrant}/trajectories_wo_gaps/trajectories_wo_gaps.trajectories.csv"
        )
    traj_file = os.path.join(
        config["working_dir"],
        "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.trajectories.csv"
        )
    # If there is no "without gaps" file, return the "trajectories" file
    if os.path.exists(traj_wo_gaps_file):
        return(traj_wo_gaps_file)
    else:
        return(traj_file)

# Get frames-per-second
def get_fps(wildcards):
    fps = int(samples_df.loc[samples_df["sample"] == wildcards.sample, "fps"])
    return(fps)

# Get relative location of reference iCab fish
def get_ref_loc(wildcards):
    if wildcards.assay == "open_field":
        target_col = "cab_coords_" + "of_" + wildcards.quadrant
        ref_loc = samples_df.loc[samples_df["sample"] == wildcards.sample, target_col]
    elif wildcards.assay == "novel_object":
        target_col = "cab_coords_" + "no_" + wildcards.quadrant
        ref_loc = samples_df.loc[samples_df["sample"] == wildcards.sample, target_col]
    return(ref_loc)

# Assign reference and test fish IDs, and filter for frames up to 10 minutes
rule assign_ref_test:
    input:
        get_final_csvs,
    output:
        os.path.join(
            config["data_store_dir"],
            "final_tracks/{assay}/{sample}_{quadrant}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/assign_ref_test/{assay}/{sample}/{quadrant}.log"
        ),
    params:
        fps = get_fps,
        ref_loc = get_ref_loc,
    resources:
        mem_mb = 500
    script:
        "../scripts/assign_ref_test.py"
    
    
    

## Create .csv with proportion of frames tracked
#rule calculate_prop_frames_tracked:
#    input:
#        get_final_csvs,
#    output:

