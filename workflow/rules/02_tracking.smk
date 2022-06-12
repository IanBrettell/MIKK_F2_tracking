# pull tracking parameters from config/samples.csv
def get_vid_length(wildcards):
    start = samples_df.loc[
        (samples_df['sample'] == wildcards.sample) & \
        (samples_df['assay'] == wildcards.assay),
        "frame_start"
    ]
    end = samples_df.loc[
        (samples_df['sample'] == wildcards.sample) & \
        (samples_df['assay'] == wildcards.assay),
        "frame_end"
    ]
    vid_length = int(list(end)[0]) - int(list(start)[0])
    return(vid_length)

def get_bgsub(wildcards):
    bgsub = samples_df.loc[
        (samples_df['sample'] == wildcards.sample) & \
        (samples_df['assay'] == wildcards.assay),
        "bgsub"
    ]
    bgsub = list(bgsub)[0]
    return(bgsub)

def get_intensity_floor(wildcards):
    int_floor = samples_df.loc[
        (samples_df['sample'] == wildcards.sample) & \
        (samples_df['assay'] == wildcards.assay),
        "int_floor"
    ]
    int_floor = list(int_floor)[0]
    return(int_floor)

def get_intensity_ceiling(wildcards):
    int_ceiling = samples_df.loc[
        (samples_df['sample'] == wildcards.sample) & \
        (samples_df['assay'] == wildcards.assay),
        "int_ceiling"
    ]
    int_ceiling = list(int_ceiling)[0]
    return(int_ceiling)

def get_area_floor(wildcards):
    area_floor = samples_df.loc[
        (samples_df['sample'] == wildcards.sample) & \
        (samples_df['assay'] == wildcards.assay),
        "area_floor"
    ]
    area_floor = list(area_floor)[0]
    return(area_floor)

def get_area_ceiling(wildcards):
    area_ceiling = samples_df.loc[
        (samples_df['sample'] == wildcards.sample) & \
        (samples_df['assay'] == wildcards.assay),
        "area_ceiling"
    ]
    area_ceiling = list(area_ceiling)[0]
    return(area_ceiling)

# Track videos with idtrackerai
## Note: `output` is set as `trajectories.npy` instead of `trajectories_wo_gaps.npy`, presumably because
## in videos where there are no crossovers, the latter file is not produced.
rule track_videos:
    input:
        rules.split_videos.output
    output:
        os.path.join(
            config["working_dir"], 
            "split/{assay}/{video}/session_{sample}/trajectories/trajectories.npy"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/track_videos/{assay}/{video}/{sample}.log"
        ),
    params:
        vid_length = get_vid_length,
        vid_name = "{sample}",
        bgsub = get_bgsub,
        intensity_floor = get_intensity_floor,
        intensity_ceiling = get_intensity_ceiling,
        area_floor = get_area_floor,
        area_ceiling = get_area_ceiling,
    resources:
        # start at 5000
        mem_mb = lambda wildcards, attempt: attempt * 20000
    container:
        config["idtrackerai"]
    shell:
        """
        idtrackerai terminal_mode \
            --_video {input} \
            --_bgsub '{params.bgsub}' \
            --_range [0,{params.vid_length}] \
            --_nblobs 2 \
            --_intensity [{params.intensity_floor},{params.intensity_ceiling}] \
            --_area [{params.area_floor},{params.area_ceiling}] \
            --_session {params.vid_name} \
            --exec track_video \
                2> {log}
        """


# Convert numpy arrays to .csv files
# Not needed if `CONVERT_TRAJECTORIES_DICT_TO_CSV_AND_JSON = True`
# is added to `local_settings.py`
#rule trajectories_to_csv:
#    input:
#        trajectories = rules.track_videos.output,
#        script = "workflow/scripts/trajectories_to_csv.py"
#    output:
#        os.path.join(config["working_dir"], "split/{assay}/session_{sample}_{quadrant}/trajectories/trajectories.trajectories.csv")
#    log:
#        os.path.join(config["working_dir"], "logs/trajectories_to_csv/{assay}/{sample}/{quadrant}.log"),
#    params:
#        in_path = os.path.join(config["working_dir"], "split/{assay}/session_{sample}_{quadrant}")
#    resources:
#        mem_mb = 100,
#    shell:
#        """
#        python {input.script} {params.in_path} \
#            2> {log}
#        """



