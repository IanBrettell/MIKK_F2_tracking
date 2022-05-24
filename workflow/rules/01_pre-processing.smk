def get_fps(wildcards):
    fps = samples_df.loc[samples_df["sample"] == wildcards.sample, "fps"]
    return(fps)

# Create copies recoded with opencv to create smaller videos
# That can be viewed with idtrackerai
rule recode_videos:
    input:
        rules.copy_videos.output,
    output:
        os.path.join(
            config["working_dir"],
            "recoded/{video}.avi"
        ),        
    log:
        os.path.join(
            config["working_dir"],
            "logs/recode_videos/{video}.log"
        ),
    container:
        config["opencv"]
    resources:
        mem_mb = 1000
    script:
        "../scripts/recode_videos.py"

## Generate single-frame grab showing coordinate of splits
#rule set_split_coords:
#    input:
#        video = rules.copy_videos.output,
#    output:
#        fig = "results/split_coord_images/{assay}/{movie}.png",
#    log:
#        os.path.join(
#            config["working_dir"],
#            "logs/set_split_coords/{assay}/{movie}.log")
#    params:
#        assay = "{assay}",
#        sample = "{movie}",
#        samples_file = lambda wildcards: config["samples_file"]
#    container:
#        config["opencv"]
#    resources:
#        mem_mb = 500
#    script:
#        "../scripts/set_split_coords.py"
#
## Split videos into quadrants and assays (1 raw video * 4 quadrants * 2 assays = 8 output videos)
#rule split_videos:
#    input:
#        rules.recode_videos.output,
#    output:
#        os.path.join(
#            config["working_dir"],
#            "split/{assay}/{sample}_{quadrant}.avi"
#        ),
#    log:
#        os.path.join(
#            config["working_dir"],
#            "logs/split_videos/{assay}/{sample}/{quadrant}.log"
#        ),
#    params:
#        assay = "{assay}",
#        sample = "{sample}",
#        quadrant = "{quadrant}",
#        samples_file = lambda wildcards: config["samples_file"]
#    container:
#        config["opencv"]
#    resources:
#        mem_mb = 1000
#    script:
#        "../scripts/split_videos.py"


