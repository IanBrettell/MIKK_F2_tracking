# Copy videos from FTP to Codon (for ease downstream)
rule copy_videos:
    input:
        os.path.join(config["raw_data_dir"], "{sample}.avi"),
    output:
        os.path.join(config["working_dir"], "raw_videos/{sample}.avi"),
    log:
        os.path.join(config["working_dir"], "logs/copy_videos/{sample}.log"),
    shell:
        """
        cp {input} {output} \
            2> {log}
        """

# Decompress videos so that they can be viewed in Fiji
# in order to determine the start and end frames for open_field and novel_object assays
# NOTE: creates huge files, so consider hashing out the rule and deleting the files once
# metadata has been collected
rule decompress_videos:
    input:
        rules.copy_videos.output,
    output:
        os.path.join(config["working_dir"], "decompressed/{sample}.avi"),
    log:
        os.path.join(config["working_dir"], "logs/decompress_videos/{sample}.log"),
    container:
        config["ffmpeg"]
    shell:
        """
        ffmpeg -i {input} -an -vcodec rawvideo -y {output} \
            2> {log}
        """

# Generate single-frame grab showing coordinate of splits
rule set_split_coords:
    input:
        video = rules.copy_videos.output,
    output:
        fig = "results/split_coord_images/{assay}/{sample}.png",
    log:
        os.path.join(config["working_dir"], "logs/set_split_coords/{assay}/{sample}.log")
    params:
        sample = "{sample}",
        assay = "{assay}",
        samples_file = lambda wildcards: config["samples_file"]
    container:
        config["opencv"]
    script:
        "../scripts/set_split_coords.py"

# Split videos into quadrants and assays (1 raw video * 4 quadrants * 2 assays = 8 output videos)
rule split_videos:
    input:
        rules.copy_videos.output,
    output:
        os.path.join(config["working_dir"], "split/{assay}/{sample}_{quadrant}.avi"),
    log:
        os.path.join(config["working_dir"], "logs/split_videos/{assay}/{sample}/{quadrant}.log"),
    params:
        assay = "{assay}",
        sample = "{sample}",
        quadrant = "{quadrant}",
        samples_file = lambda wildcards: config["samples_file"]
    container:
        config["opencv"]
    script:
        "../scripts/split_videos.py"


