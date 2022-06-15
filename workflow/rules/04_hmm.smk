rule movement_metrics:
    input:
        rules.assign_ref_test.output,
    output:
        os.path.join(
            config["working_dir"],
            "with_metrics/{assay}/{video}/{sample}/{interval}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/movement_metrics/{assay}/{video}/{sample}/{interval}.log"
        ),
    params:
        seconds_interval = "{interval}",
        source_file = "workflow/scripts/movement_metrics_source.R",
    resources:
        mem_mb = 500
    container:
        config["rocker_tidyverse"]
    script:
        "../scripts/movement_metrics.R"

rule merge_csvs:
    input:
        expand(os.path.join(
            config["working_dir"],
            "with_metrics/{assay}/{video}/{sample}/{{interval}}.csv"
            ),
                zip,
                assay = ASSAYS_ZIP_EX_TRK,
                video = VIDEOS_ZIP_EX_TRK,
                sample = SAMPLES_ZIP_EX_TRK,
        ),
    output:
        os.path.join(
            config["working_dir"],
            "merged/{interval}.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/merge_csvs/{interval}.log"
        ),
    resources:
        # start at 15000
        mem_mb = lambda wildcards, attempt: attempt * 15000
    container:
        config["rocker_tidyverse"]
    script:
        "../scripts/merge_csvs.R"

rule send_to_shared_folder:
    input:
        rules.merge_csvs.output,
    output:
        os.path.join(
            config['shared_dir'],
            "{interval}/F2.csv"
        ),
    log:
        os.path.join(
            config["working_dir"],
            "logs/send_to_shared_folder/{interval}.log"
        ),
    resources:
        mem_mb = 200
    shell:
        """
        cp {input} {output}
        """