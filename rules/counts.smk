if "counts" in config:

    # this rule only applies if counts is defined in the config file
    rule load_counts:
        input:
            counts=config["counts"]["path"],
            cells=config["cells"],
        output:
            "analysis/all.rds",
        params:
            spike_pattern=config["spike-ins"]["pattern"],
            species=config["species"],
            feature_ids=config["counts"]["feature_ids"],
        log:
            "logs/load-counts.log",
        #conda:
        #    "../envs/eval.yaml"
        singularity:
            "/share/work/HPC/work_tmp/liangyong/github/Single_Cell_Pipeline_update2022/singularity_test/recipe/eval.sif"        

        script:
            "../scripts/load-counts.R"
