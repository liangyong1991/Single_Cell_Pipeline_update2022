# Copyright 2018 Johannes Köster.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule cell_cycle:
    input:
        "analysis/filtered.rds",
    output:
        "analysis/cell-cycle-assignments.rds",
    params:
        species=config["species"],
    log:
        "logs/cell-cycle.log",
    #conda:
    #    "../envs/eval.yaml"
    singularity:
            "/share/work/HPC/work_tmp/liangyong/github/Single_Cell_Pipeline_update2022/singularity_test/recipe/eval.sif"
    script:
        "../scripts/cell-cycle.R"


rule cell_cycle_scores:
    input:
        rds="analysis/cell-cycle-assignments.rds",
        cells=config["cells"],
    output:
        report(
            "plots/cycle-scores.{covariate}.pdf",
            caption="../report/cycle-scores.rst",
            category="Quality Control",
        ),
    log:
        "logs/cell-cycle-scores.{covariate}.log",
    #conda:
    #    "../envs/eval.yaml"
    singularity:
            "/share/work/HPC/work_tmp/liangyong/github/Single_Cell_Pipeline_update2022/singularity_test/recipe/eval.sif"
    script:
        "../scripts/cell-cycle-scores.R"
