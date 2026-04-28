# Docker container based on a minimal Ubuntu installation that includes conda-forge's mambaforge installer.
container: "docker://condaforge/mambaforge"

import pandas as pd
from snakemake.utils import validate, min_version
import os
srcdir = workflow.basedir

##### set minimum snakemake version #####
min_version("5.1.2")


##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="06_Schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index(["project", "condition", "sample"], drop=False)
validate(samples, schema="06_Schemas/samples.schema.yaml")

coldata = pd.read_table(config["coldata"]).set_index(["project", "condition", "type"], drop=False)
validate(coldata, schema="06_Schemas/coldata.schema.yaml")

condition = pd.read_table(config["condition"]).set_index(["condition"], drop=False)
validate(coldata, schema="06_Schemas/condition.schema.yaml")

##### Set variables ####
ROOTDIR = os.getcwd()
RAWDATA = os.path.join("/shared/projects/tp_2626_biological_data_183961/data/")
REF = os.path.join("/shared/projects/tp_2626_biological_data_183961/ref/")
CONTAINER = os.path.join(ROOTDIR, "02_Container") + "/"
SCRIPTDIR = os.path.join(ROOTDIR, "03_Script") + "/"
ENVDIR = os.path.join(ROOTDIR, "04_Workflow") + "/"
OUTPUTDIR = os.path.join(ROOTDIR, "05_Output") + "/"
REPORT = os.path.join(ROOTDIR, "07_Report") + "/"

# ----------------------------------------------
# Target rules
# ----------------------------------------------
# SAMPLES = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
SAMPLES = expand("{samples.sample}",samples=samples.itertuples())
RUN =  config["run"]["type"].split(',')
EXT = config["run"]["ext"]
ref_level = config["diffexp"]["ref_level"]
genome = config["ref"]["genome"]
index = config["ref"]["index"]
annotation = config["ref"]["annotation"]
# RUN_ID = expand("{samples.project}_{samples.condition}_{samples.sample}",samples=samples.itertuples())
RUN_ID = expand("{samples.sample}",samples=samples.itertuples())

rule all:
  input:
    ## clean.smk ##
    # expand( OUTPUTDIR + "03_fastqc/{samples}_{run}.trimmed_fastqc.html", samples=SAMPLES, run=RUN),
    ## count.smk ##
    expand(OUTPUTDIR + "06_featurecounts/{samples}_count.txt", samples=SAMPLES)
    ## diffexp.smk ##
    # html_report = OUTPUTDIR + "09_differential_expression/diffexp.html",

# ----------------------------------------------
# setup report
# ----------------------------------------------

report: "report/workflow.rst"

# ----------------------------------------------
# Load rules 
# ----------------------------------------------

include: ENVDIR + "clean.smk"
include: ENVDIR + "count.smk"
include: ENVDIR + "diffexp.smk"
