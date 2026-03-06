SAMPLES = expand("{samples.sample}",samples=samples.itertuples())
RUN =  config["run"]["type"].split(',')
EXT = config["run"]["ext"]

# ----------------------------------------------
# Trimmomatic: trimming reads and removing adapter sequences
# ----------------------------------------------

rule trimmomatic:
  input:
    sample=expand(RAWDATA + "{samples}_{run}.{ext}", samples=SAMPLES, run=RUN, ext=EXT)
  output:
    sample_trimmed=expand(OUTPUTDIR + "02_trimmomatic/{samples}_{run}.trimmed.fastq", samples=SAMPLES, run=RUN),
    sample_untrimmed=expand(OUTPUTDIR + "02_trimmomatic/{samples}_{run}un.trimmed.fastq", samples=SAMPLES, run=RUN)
  conda:
    CONTAINER + "trimmomatic.yaml"
  shell:
    """
    sample=({input.sample})
    sample_trimmed=({output.sample_trimmed})
    sample_untrimmed=({output.sample_untrimmed})
    len=${{#sample[@]}}
    for (( i=0; i<$len; i=i+2 ))
    do trimmomatic PE -threads 4 ${{sample[$i]}} ${{sample[$i+1]}} ${{sample_trimmed[$i]}} ${{sample_untrimmed[$i]}} ${{sample_trimmed[$i+1]}} ${{sample_untrimmed[$i+1]}} LEADING:20 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36
    done
    """


# ----------------------------------------------
# FastQC to check the reads trimmed quality
# ----------------------------------------------
rule fastqc_trimmed:
  input:
    expand(OUTPUTDIR + "02_trimmomatic/{samples}_{run}.trimmed.fastq", samples=SAMPLES, run=RUN)
  output:
    expand(OUTPUTDIR + "03_fastqc/{samples}_{run}.trimmed_fastqc.html", samples=SAMPLES, run=RUN),
    expand(OUTPUTDIR + "03_fastqc/{samples}_{run}.trimmed_fastqc.zip", samples=SAMPLES, run=RUN)
  conda:
    CONTAINER + "fastqc.yaml"
  shell:
    "fastqc --outdir {OUTPUTDIR}03_fastqc/ {input}"

