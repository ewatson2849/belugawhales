import pandas as pd

configfile: "config/assembly.yaml"

samples_df = pd.read_csv("tsv/beluga_raw_reads.tsv", sep="\t")
SAMPLES = samples_df["sample_id"].tolist()
READS = {row.sample_id: {"r1": row.r1, "r2": row.r2} for row in samples_df.itertuples()}

rule all:
    input:
        expand("results/{genera}/fastqc/{sample}/{sample}_R1_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/fastqc/{sample}/{sample}_R2_fastqc.html", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/fastqc/{sample}/{sample}_R1_fastqc.zip", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/fastqc/{sample}/{sample}_R2_fastqc.zip", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/trim_galore/{sample}/{sample}_val_1.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/trim_galore/{sample}/{sample}_val_2.fq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/bowtie_index/{genera}_indexed_ref.1.bt2", genera=config["genera"]),
        expand("results/{genera}/bowtie_index/{genera}_indexed_ref.2.bt2", genera=config["genera"]),
        expand("results/{genera}/bowtie_index/{genera}_indexed_ref.3.bt2", genera=config["genera"]),
        expand("results/{genera}/bowtie_index/{genera}_indexed_ref.4.bt2", genera=config["genera"]),
        expand("results/{genera}/bowtie_index/{genera}_indexed_ref.rev.1.bt2", genera=config["genera"]),
        expand("results/{genera}/bowtie_index/{genera}_indexed_ref.rev.2.bt2", genera=config["genera"]),
        expand("results/{genera}/bowtie_align/{sample}/{sample}_host_removed_R1.fastq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/bowtie_align/{sample}/{sample}_host_removed_R2.fastq.gz", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/spades/{sample}/contigs.fasta", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/virsorter/{sample}/final-viral-combined.fa", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/virsorter/{sample}/final-viral-score.tsv", sample=SAMPLES, genera=config["genera"]),
        expand("results/{genera}/pharokka/{sample}/pharokka.gff", sample=SAMPLES, genera=config["genera"]),
        #expand("results/{genera}/pharokka/{sample}/pharokka_plot.png", sample=SAMPLES, genera=config["genera"])
        expand("results/{genera}/vibrant/{sample}/{sample}.faa", sample = SAMPLES, genera = config["genera"])

rule fastqc:
    """
    Run fastqc on beluga raw reads 
    """
    input:
        fwd = lambda wildcards: READS[wildcards.sample]["r1"],
        rev = lambda wildcards: READS[wildcards.sample]["r2"]
    output:
        "results/{genera}/fastqc/{sample}/{sample}_R1_fastqc.html",
        "results/{genera}/fastqc/{sample}/{sample}_R2_fastqc.html",
        "results/{genera}/fastqc/{sample}/{sample}_R1_fastqc.zip",
        "results/{genera}/fastqc/{sample}/{sample}_R2_fastqc.zip"
    params:
        threads = 4, 
        genera = config["genera"]
    shell:
        """
        module unload miniconda
        module load FastQC/0.12.1-Java-11
        fastqc {input.fwd} {input.rev} -t {threads} -o results/{params.genera}/fastqc/{wildcards.sample}
        mv results/{params.genera}/fastqc/{wildcards.sample}/*R1*.html results/{params.genera}/fastqc/{wildcards.sample}/{wildcards.sample}_R1_fastqc.html
        mv results/{params.genera}/fastqc/{wildcards.sample}/*R2*.html results/{params.genera}/fastqc/{wildcards.sample}/{wildcards.sample}_R2_fastqc.html
        mv results/{params.genera}/fastqc/{wildcards.sample}/*R1*.zip results/{params.genera}/fastqc/{wildcards.sample}/{wildcards.sample}_R1_fastqc.zip
        mv results/{params.genera}/fastqc/{wildcards.sample}/*R2*.zip results/{params.genera}/fastqc/{wildcards.sample}/{wildcards.sample}_R2_fastqc.zip
        """

rule qc:
    """
    Run trim galore on beluga raw reads
    """
    input:
        fwd = lambda wildcards: READS[wildcards.sample]["r1"],
        rev = lambda wildcards: READS[wildcards.sample]["r2"]
    output:
        "results/{genera}/trim_galore/{sample}/{sample}_val_1.fq.gz",
        "results/{genera}/trim_galore/{sample}/{sample}_val_2.fq.gz"
    params:
        genera = config["genera"]
    shell:
        """
        module unload miniconda
        module load Trim_Galore/0.6.7-GCCcore-10.2.0
        trim_galore --paired --fastqc --length 125 --cores 4 {input.fwd} {input.rev} --output_dir results/{params.genera}/trim_galore/{wildcards.sample} --basename {wildcards.sample}
        """

rule bowtie_build:
    """
    Create index of beluga whale genome
    """
    input:
        beluga_ref = "input/beluga_genome_sequence.fna"
    output:
        "results/{genera}/bowtie_index/{genera}_indexed_ref.1.bt2",
        "results/{genera}/bowtie_index/{genera}_indexed_ref.2.bt2",
        "results/{genera}/bowtie_index/{genera}_indexed_ref.3.bt2",
        "results/{genera}/bowtie_index/{genera}_indexed_ref.4.bt2",
        "results/{genera}/bowtie_index/{genera}_indexed_ref.rev.1.bt2",
        "results/{genera}/bowtie_index/{genera}_indexed_ref.rev.2.bt2"
    conda:
        "envs/read_aln.yaml"
    params:
        genera=config["genera"]
    log:
        "results/logs/bowtie_index/{genera}_bowtie_index.log"
    shell:
        """
        bowtie2-build -f {input.beluga_ref} results/{params.genera}/bowtie_index/{params.genera}_indexed_ref > {log}
        """

rule bowtie_align:
    """
    Filter out beluga fecal raw reads that map back to beluga reference genome
    """
    input:
        r1 = "results/{genera}/trim_galore/{sample}/{sample}_val_1.fq.gz",
        r2 = "results/{genera}/trim_galore/{sample}/{sample}_val_2.fq.gz",
        idx = [
            "results/{genera}/bowtie_index/{genera}_indexed_ref.1.bt2",
            "results/{genera}/bowtie_index/{genera}_indexed_ref.2.bt2",
            "results/{genera}/bowtie_index/{genera}_indexed_ref.3.bt2",
            "results/{genera}/bowtie_index/{genera}_indexed_ref.4.bt2",
            "results/{genera}/bowtie_index/{genera}_indexed_ref.rev.1.bt2",
            "results/{genera}/bowtie_index/{genera}_indexed_ref.rev.2.bt2"
        ]
    output:
        "results/{genera}/bowtie_align/{sample}/{sample}_host_removed_R1.fastq.gz",
        "results/{genera}/bowtie_align/{sample}/{sample}_host_removed_R2.fastq.gz"
    conda:
        "envs/read_aln.yaml"
    params:
        genera=config["genera"]
    shell:
        """
        bowtie2 -p 4 -x results/{wildcards.genera}/bowtie_index/{wildcards.genera}_indexed_ref -1 {input.r1} -2 {input.r2} \
        --very-sensitive-local --un-conc-gz results/{wildcards.genera}/bowtie_align/{wildcards.sample}/{wildcards.sample}_host_removed > results/{wildcards.genera}/bowtie_align/{wildcards.sample}/{wildcards.sample}_mapped_and_unmapped.sam
        mv results/{wildcards.genera}/bowtie_align/{wildcards.sample}/{wildcards.sample}_host_removed.1 results/{wildcards.genera}/bowtie_align/{wildcards.sample}/{wildcards.sample}_host_removed_R1.fastq.gz
        mv results/{wildcards.genera}/bowtie_align/{wildcards.sample}/{wildcards.sample}_host_removed.2 results/{wildcards.genera}/bowtie_align/{wildcards.sample}/{wildcards.sample}_host_removed_R2.fastq.gz
        """

rule spades:
    """
    Assemble reads with SPAdes
    """
    input:
        fwd = "results/{genera}/bowtie_align/{sample}/{sample}_host_removed_R1.fastq.gz",
        rev = "results/{genera}/bowtie_align/{sample}/{sample}_host_removed_R2.fastq.gz"
    output:
        "results/{genera}/spades/{sample}/contigs.fasta"
    params:
        genera=config["genera"]
    resources:
        mem_mb=1000000,
        threads=4
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/spades
        spades.py --meta --threads {resources.threads} -1 {input.fwd} -2 {input.rev} -o results/{params.genera}/spades/{wildcards.sample}
        """

rule virsorter:
    """
    Run VirSorter on SPAdes assembly
    """
    input:
        "results/{genera}/spades/{sample}/contigs.fasta"
    output:
        "results/{genera}/virsorter/{sample}/final-viral-combined.fa",
        "results/{genera}/virsorter/{sample}/final-viral-score.tsv"
    params:
        genera=config["genera"]
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/vs2
        virsorter run -w results/{params.genera}/virsorter/{wildcards.sample} -i {input} --min-length 1500 -j 4 all
        """

rule pharokka:
    """
    Annotate phage genome
    """
    input:
        "input/NODE_59.fasta"
    output:
        "results/{genera}/pharokka/{sample}/pharokka.gff"
        #"results/{genera}/pharokka/{sample}/pharokka_plot.png"
    params:
        genera=config["genera"],
        outdir = "results/{genera}/pharokka/{sample}"
    shell:
        """
        module unload miniconda
        source activate /home/flg9/.conda/envs/pharokka_env
        mkdir -p {params.outdir}
        pharokka.py -i {input} -o {params.outdir} -d pharokka_db -t 4 -p pharokka -f
        #pharokka_plotter.py -i {input} -n pharokka_plot -o {params.outdir} -p pharokka -f
        """

#pharokka_plotter.py -i pharokka_plot/NODE_59.fasta -n pharokka_plot -o pharokka_plot -p pharokka -t 'Beluga ${Clostridium}$ Phage NODE_58' --interval 5000 --truncate 40 --label_size 10 --force

rule VIBRANT:
    """
    Cross-verify phages Virsorter identified with VIBRANT
    """
    input:
        "results/{genera}/spades/{sample}/contigs.fasta"
    output:
        "results/{genera}/vibrant/{sample}/{sample}.faa"
    params:
        outdir = "results/{genera}/vibrant/{sample}/"
    shell:
        """
        module unload miniconda 
        source activate /home/flg9/.conda/envs/VIBRANT
        mkdir -p {params.outdir}
        python3 VIBRANT/VIBRANT_run.py -i {input} -folder {params.outdir} -t 4 -d /vast/palmer/pi/turner/flg9/turner_lab/VIBRANT/databases/ -m /vast/palmer/pi/turner/flg9/turner_lab/VIBRANT/files/
        """