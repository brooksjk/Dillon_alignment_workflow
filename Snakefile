configfile: "config.yaml"
import os
import re

# Folder containing the files
folder_path = config["SAMPLES"]

# Regular expression pattern to extract sample names up to the last "_1" or "_2"
pattern = r'^(.+?)(?:_[12])?(?![\d\w])'

# Store unique sample names
sample_names = set()

# Iterate over files in the folder
for file_name in os.listdir(folder_path):
    # Full path of the file
    file_path = os.path.join(folder_path, file_name)
    
    # Check if it's a file
    if os.path.isfile(file_path):
        # Extract sample name using regular expression
        match = re.match(pattern, file_name)
        if match:
            sample_names.add(match.group(1))

# Convert set to list if needed
sample_names_list = list(sample_names)

print(sample_names_list)

rule all:
    input:
        expand(config["WORKDIR"] + "/{sample}_bwa_wb.dup.bam", sample=sample_names_list)

# Rule for processing ancient samples
rule process_samples:
    input:
        software=config["ADAPTREMOVE"],
        ref=config["REF"]
    output:
        bam=config["WORKDIR"] + "/{sample}_bwa_wb.bam",
        sorted_bam=config["WORKDIR"] + "/{sample}_bwa_wb.srt.bam",
        dedup_bam=config["WORKDIR"] + "/{sample}_bwa_wb.dup.bam",
        dedup_bam_index=config["WORKDIR"] + "/{sample}_bwa_wb.dup.bam.bai"
    params:
        reads=config["SAMPLES"] + "/{sample}",
        sample=config["WORKDIR"] + "/{sample}"
    resources: cpus=40,mem_mb=40000,time_min=10080
    shell:
        """
        module load bwa
        module load samtools

        if [[ "{wildcards.sample}" =~ ^Ancient_ ]]; then
            {input.software} --collapse --trimns --trimqualities --mm 5 --minlength 25 --file1 {params.reads}_1.fastq.gz --file2 {params.reads}_1.fastq.gz --basename {wildcards.sample}_trimmed --threads 40
            bwa aln -l 1024 -n 0.01 -o 2 -t 40 {input.ref} {params.sample}_trimmed.pair1.truncated > {wildcards.sample}_bwa_wb_1.sai
            bwa aln -l 1024 -n 0.01 -o 2 -t 40 {input.ref} {params.sample}_trimmed.pair2.truncated > {wildcards.sample}_bwa_wb_2.sai
            bwa sampe {input.ref} {wildcards.sample}_bwa_wb_1.sai {wildcards.sample}_bwa_wb_2.sai {params.sample}_trimmed.pair1.truncated {params.sample}_trimmed.pair2.truncated | samtools view -b -S -q 25 -F 4 > {output.bam}
            samtools sort {output.bam} -o {output.sorted_bam}
            samtools rmdup -s {output.sorted_bam} {output.dedup_bam}
            samtools index {output.dedup_bam}
        elif [[ "{wildcards.sample}" =~ ^Modern_ ]]; then
            {input.software} --file1 {params.reads}_1.fastq.gz --file2 {params.reads}_1.fastq.gz --basename {wildcards.sample}_trimmed
            bwa mem {input.ref} {params.sample}_trimmed.pair1.truncated {params.sample}_trimmed.pair2.truncated > {wildcards.sample}_bwa_wb.sam
            samtools view -b -q 30 -F 4 {wildcards.sample}_bwa_wb.sam > {output.bam}
            samtools sort {output.bam} -o {output.sorted_bam}
            samtools rmdup -s {output.sorted_bam} {output.dedup_bam}
            samtools index {output.dedup_bam}
        elif [[ "{wildcards.sample}" =~ ^SE_Ancient_ ]]; then
            {input.software} --file1 {params.reads}_*.fastq.gz --basename {wildcards.sample}_trimmed
            bwa aln -l 1024 -n 0.01 -o 2 -t 4 {input.ref} {params.sample}_trimmed.truncated > {wildcards.sample}_bwa_wb.sai
            bwa samse {input.ref} {wildcards.sample}_bwa_wb.sai {wildcards.sample}_trimmed.truncated | samtools view -b -S -q 25 -F 4 > {output.bam}
            samtools sort {output.bam} -o {output.sorted_bam}
            samtools rmdup -s {output.sorted_bam}  {output.dedup_bam}
            samtools index {output.dedup_bam}
        else
            echo "No sample tag provided."
        fi
        """
