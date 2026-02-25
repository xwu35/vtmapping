localrules:
    rename_raw_reads,
    download_phix_genome,
    download_human_genome

rule rename_raw_reads:
    """rename raw reads to make sure the extension works for fastqc"""
    input:
        R1=lambda wildcards: R1_MAP[wildcards.sample],
        R2=lambda wildcards: R2_MAP[wildcards.sample]
    output:
        R1=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R1.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R2.fastq.gz")
    shell:
        """
        ln -s {input.R1} {output.R1} 
        ln -s {input.R2} {output.R2} 
        """

rule fastqc_raw_reads:
    """run fastqc on raw reads"""
    input:
        R1=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R1.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "renamed_raw_reads", "{sample}_R2.fastq.gz")
    output:
        R1_zip=os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads", "{sample}_R1_fastqc.zip"),
        R2_zip=os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads", "{sample}_R2_fastqc.zip")
    params:
        dir=lambda w, output: os.path.dirname(output.R1_zip)
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        fastqc {input.R1} -o {params.dir} -t {threads}
        fastqc {input.R2} -o {params.dir} -t {threads}
        """

rule multiqc_raw_reads:
    """Aggregate fastqc results from raw reads"""
    input:
        R1_zip=expand(os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads", "{sample}_R1_fastqc.zip"), sample=SAMPLE),
        R2_zip=expand(os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads", "{sample}_R2_fastqc.zip"), sample=SAMPLE)
    output:
        html=os.path.join(dir["output"]["reads_processing"], "fastqc", "multiqc_raw_reads", "multiqc_report.html")
    params:
        in_dir=directory(os.path.join(dir["output"]["reads_processing"], "fastqc", "raw_reads")),
        out_dir=lambda w, output: os.path.dirname(output.html)
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        multiqc {params.in_dir} --outdir {params.out_dir} 
        """

rule trimmomatic:
    """remove adapters and low quality reads"""
    input:
        R1=lambda wildcards: R1_MAP[wildcards.sample],
        R2=lambda wildcards: R2_MAP[wildcards.sample]
    output:
        R1P=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"),
        R1U=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R1.unpaired.fastq.gz"),
        R2P=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz"),
        R2U=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R2.unpaired.fastq.gz")
    params:
        seqType = config["trimmomatic"]["seqType"],
        phred = config["trimmomatic"]["phred"],
        adapter = config["adapter"],
        adapter_params = config["trimmomatic"]["adapter_params"],
        post_adapter_params = config["trimmomatic"]["post_adapter_params"],
    log:
        os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}.trimmomatic.log"),
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["med_mem"]
    conda:
        os.path.join(dir["env"], "trimmomatic.yml")
    shell:
        """
        trimmomatic {params.seqType} {params.phred} -threads {threads} \
            {input.R1} {input.R2} \
            {output.R1P} {output.R1U} {output.R2P} {output.R2U} \
            ILLUMINACLIP:{params.adapter}:{params.adapter_params} \
            {params.post_adapter_params} 2>{log} 
        """

rule fastqc_trimmed_reads:
    """run fastqc on trimmed reads after trimmomatic to check if adapters were removed"""
    input:
        R1P=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"),
        R2P=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz")
    output:
        R1_zip=os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic", "{sample}_R1.trimmomatic_fastqc.zip"),
        R2_zip=os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic", "{sample}_R2.trimmomatic_fastqc.zip")
    params:
        dir=lambda w, output: os.path.dirname(output.R1_zip)
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        fastqc {input.R1P} -o {params.dir} -t {threads}
        fastqc {input.R2P} -o {params.dir} -t {threads}
        """

rule multiqc_trimmed_reads:
    """aggregate fastqc results from trimmed reads after trimmomatic"""
    input:
        R1_zip=expand(os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic", "{sample}_R1.trimmomatic_fastqc.zip"), sample=SAMPLE),
        R2_zip=expand(os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic", "{sample}_R2.trimmomatic_fastqc.zip"), sample=SAMPLE)
    output:
        html=os.path.join(dir["output"]["reads_processing"], "fastqc", "multiqc_after_trimmomatic", "multiqc_report.html")
    params:
        in_dir=directory(os.path.join(dir["output"]["reads_processing"], "fastqc", "after_trimmomatic")),
        out_dir=lambda w, output: os.path.dirname(output.html)
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "qc.yml")
    shell:
        """
        multiqc {params.in_dir} --outdir {params.out_dir} 
        """

rule download_phix_genome:
    """download the phix genome"""
    output:
        os.path.join(dir["db"], "phix_genome", "phix.fasta")
    log:
        os.path.join(dir["db"], "phix_genome", "phix.log")
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
    shell:
        """
        # if the phix.fasta file doesn't exist, then download
        if [[ ! -f "{output}" ]]; then
            curl -L -o {output}.gz -s {params.url} > {log} 2>&1
            gunzip {output}.gz
        fi
        """

if config["mapper"]=="bowtie2":
    rule bowtie2_build_phix_genome:
        """build bowtie2 index for phiX genome"""
        input:
            phix_fasta=os.path.join(dir["db"], "phix_genome", "phix.fasta")
        output:
            multiext(os.path.join(dir["db"], "phix_genome", "phix"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params: 
            prefix=lambda w, input: os.path.splitext(input.phix_fasta)[0]
        threads: 
            config["resources"]["small_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2-build \
                --threads {threads} \
                {input.phix_fasta} \
                {params.prefix} 
            """

    rule phix_genome_mapping_bowtie2:
        """map reads back to phiX genome using bowtie2"""
        input:
            bt_index=multiext(os.path.join(dir["db"], "phix_genome", "phix"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
            R1P=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"),
            R2P=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}.sam"))
        params:
            setting=config["bowtie2"]["extra_settings"],
            prefix=os.path.join(dir["db"], "phix_genome", "phix")
        log:
            os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}.bowtie2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2 {params.setting} -p {threads} \
                -x {params.prefix} -1 {input.R1P} -2 {input.R2P} \
                -S {output.sam} 2> {log} 
            """

elif config["mapper"]=="minimap2":
    rule phix_genome_mapping_minimap2:
        """map reads back to phix genome using minimap2"""
        input:
            phix_fasta=os.path.join(dir["db"], "phix_genome", "phix.fasta"),
            R1=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R1.trimmomatic.fastq.gz"),
            R2=os.path.join(dir["output"]["reads_processing"], "trimmomatic", "{sample}_R2.trimmomatic.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}.sam"))
        params:
            setting=config["minimap2"]["settings"]
        log:
            os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}.minimap2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "coverm.yml")
        shell:
            """
            minimap2 -ax {params.setting} \
                -t {threads} \
                --secondary=no \
                {input.phix_fasta} \
                {input.R1} \
                {input.R2} > {output.sam} 2>{log}
            """

rule phix_reads_removal:
    """extract unmapped paired reads"""
    input:
        sam=os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}.sam")
    output:
        R1=os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}_R1.phixfilt.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}_R2.phixfilt.fastq.gz")
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools view -bS -f 12 {input.sam} | \
        samtools sort -n -@ {threads} -T /tmp -o - | \
        samtools fastq \
            -1 {output.R1} \
            -2 {output.R2} \
            -0 /dev/null \
            -s /dev/null \
            -@ {threads} \
            -N -  
        """

rule download_human_genome:
    """download the human genome"""
    output:
        os.path.join(dir["db"], "human_genome", "GRCh38.fasta")
    log:
        os.path.join(dir["db"], "human_genome", "GRCh38.log")
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    shell:
        """
        # if the GRCh38.fasta file doesn't exist, then download
        if [[ ! -f "{output}" ]]; then
            curl -L -o {output}.gz -s {params.url} > {log} 2>&1
            gunzip {output}.gz
        fi
        """

if config["mapper"]=="bowtie2":
    rule bowtie2_build_human_genome:
        """build bowtie2 index for human genome"""
        input:
            human_fasta=os.path.join(dir["db"], "human_genome", "GRCh38.fasta")
        output:
            multiext(os.path.join(dir["db"], "human_genome", "GRCh38"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params: 
            prefix=lambda w, input: os.path.splitext(input.human_fasta)[0]
        threads: 
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2-build \
                --threads {threads} \
                {input.human_fasta} \
                {params.prefix} 
            """

    rule human_genome_mapping_bowtie2:
        """map reads back to human genome using bowtie2"""
        input:
            bt_index=multiext(os.path.join(dir["db"], "human_genome", "GRCh38"), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
            R1=os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}_R1.phixfilt.fastq.gz"),
            R2=os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}_R2.phixfilt.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}.sam"))
        params:
            setting=config["bowtie2"]["extra_settings"],
            prefix=os.path.join(dir["db"], "human_genome", "GRCh38")
        log:
            os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}.bowtie2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2 {params.setting} -p {threads} \
                -x {params.prefix} \
                -1 {input.R1} -2 {input.R2} \
                -S {output.sam} 2> {log} 
            """

elif config["mapper"]=="minimap2":
    rule human_genome_mapping_minimap2:
        """map reads back to human genome using minimap2"""
        input:
            human_fasta=os.path.join(dir["db"], "human_genome", "GRCh38.fasta"),
            R1=os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}_R1.phixfilt.fastq.gz"),
            R2=os.path.join(dir["output"]["reads_processing"], "phix_filtered", "{sample}_R2.phixfilt.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}.sam"))
        params:
            setting=config["minimap2"]["settings"]
        log:
            os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}.minimap2.log")
        threads:
            config["resources"]["med_cpu"]
        resources:
            mem_mb=config["resources"]["small_mem"]
        conda:
            os.path.join(dir["env"], "coverm.yml")
        shell:
            """
            minimap2 -ax {params.setting} \
                -t {threads} \
                --secondary=no \
                {input.human_fasta} \
                {input.R1} \
                {input.R2} > {output.sam} 2>{log}
            """

rule human_reads_removal:
    """extract unmapped paired reads"""
    input:
        sam=os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}.sam")
    output:
        R1=os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}_R1.humanfilt.fastq.gz"),
        R2=os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}_R2.humanfilt.fastq.gz")
    threads:
        config["resources"]["med_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools view -bS -f 12 {input.sam} | \
        samtools sort -n -@ {threads} -T /tmp -o - | \
        samtools fastq \
            -1 {output.R1} \
            -2 {output.R2} \
            -0 /dev/null \
            -s /dev/null \
            -@ {threads} \
            -N -
        """

rule count_reads_after_human_reads_removal:
    """check read counts after removing human contaminated reads"""
    input:
        R1=expand(os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}_R1.humanfilt.fastq.gz"), sample=SAMPLE)
    output:
        R1_stats=os.path.join(dir["output"]["reads_statistics"], "after_human_reads_removal", "R1_stats.tsv")
    conda:
        os.path.join(dir["env"], "seqkit.yml")
    threads:
        config["resources"]["small_cpu"]
    shell:
        """
        seqkit stats -j {threads} -To {output.R1_stats} {input.R1}
        """

if config["mapper"]=="bowtie2":
    rule bowtie2_build_host_genome:
        """build bowtie2 index for host genome"""
        input:
            ref=config["host_genome"]
        output:
            multiext(config["host_genome"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params: 
            prefix=config["host_genome"]
        threads: 
            config["resources"]["small_cpu"]
        conda:
           os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2-build \
                --threads {threads} \
                {input.ref} \
                {params.prefix} 
            """

    rule host_genome_mapping_bowtie2:
        """map reads back to host genome using bowtie2"""
        input:
            bt_index=multiext(config["host_genome"], ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
            R1=os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}_R1.humanfilt.fastq.gz"),
            R2=os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}_R2.humanfilt.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}.sam"))
        params:
            setting=config["bowtie2"]["extra_settings"],
            prefix=config["host_genome"]
        log:
            os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}.bowtie2.log")
        threads:
            config["resources"]["med_cpu"]
        conda:
            os.path.join(dir["env"], "bowtie2.yml")
        shell:
            """
            bowtie2 {params.setting} -p {threads} \
                -x {params.prefix} \
                -1 {input.R1} -2 {input.R2} \
                -S {output.sam} 2> {log} 
            """
            
elif config["mapper"]=="minimap2":
    rule host_genome_mapping_minimap2:
        """map reads back to host genome using minimap2"""
        input:
            R1=os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}_R1.humanfilt.fastq.gz"),
            R2=os.path.join(dir["output"]["reads_processing"], "human_filtered", "{sample}_R2.humanfilt.fastq.gz")
        output:
            sam=temp(os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}.sam"))
        params:
            setting=config["minimap2"]["settings"],
            host_fasta=config["host_genome"]
        log:
            os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}.minmap2.log")
        threads:
            config["resources"]["med_cpu"]
        resources:
            mem_mb=config["resources"]["small_mem"]
        conda:
            os.path.join(dir["env"], "coverm.yml")
        shell:
            """
            minimap2 -ax {params.setting} \
                -t {threads} \
                --secondary=no {params.host_fasta} \
                {input.R1} {input.R2} > {output.sam} 2> {log}
            """

rule host_separate_mapped_and_umapped_reads:
    """separate mapped and umapped reads"""
    input:
        sam=os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}.sam")
    output:
        mapped=temp(os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}_mapped_sorted.bam")),
        unmapped=temp(os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}_unmapped_sorted.bam"))
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        # separate mapped and unmapped reads and sort the bam file by reference coordinates 
        samtools view -bS -F 4 {input.sam} | samtools sort - -o {output.mapped} -@ {threads} -T /tmp 
        samtools view -bS -f 4 {input.sam} | samtools sort - -o {output.unmapped} -@ {threads} -T /tmp 
        """

rule host_filter_mapped:
    """remove mapped reads with >= the given identity"""
    input:
        mapped=os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}_mapped_sorted.bam"),
    output:
        filtered=temp(os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}_mapped_sorted_filtered.bam"))
    threads:
        config["resources"]["med_cpu"]
    params:
        pid=config["identity"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """                    
        coverm filter -b {input.mapped} \
            -o {output.filtered} \
            --inverse \
            --min-read-percent-identity {params.pid} \
            --threads {threads}
        """

rule host_merge_bam_files:
    """merge unmapped and mapped filtered bam files"""
    input:
        unmapped=os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}_unmapped_sorted.bam"),
        mapped_filtered=os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}_mapped_sorted_filtered.bam")
    output:
        kept=os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}_kept.bam")
    threads:
        config["resources"]["small_cpu"]
    resources:
        mem_mb=config["resources"]["small_mem"]
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """
        samtools merge \
            -o {output.kept} \
            -@ {threads} \
            {input.unmapped} \
            {input.mapped_filtered} 
        """

rule count_kept_reads:
    input:
        kept=expand(os.path.join(dir["output"]["reads_processing"], "host_filtered", "{sample}_kept.bam"), sample=SAMPLE)
    output:
        os.path.join(dir["output"]["reads_statistics"], "after_host_reads_removal", "read_counts_kept.txt")
    params:
        dir=os.path.join(dir["output"]["reads_processing"], "host_filtered"),
        script=os.path.join(dir["scripts"], "count_primary_alignments.sh")
    conda:
        os.path.join(dir["env"], "coverm.yml")
    shell:
        """                    
        {params.script} -i {params.dir} -e _kept.bam -o {output}
        """

rule plot_reads_proportion:
    """combine reads counts and calculate the precentage of unmapped reads"""
    input:
        after_human=os.path.join(dir["output"]["reads_statistics"], "after_human_reads_removal", "R1_stats.tsv"),
        after_host=os.path.join(dir["output"]["reads_statistics"], "after_host_reads_removal", "read_counts_kept.txt")
    output:
        table=os.path.join(dir["output"]["reads_statistics"], "combined_read_counts.txt"),
        figure=os.path.join(dir["output"]["reads_statistics"], "reads_composition_barplot.svg")
    params:
        script=os.path.join(dir["scripts"], "combine_read_counts_and_plot.R")
    conda:
        os.path.join(dir["env"], "R.yml")
    script:
        "{params.script}"