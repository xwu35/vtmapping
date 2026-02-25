# VTmapping

A workflow for plotting the proportion of mapped and unmapped reads to the host bacteria, providing a quick initial examination of the viral tagging data.

## Set up environment

### Clone the repository 

```bash
git clone https://github.com/xwu35/vtmapping.git
```

### Install Snakemake

vtmapping is built for Snakemake version 7. Version 8 and above introduce breaking changes and deprecations and have not been tested. It may not function correctly with newer versions. Please install Snakemake version 7 using the script below.

```bash
# option 1 using conda
conda env create -n snakemake -f vtmapping/snakemake_env.yml

# option 2 using mamba if it's installed
mamba env create -n snakemake -f vtmapping/snakemake_env.yml
```

### Download snakemake profile

The profile is required to run the workflow on HPC. Skip this step if you already have a SLURM profile in `~/.config/snakemake`.

```bash
# download the profile
git clone https://github.com/xwu35/slurm

# move the profile to the right directory
mv slurm ~/.config/snakemake 
```

## Sample information table

The sample information table should look like this:

| sample   | R1                                                                                         | R2                                                                                         |
|----------|--------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| Fp22_10A | Baldridge_10A_Fp22_CD_SIC_index_0496_SIC_index_0550_TCACTACAG_CGCTACTA_S74_R1_001.fastq.gz | Baldridge_10A_Fp22_CD_SIC_index_0496_SIC_index_0550_TCACTACAG_CGCTACTA_S74_R2_001.fastq.gz | 
| Fp22_3C  | Baldridge_3C_Fp22_HHC_SIC_index_0498_SIC_index_0543_TTGAGAGTG_CAAACCTG_S20_R1_001.fastq.gz | Baldridge_3C_Fp22_HHC_SIC_index_0498_SIC_index_0543_TTGAGAGTG_CAAACCTG_S20_R2_001.fastq.gz | 
| Fp22_7B  | Baldridge_7B_Fp22_CD_SIC_index_0497_SIC_index_0547_GAGATGAAA_GTGTAGAG_S51_R1_001.fastq.gz  | Baldridge_7B_Fp22_CD_SIC_index_0497_SIC_index_0547_GAGATGAAA_GTGTAGAG_S51_R2_001.fastq.gz  | 
| Fp22_8G  | Baldridge_8G_Fp22_CD_SIC_index_0539_SIC_index_0548_TGTGAGGCT_AGGTCGCA_S64_R1_001.fastq.gz  | Baldridge_8G_Fp22_CD_SIC_index_0539_SIC_index_0548_TGTGAGGCT_AGGTCGCA_S64_R2_001.fastq.gz  | 

### Export PATH

Add the vtmapping path to your environment variable so you can run `vtreads.py` without the full path.

```bash
echo 'export PATH="/your/path/to/vtmapping:$PATH"' >> ~/.bashrc
```

To make the changes take effect, you can either log out and log back in, or you can source the `.bashrc` file using the following command:

```bash
source ~/.bashrc
```

## Usage

vtmapping supports two mapping software options (bowtie2 and minimap2), with minimap2 used by default. Detailed usage information can be viewed using the -h or --help flags `python vtmapping.py -h`.

Do not run this on the login node. Submit it as an sbatch job on the HPC using `sbatch run_vtmapping.sh`. Make sure to update the --mail-user field before submitting the job. A dry-run can be performed to check which rules will be executed and which files will be produced. 

If you did not export the path as shown above, you will need to specify the full path to vtmapping.py

```bash
conda activate snakemake

# test run
vtmapping.py --test

# general usage
vtmapping.py 
    --reads_dir /path/to/your/raw/reads \
    --sample_info /path/to/your/sample/info/table \
    --output_dir /name/your/output \
    --reference_genome /path/to/your/genome/fasta 
```

Options
```
$ vtmapping.py -h

Options:
  --reads_dir PATH         Reads directory
  --sample_info FILE       Sample information table (tab separated). The table must contain three columns (sample, R1, R2)
  --output_dir PATH        Output directory  [default: OUTPUT]
  --reference_genome FILE  Reference genome sequences
  --identity INTEGER       Percent of identity to remove reads mapped to the host genome  [default: 97]
  --adapter TEXT           Adapter sequences. By default, the adapters.fna file within `vtmapping/db` is used  [default: ""]
  --mapper TEXT            Mapping software; available options are: minimap2, bowtie2  [default: minimap2]
  --step TEXT              Steps to run; available options are: fastqc, preprocess  [default: preprocess]
  --test                   test run
  --dryrun                 Check rules to run and files to produce
  --conda_envs TEXT        Directory to store conda environments. By default, the "conda_env" directory within vtmapping is used  [default: ""]
  --profile TEXT           Snakemake profile for cluster execution  [default: slurm]
  -h, --help               Show this message and exit.
```

### Specific steps

Specific steps can be run using the `--step` flag. 

- **fastqc**: QC on raw reads. This is useful for checking your sequence quality and adapter content.
- **mapping**: remove adapters, low quality reads, phiX and human contamination reads and map remaining sequecnes to the host genome (QC on raw reads is not included in this step)

vtmapping runs mapping by default.

## Output description

- **Quality control results**: output_dir/reads_processing/fastqc
- **Read counts**: output_dir/reads_statistics/{combined_read_counts.txt,reads_composition_barplot.svg}