#------------ SET UP THE DIRECTORIES
dir = dict()
dir["output"] = dict()

# WORKFLOW DIRs
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")
dir["db"] = os.path.join(workflow.basedir, "..", "db")

# OUTPUT DIRs
dir["output"]["base"] = RESULTS_DIR
dir["output"]["fastqc"] = os.path.join(dir["output"]["base"], "fastqc")
dir["output"]["reads_processing"] = os.path.join(dir["output"]["base"], "reads_processing")
dir["output"]["reads_statistics"] = os.path.join(dir["output"]["base"], "reads_statistics")

#------------ SET UP THE OUTPUT
# raw reads: read counts and fastqc
fastqc_input = [
    os.path.join(dir["output"]["fastqc"], "multiqc_raw_reads", "multiqc_report.html")
]

# mapping: remove adaptors, phix, human and host contmaination
mapping_input = [
    os.path.join(dir["output"]["fastqc"], "multiqc_after_trimmomatic", "multiqc_report.html"),
    os.path.join(dir["output"]["reads_statistics"], "combined_read_counts.txt"),
    os.path.join(dir["output"]["reads_statistics"], "reads_composition_barplot.svg")] 