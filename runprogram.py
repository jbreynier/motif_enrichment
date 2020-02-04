import os
import subprocess
import sys
import exceptions
import pipeline
import refgenome

def bedtools(motif_pipeline, genome):
    file_prefix = motif_pipeline.output_dir + motif_pipeline.string_name()
    for file_path in [file_prefix + "_rand.bed", file_prefix + "_sv.bed"]:
        if not os.path.isfile(file_path):
            message = ("Error: the following file path <{path}> "
                    "is missing").format(path=file_path)
            raise exceptions.IncorrectPathError(message)
    
    for rand_sv in ["rand", "sv"]:
        bedtools_script = ("bedtools getfasta -name -fo {output_fasta} "
                            "-fi {genome_fasta} "
                            "-bed {input_bed}").format(output_fasta=file_prefix+"_"+rand_sv+".fasta",
                                                        genome_fasta=genome.path_fasta,
                                                        input_bed=file_prefix+"_"+rand_sv+".bed")
        bedtools_sub = subprocess.Popen(bedtools_script.split(),
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE,
                        universal_newlines = True)
        motif_pipeline.write_log("running bedtools for {rand_sv}".format(rand_sv=rand_sv))
        outs, errs = bedtools_sub.communicate()
        print(outs)
        print(errs)
        if bedtools_sub.returncode != 0:
            print("An error occurred, the program did not run to completion.")
            sys.exit()

def FIMO(motif_pipeline):
    file_prefix = motif_pipeline.output_dir + motif_pipeline.string_name()
    for file_path in [file_prefix + "_rand.fasta", file_prefix + "_sv.fasta"]:
        if not os.path.isfile(file_path):
            message = ("Error: the following file <{path}> "
                    "is missing").format(path=file_path)
            raise exceptions.IncorrectPathError(message)
    for rand_sv in ["rand", "sv"]:
        FIMO_script = ("fimo --oc {output_dir} --thresh {FIMO_thresh} {meme_file} "
                            "{fasta_file}").format(output_dir=motif_pipeline.output_dir+"FIMO_"+rand_sv,
                                                    meme_file=motif_pipeline.motif_path,
                                                    fasta_file=file_prefix+"_"+rand_sv+".fasta",
                                                    FIMO_thresh=motif_pipeline.FIMO_thresh)
        FIMO_sub = subprocess.Popen(FIMO_script.split(),
                        stdout = subprocess.PIPE,
                        stderr = subprocess.PIPE,
                        universal_newlines = True)
        motif_pipeline.write_log("running FIMO for {rand_sv}".format(rand_sv=rand_sv))
        outs, errs = FIMO_sub.communicate()
        print(outs)
        print(errs)
        if FIMO_sub.returncode != 0:
            print("An error occurred, the program did not run to completion.")
            sys.exit()

def AME(motif_pipeline):
    file_prefix = motif_pipeline.output_dir + motif_pipeline.string_name()
    for file_path in [file_prefix + "_rand.fasta", file_prefix + "_sv.fasta"]:
        if not os.path.isfile(file_path):
            message = ("Error: the following file path <{path}> "
                    "is missing").format(path=file_path)
            raise exceptions.IncorrectPathError(message)
    AME_script = ("ame --verbose 5 --oc {output_dir} "
                    "--scoring {method} "
                    "--control {fasta_rand} "
                    "{fasta_sv} {meme_file} ").format(output_dir=motif_pipeline.output_dir+"AME",
                                                    fasta_rand=file_prefix+"_rand.fasta",
                                                    fasta_sv=file_prefix+"_sv.fasta",
                                                    meme_file=motif_pipeline.motif_path,
                                                    method=motif_pipeline.AME_scoring)
    with open(file_prefix + "_AME_results.csv", "w") as output_file:
        AME_sub = subprocess.Popen(AME_script.split(),
                    stdout = subprocess.PIPE,
                    stderr = output_file,
                    universal_newlines = True)
    motif_pipeline.write_log("running AME")
    outs, _ = AME_sub.communicate()
    print(outs)
    if AME_sub.returncode != 0:
        print("An error occurred, the program did not run to completion.")
        sys.exit()