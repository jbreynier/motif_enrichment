import os
import subprocess
import sys
import exceptions
import pipeline
import refgenome

def bedtools(motif_pipeline, genome):
    '''Extract fasta sequences from genome using bedtools'''
    file_prefix = motif_pipeline.subdir_name + motif_pipeline.prefix
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

def merge(motif_pipeline):
    '''Merge all the random/normal SV breakpoint files'''
    list_concat_sv = []
    list_concat_rand = []
    for file_name in motif_pipeline.list_bedpe:
        for sv_type in motif_pipeline.SV_types:
            sv_bed = motif_pipeline.output_dir + "bed_files/" + file_name + "_" + sv_type + "_sv.bed"
            rand_bed = motif_pipeline.output_dir + "bed_files/" + file_name + "_" + sv_type + "_rand" + str(motif_pipeline.rand_sv_ratio) + ".bed"
            if os.path.isfile(sv_bed):
                list_concat_sv.append(sv_bed)
            if os.path.isfile(rand_bed):
                list_concat_rand.append(rand_bed)
    with open(motif_pipeline.subdir_name + motif_pipeline.prefix + "_rand.bed", 'w') as out_rand:
        for rand_path in list_concat_rand:
            with open(rand_path) as rand_file:
                out_rand.write(rand_file.read())
    with open(motif_pipeline.subdir_name + motif_pipeline.prefix + "_sv.bed", 'w') as out_sv:
        for sv_path in list_concat_sv:
            with open(sv_path) as sv_file:
                out_sv.write(sv_file.read())

def FIMO(motif_pipeline):
    '''Runs FIMO (Find Individual Motif Occurrences) program on samples'''
    file_prefix = motif_pipeline.subdir_name + motif_pipeline.prefix
    for file_path in [file_prefix + "_rand.fasta", file_prefix + "_sv.fasta"]:
        if not os.path.isfile(file_path):
            message = ("Error: the following file <{path}> "
                    "is missing").format(path=file_path)
            raise exceptions.IncorrectPathError(message)
    for rand_sv in ["rand", "sv"]:
        FIMO_script = ("fimo --oc {output_dir} --thresh {FIMO_thresh} --max-stored-scores 100000000 {meme_file} "
                            "{fasta_file}").format(output_dir=motif_pipeline.subdir_name+"FIMO_"+rand_sv,
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
    '''Run AME (Analysis of Motif Enrichment) program on samples'''
    file_prefix = motif_pipeline.subdir_name + motif_pipeline.prefix
    for file_path in [file_prefix + "_rand.fasta", file_prefix + "_sv.fasta"]:
        if not os.path.isfile(file_path):
            message = ("Error: the following file path <{path}> "
                    "is missing").format(path=file_path)
            raise exceptions.IncorrectPathError(message)
    AME_script = ("ame --verbose 5 --oc {output_dir} "
                    "--scoring {method} "
                    "--method ranksum "
                    "--control {fasta_rand} "
                    "{fasta_sv} {meme_file} ").format(output_dir=motif_pipeline.subdir_name+"AME",
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