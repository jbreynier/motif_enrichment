import glob
import csv
import sys
import pandas as pd
import random
import refgenome
import pipeline
from scipy.stats import hypergeom
from datetime import datetime
from parsl import python_app


# BUT WHAT IF NONE FOR A CERTAIN SV?? TRY TO RUN AGAIN??
# SOLUTION: don't create files and only concatenate once everything is done running
# + add warning "no __ SVs for sample ___"

# DON'T FORGET INT

# @bash_app
# def print_hello():
#     return "echo hello"


@python_app
def bedpe_to_bed(genome, motif_pipeline, sample_name, SV_types):
    '''
    Extract SVs from bedpe file and generate random SVs.
    Write real and random SVs for each SV type to separate files.
    '''
    import random
    import csv
    import glob
    dict_sv = {"DEL": "del", "INV": "inv", "TRA": "tra", "DUP": "dup"}
    dict_output_sv = {"del": [], "inv": [], "tra": [], "dup": []}
    dict_output_rand = {"del": [], "inv": [], "tra": [], "dup": []}
    file_prefix = motif_pipeline.output_dir + "bed_files/" + sample_name
    df_include_intervals = pd.read_csv(genome.path_include, sep="\t", names=["chrom", "start", "end"], usecols=[0,1,2], header=0)
    bedpe_file = glob.glob(motif_pipeline.input_dir+sample_name+"*.bedpe")[0]
    with open(bedpe_file) as csv_file:
        reader = csv.DictReader(csv_file, delimiter='\t')
        # list_rows_sv = []
        # list_rows_rand = []
        indexes = ["1", "2"]
        for row in reader:
            sv_type = dict_sv[row["svclass"][-3:]]
            if sv_type in motif_pipeline.SV_types:
                for index in indexes:
                    chrom = "chr" + row["chrom" + index]
                    df_include_chrom = df_include_intervals[df_include_intervals["chrom"] == chrom]
                    start = int(float(row["start" + index])) - 50
                    end = int(float(row["start" + index])) + 50
                    if df_include_chrom.apply(lambda x: x["start"] <= start and x["end"] >= end, axis=1).any():
                        # print("Included:")
                        # print(start, end)
                        # print(df_include_chrom.loc[(df_include_chrom["start"] <= start) & (df_include_chrom["end"] >= end)])
                        dict_output_sv[sv_type].append([row["chrom" + index],
                            str(start),
                            str(end),
                            sample_name + "_" + "chr" + row["chrom" + index]
                            + ":" + str(int(row["start" + index]) - 50) 
                            + "-" + str(int(row["start" + index]) + 50) 
                            + "_" + sv_type,
                            "100",
                            row["strand" + index]
                        ])
                        for _ in range(motif_pipeline.rand_sv_ratio):
                            rand_position = random.randint(1, genome.length_dict[chrom])
                            while not df_include_chrom.apply(lambda x: x["start"] <= rand_position and x["end"] >= rand_position + 100, axis=1).any():
                                rand_position = random.randint(1, genome.length_dict[chrom])
                            dict_output_rand[sv_type].append([row["chrom" + index],
                                str(rand_position),
                                str(rand_position + 100),
                                sample_name + "_" + "chr" + row["chrom" + index]
                                + ":" + str(rand_position) + "-" + str(rand_position + 100)
                                + "_" + sv_type + "_rand",
                                "100",
                                row["strand" + index]
                            ])
    for sv_type in SV_types:
        sv_path = file_prefix+"_"+sv_type+"_sv.bed"
        rand_path = file_prefix+"_"+sv_type+"_rand"+str(motif_pipeline.rand_sv_ratio)+".bed"
        with open(sv_path, "w") as output_sv, open(rand_path, "w") as output_rand:
            writer_sv = csv.writer(output_sv, delimiter='\t')
            writer_rand = csv.writer(output_rand, delimiter='\t')
            writer_sv.writerows(dict_output_sv[sv_type])
            writer_rand.writerows(dict_output_rand[sv_type])

def extract_list_sequences_AME(motif_pipeline):
    '''Extract ranked list of matches from AME error log'''
    print("extracting AME csv")
    list_data = []
    file_prefix = motif_pipeline.subdir_name + motif_pipeline.prefix
    with open(file_prefix + "_AME_results.csv", "r") as results_file:
            line = results_file.readline()
            while line:
                if not line.startswith("M2: "):
                    line = results_file.readline()
                else:
                    seq_name = line.split()[4]
                    seq_score = float(line.split()[8])
                    seq_rank = int(float(line.split()[10]))
                    list_data.append([seq_name, seq_score, seq_rank])
                    line = results_file.readline()
    df_results = pd.DataFrame.from_records(list_data, columns=["name", "score", "rank"])
    df_results.to_csv(file_prefix + "_AME_results.csv", index=False)

def extract_output_FIMO(motif_pipeline):
    '''Extract results from FIMO output'''
    print("extracting FIMO")
    df_results_sv = pd.read_csv(motif_pipeline.subdir_name + "FIMO_sv/fimo.tsv", sep="\t")
    df_results_rand = pd.read_csv(motif_pipeline.subdir_name + "FIMO_rand/fimo.tsv", sep="\t")
    significant_sv = df_results_sv["sequence_name"].nunique()
    significant_rand = df_results_rand["sequence_name"].nunique()
    results_hypergeom = hypergeom.sf(significant_sv - 1,
                            motif_pipeline.num_SV_breakpoints * (motif_pipeline.rand_sv_ratio + 1),
                            significant_sv + significant_rand,
                            motif_pipeline.num_SV_breakpoints)
    file_prefix = motif_pipeline.subdir_name + motif_pipeline.prefix
    with open(file_prefix + "_results_summary.txt", "a+") as summary:
        summary.write("Hyper-geometric test of FIMO results:\n")
        summary.write("p-value: " + str(results_hypergeom) + "\n")

def extract_output_AME(motif_pipeline):
    '''Extract results from AME output'''
    print("extracting AME output")
    df_results_AME = pd.read_csv(motif_pipeline.subdir_name + "AME/ame.tsv", sep="\t", comment="#")
    file_prefix = motif_pipeline.subdir_name + motif_pipeline.prefix
    with open(file_prefix + "_results_summary.txt", "a+") as summary:
        summary.write("One-tailed Wilcoxon rank-sum test of AME results:\n")
        summary.write("p-value: " + str(df_results_AME["p-value"].iloc[0]) + "\n")
