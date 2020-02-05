import glob
import csv
import numpy as np
import sys
import pandas as pd
import random
import refgenome
import pipeline
from scipy.stats import hypergeom

def bedpe_to_bed(genome, motif_pipeline):
    dict_sv = {"DEL": "del", "INV": "inv", "TRA": "tra", "DUP": "dup"}
    file_prefix = motif_pipeline.output_dir + motif_pipeline.string_name()
    df_include_intervals = pd.read_csv(genome.path_include, sep="\t", names=["chrom", "start", "end"], usecols=[0,1,2], header=0)
    with open(file_prefix + "_sv.bed", "w") as output_sv, open(file_prefix + "_rand.bed", "w") as output_rand:
        writer_sv = csv.writer(output_sv, delimiter='\t')
        writer_rand = csv.writer(output_rand, delimiter='\t')
        for file in glob.glob(motif_pipeline.input_dir + "*.bedpe"):
            file_name = (file.split("/")[-1]).split(".")[0]
            if file_name in motif_pipeline.list_bedpe:
                with open(file, "r") as csv_file:
                    file_name = (file.split("/")[-1]).split(".")[0]
                    reader = csv.DictReader(csv_file, delimiter='\t')
                    list_rows_sv = []
                    list_rows_rand = []
                    indexes = ["1", "2"]
                    for row in reader:
                        if dict_sv[row["svclass"][-3:]] in motif_pipeline.SV_types:
                            for index in indexes:
                                chrom = "chr" + row["chrom" + index]
                                df_include_chrom = df_include_intervals[df_include_intervals["chrom"] == chrom]
                                start = int(row["start" + index]) - 50
                                end = int(row["start" + index]) + 50
                                if df_include_chrom.apply(lambda x: x["start"] <= start and x["end"] >= end, axis=1).any():
                                    # print("Included:")
                                    # print(start, end)
                                    # print(df_include_chrom.loc[(df_include_chrom["start"] <= start) & (df_include_chrom["end"] >= end)])
                                    list_rows_sv.append([row["chrom" + index],
                                        str(start),
                                        str(end),
                                        file_name + "_" + "chr" + row["chrom" + index]
                                        + ":" + str(int(row["start" + index]) - 50) 
                                        + "-" + str(int(row["start" + index]) + 50) 
                                        + "_" + dict_sv[row["svclass"][-3:]],
                                        "100",
                                        row["strand" + index]
                                    ])
                                    for _ in range(motif_pipeline.rand_sv_ratio):
                                        rand_position = random.randint(1, genome.length_dict[chrom])
                                        while not df_include_chrom.apply(lambda x: x["start"] <= rand_position and x["end"] >= rand_position + 100, axis=1).any():
                                            rand_position = random.randint(1, genome.length_dict[chrom])
                                        list_rows_rand.append([row["chrom" + index],
                                            str(rand_position),
                                            str(rand_position + 100),
                                            file_name + "_" + "chr" + row["chrom" + index]
                                            + ":" + str(rand_position) + "-" + str(rand_position + 100)
                                            + "_" + dict_sv[row["svclass"][-3:]] + "_rand",
                                            "100",
                                            row["strand" + index]
                                        ])
                                # else:
                                    # print("Not included:")
                                    # print(df_include_chrom.loc[(df_include_chrom["start"] <= start)].iloc[[-1]])
                                    # print(start, end)
                    writer_sv.writerows(list_rows_sv)
                    writer_rand.writerows(list_rows_rand)

def extract_list_sequences_AME(motif_pipeline):
    list_data = []
    file_prefix = motif_pipeline.output_dir + motif_pipeline.string_name()
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
    df_results_sv = pd.read_csv(motif_pipeline.output_dir + "FIMO_sv/fimo.tsv", sep="\t")
    df_results_rand = pd.read_csv(motif_pipeline.output_dir + "FIMO_rand/fimo.tsv", sep="\t")
    significant_sv = df_results_sv["sequence_name"].nunique()
    significant_rand = df_results_rand["sequence_name"].nunique()
    results_hypergeom = hypergeom.sf(significant_sv - 1,
                            motif_pipeline.num_SV_breakpoints * (motif_pipeline.rand_sv_ratio + 1),
                            significant_sv + significant_rand,
                            motif_pipeline.num_SV_breakpoints)
    file_prefix = motif_pipeline.output_dir + motif_pipeline.string_name()
    with open(file_prefix + "_results_summary.txt", "a+") as summary:
        summary.write("Hyper-geometric test of FIMO results:\n")
        summary.write("p-value: " + str(results_hypergeom) + "\n")

def extract_output_AME(motif_pipeline):
    df_results_AME = pd.read_csv(motif_pipeline.output_dir + "AME/ame.tsv", sep="\t", comment="#")
    file_prefix = motif_pipeline.output_dir + motif_pipeline.string_name()
    with open(file_prefix + "_results_summary.txt", "a+") as summary:
        summary.write("One-tailed Wilcoxon rank-sum test of AME results:\n")
        summary.write("p-value: " + str(df_results_AME["p-value"].iloc[0]) + "\n")
