# import matplotlib.pyplot as plt
import exceptions
import pipeline
import pandas as pd
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def generate_histogram(motif_pipeline):
    file_prefix = motif_pipeline.output_dir + motif_pipeline.string_name()
    df_AME_results = pd.read_csv(file_prefix + "_AME_results.csv")
    total_count = df_AME_results.shape[0]
    df_AME_results["relative_rank"] = df_AME_results["rank"] / total_count
    list_rand = df_AME_results[df_AME_results["name"].str.endswith("rand")]["relative_rank"].tolist()
    list_SV = df_AME_results[~(df_AME_results["name"].str.endswith("rand"))]["relative_rank"].tolist()
    df_results_SV = {'SV': robjects.FloatVector(list_SV)}
    r_df_results_SV = robjects.DataFrame(df_results_SV)
    df_results_rand = {'rand': robjects.FloatVector(list_rand)}
    r_df_results_rand = robjects.DataFrame(df_results_rand)
    gp = ggplot2.ggplot(r_df_results_SV)
    pp = gp + ggplot2.aes_string(x='SV') + ggplot2.geom_histogram(bins=100)
    robjects.r.ggsave(filename=file_prefix+"_plot_SV.pdf", plot=pp, width=150, height=150, unit='mm')
    gp = ggplot2.ggplot(r_df_results_rand)
    pp = gp + ggplot2.aes_string(x='rand') + ggplot2.geom_histogram(bins=100)
    robjects.r.ggsave(filename=file_prefix+"_plot_rand.pdf", plot=pp, width=150, height=150, unit='mm')