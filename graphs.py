# import matplotlib.pyplot as plt
import exceptions
import pipeline
import pandas as pd
import matplotlib.pyplot as plt
# import rpy2.robjects.lib.ggplot2 as ggplot2
# from rpy2 import robjects
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.conversion import localconverter

def generate_histogram(motif_pipeline):
    '''Create histogram using AME ranking of sequences'''
    file_prefix = motif_pipeline.subdir_name + motif_pipeline.prefix
    df_AME_results = pd.read_csv(file_prefix + "_AME_results.csv")
    total_count = df_AME_results.shape[0]
    df_AME_results["relative_rank"] = df_AME_results["rank"] / total_count
    list_rand = df_AME_results[df_AME_results["name"].str.endswith("rand")]["relative_rank"].tolist()
    list_SV = df_AME_results[~(df_AME_results["name"].str.endswith("rand"))]["relative_rank"].tolist()
    fig, axs = plt.subplots(2, 1, constrained_layout=True, figsize=(4, 8))
    axs[0].hist(list_SV, bins=20, density=False)
    axs[1].hist(list_rand, bins=20, density=False)
    axs[0].title.set_text('Real SVs')
    axs[1].title.set_text('Random SVs')
    axs[0].set_xlabel('Relative rank')
    axs[0].set_ylabel('SV count')
    axs[1].set_xlabel('Relative rank')
    axs[1].set_ylabel('SV count')
    fig.suptitle("{prefix} AME results".format(prefix=motif_pipeline.prefix), fontsize=16)
    plt.savefig("{file_prefix}_AME_rank.pdf".format(file_prefix=file_prefix))