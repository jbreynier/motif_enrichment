import glob
import pandas as pd

lines_out = []

for file in glob.glob("/Users/jbreynier/Desktop/Research/Yang Lab/Fall Quarter/RSS_analysis/5_19_results/*.txt"):
    if not file.split("/")[-1].startswith("COAD_READ"):
        subtype = (file.split("/")[-1]).split("_")[0]
        motif = (file.split("/")[-1]).split("_")[1]
        SV = (file.split("/")[-1]).split("_")[2]
    else:
        subtype = "COAD_READ"
        motif = (file.split("/")[-1]).split("_")[2]
        SV = (file.split("/")[-1]).split("_")[3]
    with open(file) as results:
        lines = results.readlines()
        try:
            geometric = lines[9].split("p-value: ")[-1].replace("\n", "")
            wilcoxon = lines[11].split("p-value: ")[-1].replace("\n", "")
            print(lines[9].split("p-value: ")[-1], lines[11].split("p-value: ")[-1])
            print([subtype, motif, SV, geometric, wilcoxon])
            lines_out.append([subtype, motif, SV, geometric, wilcoxon])
        except:
            print(file)
    #     geometric = lines[9].split(" ")[1]
    #     wilcoxon = lines[11].split(" ")[1]
    # lines_out.append([subtype, motif, SV, geometric, wilcoxon])

df_out = pd.DataFrame(lines_out, columns=["subtype", "motif", "SV", "FIMO_geometric", "AME_wilcoxon"])
df_out.to_csv("list_results_RSS.csv", index=False)