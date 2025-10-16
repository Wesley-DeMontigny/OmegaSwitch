import os
import pandas as pd

dirs = os.listdir("./")
sim_files = [f for f in dirs if os.path.isfile(f) and "simulation." in f]

true_labels = []
posterior_probs_1 = []
posterior_probs_2 = []
posterior_probs_3 = []
dNdS_cov = []

for f in sim_files:
    analysis_file = f.replace("simulation", "analysis")
    sim_df = pd.read_csv(f, sep="\t")
    analysis_df = pd.read_csv(analysis_file, sep="\t")

    true_omega = sim_df["OmegaCount"][0]
    true_labels.append(true_omega)
    probs = [analysis_df.loc[analysis_df["OmegaCount"] == i].shape[0] / analysis_df.shape[0] for i in range (1,4)]
    posterior_probs_1.append(probs[0])
    posterior_probs_2.append(probs[1])
    posterior_probs_3.append(probs[2])

    


out_df = pd.DataFrame({"true_labels": true_labels, "probs_1": posterior_probs_1, "probs_2": posterior_probs_2, "probs_3": posterior_probs_3})

out_df.to_csv("posterior_stats.csv")
