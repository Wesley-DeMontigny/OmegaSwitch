import pandas as pd
import numpy as np
from scipy.stats import norm
import sys
import glob

def rank_normalize(file):
    print(f"Normalizing {file}")
    df = pd.read_csv(file, sep="\t")
    for index, name in enumerate(df):
            if(name != "Iteration"):
                    print(f"Rank-Normalizing {name}")
                    vals = df[name].to_numpy()
                    vals = np.sort(vals)
                    
                    ranks = np.arange(1, vals.shape[0]+1)
                    unique_vals = np.unique(vals)
                    group_indices = {group: np.where(vals == group)[0] for group in unique_vals}
                    averaged_ranks = {group: np.mean(ranks[indices]) for group, indices in group_indices.items()}
                    
                    max_rank = max(averaged_ranks.values())
                    
                    for r in averaged_ranks:
                            averaged_ranks[r] -= 0.5
                            averaged_ranks[r] /= max_rank
                            averaged_ranks[r] = norm.ppf(averaged_ranks[r])
                    
                    df[name] = df[name].replace(averaged_ranks)

    df.to_csv("rank_normalized_" + file, index=False, sep="\t")


analysis_trace_files = glob.glob("analysis.*")
branch_trace_files = glob.glob("branches.*")

for f in analysis_trace_files:
    rank_normalize(f)
for f in branch_trace_files:
    rank_normalize(f)