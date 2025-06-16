import pandas as pd
import xarray as xr
import numpy as np
import arviz as az
import glob

def convert_to_arviz(df):
    chains = df["Chain"].unique()
    param_df = df.drop(columns=["Iteration", "Chain", "Posterior", "Prior", "Likelihood"], errors="ignore")

    trace_dict = {
        param: np.array([
            df[df["Chain"] == chain][param].to_numpy() for chain in chains
        ])
        for param in param_df.columns
    }

    return az.from_dict(posterior=trace_dict)

bayes_ess_dict = {}
bayes_files = glob.glob("rank_normalized*_Bayes")
for i, filename in enumerate(bayes_files):
    print(f"Computing ESS for {filename}")
    t_df = pd.read_csv(filename, sep="\t")
    t_df["Chain"] = 0
    idata = convert_to_arviz(t_df)
    ess_bulk = az.ess(idata, method="bulk")
    bayes_ess_dict[filename.replace("_Bayes", "")] = ess_bulk

classic_ess_dict = {}
classic_files = glob.glob("rank_normalized*_Classic")
for i, filename in enumerate(classic_files):
    print(f"Computing ESS for {filename}")
    t_df = pd.read_csv(filename, sep="\t")
    t_df["Chain"] = 0
    idata = convert_to_arviz(t_df)
    ess_bulk = az.ess(idata, method="bulk")
    classic_ess_dict[filename.replace("_Classic", "")] = ess_bulk

ess_sum = 0.0
ess_count = 0
for i in classic_ess_dict.keys():
    ratio = bayes_ess_dict[i]/classic_ess_dict[i]
    ess_sum += sum(ratio[var] for var in ratio.data_vars)
    ess_count += len(bayes_ess_dict[i].data_vars)
    
print(f"Mean Bayes to Classic Ratio: {(ess_sum/ess_count).to_numpy()}")