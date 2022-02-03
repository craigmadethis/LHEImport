from LHEImport.LHEImport2 import read_lhe, tohdf5,extractparams
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import gc
dir = "./mil_lhe_files/"
# FILE = DIR+"run_1.lhe"
# FILE = "./event_files/lhe/eft/run_09.lhe"
# for file in os.listdir(DIR):
#     data = read_lhe(DIR+file)
#     tohdf5(data, "./mil_lhe_files/mil", key=f"{file.split('.')[0]}")



# data = read_lhe(dir+'run_1.lhe')

# for file in os.listdir(dir):
#     # print(dir+file)
#     data=read_lhe(dir+file)
#     extractparams(data, "big", key=f"{file.split('.')[0]}")
#     gc.collect()

# df = pd.read_hdf("big.h5", key="run_1")
# appended_data = []
# for file in os.listdir("./mil_lhe_files/"):
#     data = pd.read_hdf("big.h5", key=file.split('.')[0])
#     appended_data.append(data)
#
# df = pd.concat(appended_data)
#
# df.to_hdf("big.h5", key="all")







df = pd.read_hdf("big.h5", key="all")
print(df.head())
print(len(df))

weighted=[d.get('rwgt_1') for d in df.weights]
reweighted=np.array([d.get('rwgt_5') for d in df.weights])
#
bins=np.linspace(0,600,50)

pseudo_data = np.random.choice(df['pt_z'], size=200000, replace =True, p=reweighted/reweighted.sum())

def unweight(data, current_points):
    rndm = np.random.uniform(0,1,1000000) * max(reweighted)
    selected = np.where(rndm<reweighted)[0]
    current_points = np.append(current_points, data.iloc[selected])
    return current_points

sampled_points = np.array([])
while sampled_points.size<200000:
    sampled_points = unweight(df['pt_z'], sampled_points)

plt.hist(df['pt_z'], bins=bins, density=True, weights=reweighted, alpha=0.3)
plt.hist(pseudo_data, bins=bins, histtype='step', density=True)
plt.hist(sampled_points, bins=bins, histtype='step', density=True)
plt.show()
