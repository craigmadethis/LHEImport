import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_hdf("./event_files/hdf5/071221_eft.h5", "run_09")

weight_array=[]
for ind, row in data.iterrows():
    if row["weights"]["ctG_1"]:
        weight_array.append(row["weights"]["ctG_1"])

print(np.array(weight_array))
weight_array=np.array(weight_array)
mean = np.mean(weight_array)
stdev = np.std(weight_array)


filtered_array = []
for weight in weight_array: 
    if weight < (mean+stdev) and weight > (mean-stdev):
        filtered_array.append(weight)

xsec_o = np.sum(filtered_array)/1000

max = np.max(filtered_array)
# print(max)
# print(xsec_o)
# print(max * np.sum(filtered_array/max))
# print(filtered_array/max)
# print(np.sum(filtered_array/max)/1000)
print(np.array(filtered_array/max))
print(max* np.sum(np.array(filtered_array/max)))

plt.hist(filtered_array, weights=filtered_array/max)
# plt.scatter(filtered_array, filtered_array/max)
plt.show()

