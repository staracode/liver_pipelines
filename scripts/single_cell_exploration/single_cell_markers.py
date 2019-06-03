python
import scanpy
#https://icb-scanpy.readthedocs-hosted.com/en/stable/basic_usage.html

ks17= scanpy.read_10x_h5("/proj_single-cell-3/output_counts/20190226-SATO_fixKS17/KS_17/outs/raw_feature_bc_matrix.h5")
ks17_egfp = ks17[:,["EGFP_Rabbit"]]

ks17[:,["EGFP_Rabbit"]].obs
ks17[:,["EGFP_Rabbit"]].X
ks17[:,["EGFP_Rabbit"]].var


#56 barcodes with eGFP expression greater than 0
len(ks17[ks17[:,["EGFP_Rabbit"]].X>0].obs)

import matplotlib.pyplot as plt
fig = plt.hist(ks17[:,["EGFP_Rabbit"]].X)
plt.title('Mean')
plt.xlabel("value")
plt.ylabel("Frequency")
plt.savefig("/proj_single-cell-3/KS17_egfp.png")

fig = plt.hist(ks17_egfp[ks17_egfp.X>0].X)
plt.title('Mean')
plt.xlabel("value")
plt.ylabel("Frequency")
plt.savefig("/proj_single-cell-3/KS17_egfp_cells_expressing_egfp.png")

#save cells that have 1 or more count for eGFP
keep.write("/proj_single-cell-3/ks17_egfp.h5ad")