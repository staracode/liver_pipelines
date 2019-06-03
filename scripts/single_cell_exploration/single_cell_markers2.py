import scanpy
import pandas as pd
#https://icb-scanpy.readthedocs-hosted.com/en/stable/basic_usage.html

inputfile="/proj_single-cell-3/output_counts/20190226-SATO_fixKS17/KS_17/outs/raw_feature_bc_matrix.h5"
outputfile="/proj_single-cell-3/ks17_egfp.h5ad"

markergenesfile ="fetal_cell_markers_genes.txt"
with open(markergenesfile) as f:
    markergenes = f.read().splitlines()

#strong markers 
sm=["Ctsj","Ctsq","Ctsr","Cts6","Prl2b1","Prl3b1","Rhox9"]
#medium 
mm=["Bex1","Cited1","Krt8","Krt18","Lepr","Prl2c5","Prl7d1","Rhox6","Rhox12","Sct","Sparc"]
#weak markers
wm=["1700011M02Rik","Fnd3c2","Fthl17a","Ghrh","Gm9112","H19","Hsd17b2","Mdk","Procr","Trap1a"]
#pbmc 
pbmcm=["B2m","H2-D1","Junb","Malat1"]
#[x in ks17.var.index.tolist() for x in pbmcm]
nm=["Cd52","Coro1a"]
egfp=["EGFP_Rabbit"]

adata = scanpy.read_10x_h5(inputfile)
adata_markers = adata[:,sm+mm+wm+pbmcm+nm+egfp]
keep=adata_markers[adata_markers.X>0]

for marker in sm+mm+wm+pbmcm+nm+egfp:
for marker in sm:
    adata_markers_keep=adata_markers[adata_markers[:,marker].X>0].obs.index.tolist()

#take the union of barcode/rows for each marker
#look at the results in browser
adata_markers2 = adata_markers[adata_markers_keep,:]

#keep = keep.var_names_make_unique()
keep.write(outputfile)

scandpy.pl.highest_expr_genes(keep, n_top=20)