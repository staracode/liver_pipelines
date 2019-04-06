import pandas as pd

metadata = pd.read_csv("Box Sync/laptop_folders/Mattis/endoderm_project/metadata_Aras_Mattis_2018_08_01v3.txt", sep="\t")
print (metadata)

p7017 = metadata[metadata['Name'].str.contains("7017")]
print (p7017[["Location"]])
p7192 = metadata[metadata['Name'].str.contains("7192")]
p = pd.concat([p7017, p7192])
#out="\n".join([ "scp "  + x.replace("/Volumes/My Passport/", "Box\ Sync/laptop_folders/Mattis/data/") + " tfriedrich@wyndt1.compbio.ucsf.edu:mattis_dataset"  for x in  p["Location"].values.tolist()] )
out="\n".join([ "cp "  + x.replace("/Volumes/My Passport/", "Box\ Sync/laptop_folders/Mattis/data/") + " Downloads/mattis_dataset"  for x in  p["Location"].values.tolist()] )
print (out)


iHep = metadata[metadata['cell type'].str.contains("iHep")]
out2="\n".join([ "cp "  + x.replace("/Volumes/My Passport/", "Box\ Sync/laptop_folders/Mattis/data/") + " Downloads/mattis_dataset_iHep"  for x in  iHep["Location"].values.tolist()] )
print (out2)


#scp -r Downloads/mattis_dataset/ tfriedrich@wyndt1.compbio.ucsf.edu: 
#scp -r Downloads/mattis_dataset_iHep/ tfriedrich@wyndt1.compbio.ucsf.edu: 

