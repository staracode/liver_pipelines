import sys
import pandas as pd

all = pd.read_csv("/Users/tfriedrich/Box Sync/laptop_folders/Mattis/Mattis_NAFLD/NAFLD/NAFLD_variant_calls_filtered.txt", sep="\t")
p7017 = pd.read_csv("/Users/tfriedrich/Box Sync/2019_03_Aras_Holger/7017_ATCACG_L005_R1_001_raw_snps_indels_genotype_filtered.VEP.txt", sep="\t")
p7191 = pd.read_csv("/Users/tfriedrich/Box Sync/2019_03_Aras_Holger/7191_CGATGT_L005_R1_001_raw_snps_indels_genotype_filtered.VEP.txt", sep="\t")
p7192 = pd.read_csv("/Users/tfriedrich/Box Sync/2019_03_Aras_Holger/7192_TTAGGC_L005_R1_001_raw_snps_indels_genotype_filtered.VEP.txt", sep="\t")

p7017 = p7017[["Location", "SYMBOL", "Gene", "Feature", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "SIFT", "PolyPhen"]]
p7191 = p7191[["Location", "SYMBOL", "Gene", "Feature", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "SIFT", "PolyPhen"]]
p7192 = p7192[["Location", "SYMBOL", "Gene", "Feature", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "SIFT", "PolyPhen"]]
p7017.columns = [ x + "_7017" for x in p7017.columns.values.tolist() ] 
p7191.columns = [ x + "_7191" for x in p7191.columns.values.tolist() ] 
p7192.columns = [ x + "_7192" for x in p7192.columns.values.tolist() ] 

set7017 = set(p7017["Existing_variation_7017"])
set7191 = set(p7191["Existing_variation_7191"])
set7192 = set(p7192["Existing_variation_7192"])

#print (len(set7017.intersection(set7191)))
#print (len(set7017.difference(set7191)))

p7017p7191 = set7017.intersection(set7191)
p7017p7192 = set7017.intersection(set7192)
p7191p7192 = set7191.intersection(set7192)

overlap_all  = p7017p7191.intersection(set7192)
print ("num of variants in 7017 + 7191 + 7192: " + str(len(overlap_all)))
unique7017 = set7017.difference(overlap_all)
unique7191 = set7191.difference(overlap_all)
unique7192 = set7192.difference(overlap_all)

print ("num of variants in 7017 - ( 7017 + 7191 + 7192 ): " + str(len(unique7017)))
print ("num of variants in 7191 - ( 7017 + 7191 + 7192 ): " + str(len(unique7191)))
print ("num of variants in 7192 - ( 7017 + 7191 + 7192 ): " + str(len(unique7192)))

print ("num of variants in ( 7017 + 7191 ) - 7192: " + str(len(p7017p7191.difference(set7192))))
print ("num of variants in ( 7017 + 7192 ) - 7191: " + str(len(p7017p7192.difference(set7191))))
print ("num of variants in ( 7191 + 7192 ) - 7017: " + str(len(p7191p7192.difference(set7017))))

verify = all[all["Existing_variation"].isin(overlap_all)]
failed = all[~all["Existing_variation"].isin(overlap_all)]
print (verify.shape)
print (failed.shape)
verify.to_csv("/Users/tfriedrich/Downloads/NAFLD_verified_in_all_three_patients.xls", sep="\t", index=False)
failed.to_csv("/Users/tfriedrich/Downloads/NAFLD_not_in_all_three_patients.xls", sep="\t", index=False)

p7017[p7017["Existing_variation_7017"].isin(unique7017)].to_csv("/Users/tfriedrich/Downloads/unique_to_7017.xls", sep="\t", index = False)
p7191[p7191["Existing_variation_7191"].isin(unique7191)].to_csv("/Users/tfriedrich/Downloads/unique_to_7191.xls", sep="\t", index = False)
p7192[p7192["Existing_variation_7192"].isin(unique7192)].to_csv("/Users/tfriedrich/Downloads/unique_to_7192.xls", sep="\t", index = False)


#p7017andp7191 = p7017.merge(p7191, left_on="Location_7017", right_on="Location_7191", how="left")
#p7017andp7191.drop_duplicates(keep="first", inplace=True)
#p7017andp7191.to_csv("Downloads/p7017andp7191.xls", sep="\t", index=False)

#p7017andp7192 = p7017.merge(p7192, left_on="Location_7017", right_on="Location_7192", how="left")
#p7017andp7192.drop_duplicates(keep="first", inplace=True)
#p7017andp7192.to_csv("Downloads/p7017andp7192.xls", sep="\t", index=False)

#p7191andp7192 = p7191.merge(p7192, left_on="Location_7191", right_on="Location_7192", how="left")
#p7191andp7192.drop_duplicates(keep="first", inplace=True)
#p7191andp7192.to_csv("Downloads/p7191andp7192.xls", sep="\t", index=False)

#p7017andp7191andp7192 = p7017andp7191.merge(p7191andp7192, left_on="Location_7191", right_on="Location_7191", how="left")
#p7191andp7192.drop_duplicates(keep="first", inplace=True)
#p7017andp7191andp7192.to_csv("Downloads/p7017andp7191andp7192.xls", sep="\t", index=False)


