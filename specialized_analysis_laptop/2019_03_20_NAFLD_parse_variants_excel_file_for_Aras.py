import pandas as pd 
genes1 = ["PNPLA3", "MBOAT7", "TM6SF2","APOB", "LIPA", "LPIN1", "FATP5", "MTTP", "GCKR", "KLF6", "IL28B", "SOD2", "TRIB1", "LYPLAL1", "PPP1R3B", "HSD17B13", "UCP2", "ENPP1", "IRS1", "MERTK", "TERT"]
genes2 = ["PCSK9", "SIRT1", "IFN-Lamda3", "IFN-Lambda4", "CREBH", "DNMTs", "TNFalpha", "ApoC3"]
genes3 = ["APOB", "APOA1", "APOE"]
genes4 = ["APOL4","APOL3", "APOL5", "ABCG5", "ABCG8", "ABCA9", "ACACB", "APOL1", "APOB"]

variants = pd.read_csv("Box Sync/laptop_folders/Mattis/Mattis_NAFLD/NAFLD/NAFLD_variant_calls_filtered.txt", sep="\t")
variants_coding = variants[variants['Codons'] != "-"]
variants_coding = variants_coding[["SYMBOL", "Gene", "Location", "Allele", "Codons", "Existing_variation", "Protein_position", "Amino_acids", "SIFT", "PolyPhen"]]
print (variants.shape)
print (variants_coding.shape)
variants_genelist1 = variants_coding[variants_coding["SYMBOL"].isin(genes1)]
variants_genelist2  = variants_coding[variants_coding["SYMBOL"].isin(genes2)]
variants_genelist3  = variants_coding[variants_coding["SYMBOL"].isin(genes3)]
variants_genelist4  = variants_coding[variants_coding["SYMBOL"].isin(genes4)]

variants_genelist1.drop_duplicates(keep="first", inplace=True)
variants_genelist2.drop_duplicates(keep="first", inplace=True)
variants_genelist3.drop_duplicates(keep="first", inplace=True)
variants_genelist4.drop_duplicates(keep="first", inplace=True)

variants_genelist1.to_csv("2019_03_20_Mattis_genelist1_strong_evidence.txt", sep="\t", index=False)
variants_genelist2.to_csv("2019_03_20_Mattis_genelist2_weak_evidence.txt", sep="\t", index=False)
variants_genelist3.to_csv("2019_03_20_Mattis_genelist3_apo_genes_only.txt", sep="\t", index=False)
variants_genelist4.to_csv("2019_03_20_Mattis_genelist4_lipid_transport.txt", sep="\t", index=False)

