import pandas as pd

# Generates metadata from the excel with sample names

## Read excels
dir_excel = "Sample Submission Form 107420_Shiva.xlsx"
dir_cell_count = "cell_numbers.xlsx"

metadata_raw = pd.read_excel(dir_excel, "DNA or RNA samples_newIDs", header=45).iloc[4:106,] # Due to added 4 rows

with pd.ExcelFile(dir_cell_count) as xls:
    # sheets to combine
    target_sheets = ['HALLO', 'Dynamo', 'Treg']
    # read and then concat
    cell_counts = [pd.read_excel(xls, s, header=0) for s in target_sheets]
    cell_counts = pd.concat(cell_counts, ignore_index=True)

# fix column names
metadata_raw = metadata_raw.iloc[:,[2,9]]
metadata_raw.columns = ["gs.id", "remarks"]
cell_counts.columns = ["individual", "organ", "subset", "cell.count"]

# split individual organ and subset in metadata_raw
metadata_raw.remarks = metadata_raw.remarks.str.strip()
metadata_raw[["individual", "organ", "subset"]] = metadata_raw.remarks.str.extract(r"(^\D*\d+)\s+(\w+)\s+(.*)")

metadata_raw.subset.unique()



# remove whitespaces in strings)in matadata
metadata_raw = metadata_raw[metadata_raw.columns].apply(lambda x: x.str.strip())
cell_counts[["organ", "individual", "subset"]] = cell_counts[["organ", "individual", "subset"]].apply(lambda x: x.str.strip())


# populate metadata with cell numbers
keys = ["individual", "organ", "subset"]
cell_counts = cell_counts.drop_duplicates(subset=keys)
merged = metadata_raw.merge(cell_counts, on = keys, how = "left", suffixes = ("", "_new"))

# change "cell type" names for compatibility
change_by = {"CD4CM": "CD4 CM",
             "CD4EM+EMRA": "CD4 EM+EMRA",
             "CD8EM+CM": "CD8 EM+CM",
             "CD8EMRA": "CD8 EMRA",
             "CD4N": "CD4 N", #LUNG & PP
             "CD8N": "CD8 N", #LUNG & PP
             'CD4CD69+': "CD4 CD69+",
             'CD4CD69-': "CD4 CD69-",
             'CD8CD69+': "CD8 CD69+",
             'CD8CD69-': "CD8 CD69-",
             'CD4Treg Naïve': "CD4 TregNaive",
             'CD4Treg Memory': "CD4 TregMemory"
             }

merged["subset"] = merged["subset"].replace(change_by)
                                                         
# fix again individual
merged.individual = merged.individual.str.split(" ").str.join("")

# group by individual, organ and cell to determine which fastq files must be merged
merged = merged.groupby(["individual", "organ", "subset"])['gs.id'].apply(list).reset_index()
breakpoint()
merged = merged.groupby(["individual", "organ", "subset"]).agg({
    'gs.id': list,           # Turn IDs into a list
    'cell.count': 'unique'   
}).reset_index()




# populate metadata with cell numbers
merged.to_csv("metadata_merged.tsv", index=False, sep="\t")
