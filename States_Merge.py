#This one is mainly for those needs merge.
#Namely:  Y-II, Na-I, K-I, Ca-II, Zn-I


#For states merge, it is quite simple.We just read in all the pre-processed .states files separately and then merge them with duplicates removed.
#The index is then re-given after the merged file is sorted by energy level

import pandas as pd
column_names = ['Index', 'E', 'g_j', 'J', 'Uncertainty', 'Life(s)', 'g_lande', 'Configuration', 'Term', 'Abbr']
element = "Zn-I"
States_to_transz = "Kurucz-" + element + "/States_Finalz.csv"
output_filenamez = "Kurucz/Kurucz" + element +"z.states"
States_to_transy = "Kurucz-" + element + "/States_Finalw.csv"
output_filenamey = "Kurucz/Kurucz" + element +"w.states"
States_to_trans = "Kurucz-" + element + "/States_Final.csv"
output_filename = "Kurucz/Kurucz" + element +".states"
# Read the file into a pandas DataFrame
df_z = pd.read_csv(output_filenamez, delim_whitespace=True, names=column_names)
df_w = pd.read_csv(output_filenamey, delim_whitespace=True, names=column_names)

combined_df = pd.concat([df_z,df_w])
combined_df = combined_df.drop_duplicates()

combined_df.drop(["Index"],axis=1,inplace=True)
combined_df.insert(0, 'Index', range(1, len(combined_df) + 1))
print(combined_df)

def format_energy(value):
    int_part = len(str(abs(int(value))))
    if int_part > 5:
        return f"{value:12.5f}"
    else:
        return f"{value:12.6f}"

if combined_df["J"].iloc[0].is_integer():
    format_str = ("{:>12d} {:>12} {:>6d} {:>7d} {:>12.6f} {:>12.4e} {:>10.6f} {:<12} {:<7} {:>2}\n")
else:
    format_str = ("{:>12d} {:>12} {:>6d} {:>7.1f} {:>12.6f} {:>12.4e} {:>10.6f} {:<12} {:<7} {:>2}\n")


output_file = output_filename
#output_file = 'Kurucz/KuruczBaII.states'
with open(output_file, 'w') as f:
    for index, row in combined_df.iterrows():
        f.write(format_str.format(
            row['Index'], format_energy(row['E']), row['g_j'], row['J'],
            row['Uncertainty'], row['Life(s)'], row['g_lande'],
            row['Configuration'], row['Term'],row["Abbr"]
        ))

print(f"Data has been written to {output_file}")


# This part also merge the csv file since the csv file is an important data when mapping those states to the trans files.
csv_df_z = pd.read_csv(States_to_transz)
csv_df_w = pd.read_csv(States_to_transy)

combined_df_csv = pd.concat([csv_df_z,csv_df_w])
combined_df_csv = combined_df_csv.drop_duplicates()

combined_df_csv.drop(["Index"],axis=1,inplace=True)
combined_df_csv.insert(0, 'Index', range(1, len(combined_df_csv) + 1))
combined_df_csv.to_csv(States_to_trans,index=False)
print(combined_df_csv)











