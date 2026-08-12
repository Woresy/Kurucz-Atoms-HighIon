import pandas as pd

#This file is to remove those rows with either nan in wavenumber/A coefficient or -1 in indexes.
#This simply prevents the invalid data inside the .trans files

#Ti skip-row 1
element = "Zr"
output_filename = "Kurucz/" + element +"_Kurucz.trans"
column_names= ["Index1","Index2","A(10base)","wn(cm-1)"]
df = pd.read_csv(output_filename, delim_whitespace=True, names=column_names)

df_filtered = df[(df['Index1'] != -1) & (df['Index2'] != -1)]

df_filtered = df_filtered.dropna(subset=['A(10base)', 'wn(cm-1)'])

print(df_filtered)

format_str = "{:>12d}{:>1}{:>12d}{:>1}{:>10.4e}{:>1}{:>15.6e}\n"

output_path = "Notes/"+output_filename
with open(output_path, 'w') as f:
    for index, row in df_filtered.iterrows():
        f.write(format_str.format(
            int(row['Index1']),'', int(row['Index2']),'',
            row['A(10base)'],'', row['wn(cm-1)']
        ))

print(f"Data has been written to {output_path}")








