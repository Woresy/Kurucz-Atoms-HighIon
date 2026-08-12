import pandas as pd

#Remember change atoms names.
#Due to the huge size of some agafgf files. Some headers cannot be removed. Therefore, atoms after Ti shall add skiprows = 1 when reading agafgf file

#Ti skip-row 1
element = "N-I"
trans_file = "Kurucz-" + element + "/Line.csv"
AGA_file = "Kurucz-" + element + "/AGAFGF.csv"
States_to_trans = "Kurucz-" + element + "/States_Final.csv"
output_filename = "Kurucz/Kurucz" + element +".trans"



column_name = ["wl","log_gf","ele","E1","J1","label1","E2","J2","label2"]
Trans = pd.read_csv(trans_file,names=column_name)

#Some of the trans files' J2 has problem of adding .1 automatically to the first J2. The commented line turn that to correct value.
# Trans.loc[0,"J2"]= 2.5
# Trans["J2"] = Trans["J2"].astype(float)

print(Trans.head())

#The following code starts to turn the E, j and lable into the corresponding states index.

States = pd.read_csv(States_to_trans)
#print(States.head())
Trans["E1"] =Trans["E1"].abs()
Trans["E2"] =Trans["E2"].abs()
Trans = Trans.merge(States[["Index","E","J","label"]],left_on=["E1","J1","label1"],
                        right_on=["E","J","label"],how="left").rename(columns = {"Index":"Index1"})

Trans = Trans.merge(States[["Index","E","J","label"]],left_on=["E2","J2","label2"],
                        right_on=["E","J","label"],how="left").rename(columns = {"Index":"Index2"})


Trans = Trans.drop(columns=['E_x', 'J_x', 'label_x', 'E_y', 'J_y', 'label_y'])


cols = list(Trans.columns)
cols.insert(cols.index('E1'), cols.pop(cols.index('Index1')))
cols.insert(cols.index('E2'), cols.pop(cols.index('Index2')))
Trans = Trans[cols]


#Then we deal with the agafgf files to attain the A coefficient. We turn Log A to A with base 10. Wavenumber is also added to the trans.

column_name = ["wl","wn(cm-1)","log_gf","log_f","log_fe","log_A"]
#column_name = ["wl","wn(cm-1)","log_gf","log_f","log_fe","log_A","nouse","ele","E1","J1","label1","E2","J2","label2"]

#A_WN = pd.read_csv(AGA_file,names=column_name,skiprows=1)
A_WN = pd.read_csv(AGA_file,names=column_name)

A_WN.drop(["wl","log_gf","log_f","log_fe"],axis=1,inplace=True)
#A_WN.drop(["wl","log_gf","log_f","log_fe","ele","nouse"],axis=1,inplace=True)

A_WN["A(10base)"] = 10 ** A_WN["log_A"]
A_WN.drop(["log_A"],axis=1,inplace=True)
# print(Trans)
# print(A_WN)

#combined_df= Trans.merge(A_WN, on=["E1","J1","label1","E2","J2","label2"], how='left')
combined_df= pd.concat([A_WN,Trans],axis=1)

combined_df.drop(["wl","log_gf","ele","E1","J1","label1","E2","J2","label2"]
                 ,axis=1,inplace=True)

#Trans files are all sorted based on the wavenumber. Then they will be saved in the Exomol Format.

order = ["Index1","Index2","A(10base)","wn(cm-1)"]
combined_df = combined_df[order]
combined_df = combined_df.sort_values(by="wn(cm-1)",ascending=True)

#-1 represent not found in level
combined_df["Index2"] = combined_df["Index2"].fillna(-1).astype(int)
combined_df["Index1"] = combined_df["Index1"].fillna(-1).astype(int)


format_str = "{:>12d}{:>1}{:>12d}{:>1}{:>10.4e}{:>1}{:>15.6e}\n"

output_file = output_filename
with open(output_file, 'w') as f:
    for index, row in combined_df.iterrows():
        f.write(format_str.format(
            int(row['Index1']),'', int(row['Index2']),'',
            row['A(10base)'],'', row['wn(cm-1)']
        ))

print(f"Data has been written to {output_file}")




