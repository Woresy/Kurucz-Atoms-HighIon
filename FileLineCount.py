import pandas as pd

#This file is used to get the number of rows in states and trans file by giving the atom name.
element = "Ba"

output_filename = "Kurucz/"+ element+"/"+ element+"/" +"Kurucz/"+element+"__Kurucz.states"
column_names= ["Index1","Index2","A(10base)","wn(cm-1)"]
df = pd.read_csv(output_filename, delim_whitespace=True, names=column_names)

print(len(df))

output_filename = "Kurucz/"+ element+"/"+ element+"/"+"Kurucz/"+element+"__Kurucz.trans"
column_names= ["Index1","Index2","A(10base)","wn(cm-1)"]
df = pd.read_csv(output_filename, delim_whitespace=True, names=column_names)

print(len(df))
