
import pandas as pd
import re

#This python file mainly serves for the .states processing based on the .gam and lifetime files.
#Remember change atom names,(remove top part in gam and lifetime file) and change the mapping table each time.
#For multiple .gam files, remember to add suffix.
element = "N-I"
states_file = "Kurucz-" + element + "/GAM.csv"
Life_file = "Kurucz-" + element + "/LIFE.csv"
States_to_trans = "Kurucz-" + element + "/States_Final.csv"
output_filename = "Kurucz/Kurucz" + element +".states"


# Here is the mapping table, this could be found on top of each .gam file in kurucz atomic database.
mapping = """
1 s2p3      2 s2p2 3p   3 s2p2 4p   4 s2p2 5p   5 s2p2 6p   6 s2p2 7p   
7 s2p2 8p   8 s2p2 9p   9 s2p2 10p  A s2p2 11p  B s2p2 12p  C s2p2 13p  
D p5        E p4 3p     F p4 4p     G p4 5p     H p4 6p     I p4 7p     
J p4 8p     K p4 9p     L p4 10p    M p4 11p    N p4 12p    O p4 13p    
P s2p2 4f   Q s2p2 5f   R s2p2 6f   S s2p2 7f   T s2p2 8f   U s2p2 9f   
V s2p2 10f  W s2p2 11f  X s2p2 12f  Y s2p2 13f  Z p4 4f     a p4 5f     
b p4 6f     c p4 7f     d p4 8f     e p4 9f     f p4 10f    g p4 11f    
h p4 12f    i p4 13f    j s2p2 6h   k s2p2 7h   l s2p2 8h   m s2p2 9h   
n p4 6h     o p4 7h     p p4 8h     q p4 9h     r s2p2 8k   s s2p2 9k   
t p4 8k     u p4 9k     v sp3 3s    w sp3 4s    x sp3 5s    y sp3 6s    
z sp3 7s    ! sp3 8s    " sp3 9s    # sp3 10s   $ sp3 3d    % sp3 4d    
& sp3 5d    ' sp3 6d    ( sp3 7d    ) sp3 8d    * sp3 9d    + sp3 10d   
,           -           .           /           :           ;           
<           =           >           ?           @           [           
\           ]           ^           _           `           {           
|           }           ~           0           
1 s2p2 3s   2 s2p2 4s   3 s2p2 5s   4 s2p2 6s   5 s2p2 7s   6 s2p2 8s   
7 s2p2 9s   8 s2p2 10s  9 s2p2 11s  A s2p2 12s  B s2p2 13s  C p4 3s     
D p4 4s     E p4 5s     F p4 6s     G p4 7s     H p4 8s     I p4 9s     
J p4 10s    K p4 11s    L p4 12s    M p4 13s    N s2p2 3d   O s2p2 4d   
P s2p2 5d   Q s2p2 6d   R s2p2 7d   S s2p2 8d   T s2p2 9d   U s2p2 10d  
V s2p2 11d  W s2p2 12d  X p4 3d     Y p4 4d     Z p4 5d     a p4 6d     
b p4 7d     c p4 8d     d p4 9d     e p4 10d    f p4 11d    g p4 12d    
h s2p2 5g   i s2p2 6g   j s2p2 7g   k s2p2 8g   l s2p2 9g   m p4 5g     
n p4 6g     o p4 7g     p p4 8g     q p4 9g     r s2p2 7i   s s2p2 8i   
t s2p2 9i   u p4 7i     v p4 8i     w p4 9i     x sp4       y sp3 3p    
z sp3 4p    ! sp3 5p    " sp3 6p    # sp3 7p    $ sp3 8p    % sp3 9p    
& sp3 10p   ' sp3 11p   ( sp3 12p   ) sp3 13p   * sp3 14p   + sp3 4f    
, sp3 5f    - sp3 6f    . sp3 7f    / sp3 8f    : sp3 9f    ;           
<           =           >           ?           @           [           
\           ]           ^           _           `           {           
|           }           ~           0           
"""

# The mapping table has two version, one is for odd and one is for even.
# If the .gam starts with an odd state, then the mapping table will also begin with odd. Vice versa.
# Use the first mapping in the bottom's mapping time as a split word.

##start with odd
eve_1 = "1 s2p2 3s"
odd_part, eve_part = mapping.split(eve_1)
eve_part = eve_1 + eve_part

##start with eve
# odd_1 = "1 s2p4 3s"
# eve_part, odd_part = mapping.split(odd_1)
# odd_part = odd_1 + odd_part




# Here we read the .gam files and then remove duplicates.
column_name = ["ELEM","Index","E","J","label","g_lande"]
states = pd.read_csv(states_file,names= column_name)

states.drop(["Index"],axis=1,inplace=True)
states.drop_duplicates(subset=["ELEM", "E", "J", "label"], inplace=True)
states.reset_index(drop=True, inplace=True)
index = range(1,len(states)+1)
index = pd.DataFrame(index,columns=["Index"])
states = pd.concat([index,states],axis = 1)

# g_j is computed and uncertainty is assigned with 0.1 as default
states["g_j"] = 2 * states["J"] + 1
states["g_j"] = states["g_j"].astype(int)
states["Uncertainty"] = 0.1

# Here we read in the lifetime file and remove duplicates. Lifetime is turned from ns to s to fit exomol format.
column_name = ["ELEM","Index","E","J","label","SUM_A","Life1","Life(ns)"]
life = pd.read_csv(Life_file,names= column_name)

life.drop(["SUM_A","Life1","ELEM","Index"],axis=1,inplace=True)

life.drop_duplicates(subset=["E", "J", "label"], inplace=True)
life.reset_index(drop=True, inplace=True)
life["Life(s)"] = life["Life(ns)"] / 1e9
life.drop(["Life(ns)"],axis=1,inplace=True)


#Lifetime is merged with the states based on E, J, label.
combined_df = states.merge(life[['E', 'J', 'label', 'Life(s)']],
                           on=['E', 'J', 'label'],
                           how = 'left')

#For some atoms that do not have lifetime provided, use nan for lifetime column.
# states["Life(s)"] = float("nan")
# combined_df = states.copy()

order = ['Index','ELEM', 'E', 'g_j', 'J', 'Uncertainty', 'Life(s)', 'g_lande', 'label']
combined_df = combined_df[order]


#Abbr is given to distinguish if the states is caculated or observed. If observed, use NI. Else, use CA.
combined_df['Abbr'] = combined_df['E'].apply(lambda x: 'CA' if x < 0 else 'NI')
combined_df["E"] = combined_df["E"].abs()

def match_table(part):
    lines = part.splitlines()
    mapping = {}
    for line in lines:
        if line.strip():
            entries = line.split()
            i = 0
            while i < len(entries):
                if len(entries[i]) == 1:
                    key = entries[i]
                    value = []
                    i += 1
                    while i < len(entries) and len(entries[i]) != 1:
                        value.append(entries[i])
                        i += 1
                    mapping[key] = ' '.join(value)
    return mapping


mapping_even = match_table(eve_part)
mapping_odd = match_table(odd_part)


# Here, we start to split the labels into Configuration and Terms.
# The detail spliting cases are given in the report and github page.
configuration_list = []
term_list = []
for index, row in combined_df.iterrows():
    label = row["label"]
    elem = row["ELEM"][-3:]
    parts = label.split()
    if len(parts) == 2:
        config, possible_term = parts
        if config.endswith("nd"):
            configuration = "unknown"
            term = "unknown"
        elif possible_term.isdigit():
            configuration = config[:-2]
            term = config[-2:]
            if elem in ['EVE', 'ERz', 'EPo']:
                configuration = mapping_even.get(configuration[0], configuration[0]) + configuration[1:]
            elif elem in ['ODD', 'ORz', 'OPo']:
                configuration = mapping_odd.get(configuration[0], configuration[0]) + configuration[1:]
        else:
            configuration = config
            term = possible_term

    elif len(parts) == 3:
        config_1, config_2,pos_term = parts
        configuration = config_1 + config_2
        term = pos_term

    else:
        if len(label)>=4 and label[-4].isupper():
            configuration = label[:-5]
            term = label[-5:-3]
            if elem in ['EVE', 'ERz', 'EPo']:
                configuration = mapping_even.get(configuration[0], configuration[0]) + configuration[1:]
            elif elem in ['ODD', 'ORz', 'OPo']:
                configuration = mapping_odd.get(configuration[0], configuration[0]) + configuration[1:]
        elif '?' in label[-4:-1]:
            configuration,term = label.split('?')[0], label.split('?')[1]
        else:
            configuration = "unknown"
            term = "unknown"

    if len(term) == 3:
        if re.match(r'^[a-zA-Z]\d[a-zA-Z]$', term):
            term = f"{term[0]}({term[1:]})"
    configuration = configuration.replace(' ','')

    configuration_list.append(configuration)
    term_list.append(term)
combined_df['Configuration'] = configuration_list
combined_df['Term'] = term_list

#Here, we remove those rows with unknown labels. States are sorted based on Energy level. The first lifetime will be infinity always.
combined_df = combined_df[(combined_df['Configuration']!='unknown') & (combined_df["Term"] != 'unknown')]
combined_df = combined_df.sort_values(by = "E", ascending= True)
combined_df.at[0, 'Life(s)'] = float('inf')

combined_df.drop(["Index"],axis=1,inplace=True)
combined_df.insert(0, 'Index', range(1, len(combined_df) + 1))
#combined_df = pd.concat([index,combined_df],axis = 1)

# The reason I save .csv for states here is that the original label is much eaiser for mapping in trans files.
combined_df_alter = combined_df.copy()
combined_df_alter.drop(["ELEM","Configuration","Term"],axis=1,inplace=True)
combined_df_alter.to_csv(States_to_trans,index=False)
#combined_df_alter.to_csv("Kurucz-Si-I/States_Final.csv",index=False)
print("Saved Done")


combined_df = combined_df.drop(columns=['label','ELEM'])

order = ['Index', 'E', 'g_j', 'J', 'Uncertainty', 'Life(s)', 'g_lande', 'Configuration','Term','Abbr']
combined_df = combined_df[order]

#The following code reformat the states data into the Exomol format.

def format_energy(value):
    # E occupies a 12-character field, but "{:12.Nf}" treats 12 as a minimum, so a
    # value with seven integer digits (E > 1e6 cm^-1, reached only in the higher
    # ionization stages) runs long and shifts every later column right.
    int_part = len(str(abs(int(value))))
    decimals = 6 if int_part <= 5 else max(0, 11 - int_part)
    text = f"{value:12.{decimals}f}"
    # Rounding can carry into a new integer digit (999999.999999 -> 1000000.0).
    while len(text) > 12 and decimals > 0:
        decimals -= 1
        text = f"{value:12.{decimals}f}"
    return text

if combined_df["J"].iloc[0].is_integer():
    format_str = ("{:>12d} {:>12} {:>6d} {:>7d} {:>12.6f} {:>12.4e} {:>10.6f} {:<12} {:<7} {:>2}\n")
else:
    format_str = ("{:>12d} {:>12} {:>6d} {:>7.1f} {:>12.6f} {:>12.4e} {:>10.6f} {:<12} {:<7} {:>2}\n")
# Write the DataFrame to a text file with the specified format
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






















