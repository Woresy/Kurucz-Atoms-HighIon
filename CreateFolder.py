import os

#Create Folder in correct format automatically
elements = ["Al","Al_p","Ar","Ar_p","B","B_p","Ba","Ba_p","Be","Be_p","C","C_p","Ca","Ca_p","Cl","Cl_p","Co","Co_p","Cr","Cr_p","Cu","Cu_p","F","F_p","Fe","Fe_p","K","K_p","Mg","Mg_p","Mn","Mn_p",
           "Mo","Mo_p","N","N_p","Na","Na_p","Nb","Nb_p","Ne","Ne_p","Ni","Ni_p","O","O_p","P","P_p","Pd","Pd_p","Rh","Rh_p","Ru","Ru_p",
           "S","S_p","Sc","Sc_p","Si","Si_p","Sr","Sr_p","Tc","Tc_p","Ti","Ti_p","V","V_p","Y","Y_p","Zn","Zn_p","Zr","Zr_p"]

for element in elements:
    path ="Kurucz/"+element+"/"+element+"/Kurucz"
    os.makedirs(path)


