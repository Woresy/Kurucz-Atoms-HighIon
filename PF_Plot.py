import pandas as pd
import matplotlib.pyplot as plt


#This python file is for ploting the partition function value


#The following code plots both the scraped Kurucz atomic partition function values and the pyexocross calculated pf values.
#This will give a direct match between the scraped one and calculated data, which proves the correctness of the data.

scrape = "Kurucz/Al/Al/Kurucz/Al__Kurucz.pf"
gen = "PyExoCross/TianyangXie/Al/27Al/Kurucz/27Al__Kurucz.pf"
column_names = ['T', 'Value']
pf_scrape = pd.read_csv(scrape, delim_whitespace=True,names=column_names)
pf_gen = pd.read_csv(gen, delim_whitespace=True,names=column_names)
filtered_data_scrape = pf_scrape[pf_scrape["T"] <= 6000]
filtered_data_gen = pf_gen[pf_gen["T"] <= 6000]

fig, axs = plt.subplots(2, 2, figsize=(10, 8))
axs[0, 0].plot(filtered_data_gen["T"], filtered_data_gen["Value"], marker='o',label="Pyexocross Generated PF")
axs[0, 0].plot(filtered_data_scrape["T"], filtered_data_scrape["Value"], marker='o',label="Kurucz Scraped PF")
axs[0, 0].set_title('Al I')
axs[0, 0].legend()
axs[0, 0].set_xlabel('Temperature')
axs[0, 0].set_ylabel('Partition Function')

scrape = "Kurucz/Cu_p/Cu_p/Kurucz/Cu_p__Kurucz.pf"
gen = "PyExoCross/TianyangXie/Cu/63Cu/Kurucz/63Cu__Kurucz.pf"
column_names = ['T', 'Value']
pf_scrape = pd.read_csv(scrape, delim_whitespace=True,names=column_names)
pf_gen = pd.read_csv(gen, delim_whitespace=True,names=column_names)
filtered_data_scrape = pf_scrape[pf_scrape["T"] <= 6000]
filtered_data_gen = pf_gen[pf_gen["T"] <= 6000]

axs[0, 1].plot(filtered_data_gen["T"], filtered_data_gen["Value"], marker='o',label="Pyexocross Generated PF")
axs[0, 1].plot(filtered_data_scrape["T"], filtered_data_scrape["Value"], marker='o',label="Kurucz Scraped PF")
axs[0, 1].set_title('Cu II')
axs[0, 1].legend()
axs[0, 1].set_xlabel('Temperature')
axs[0, 1].set_ylabel('Partition Function')

scrape = "Kurucz/Mg/Mg/Kurucz/Mg__Kurucz.pf"
gen = "PyExoCross/TianyangXie/Mg/24Mg/Kurucz/24Mg__Kurucz.pf"
column_names = ['T', 'Value']
pf_scrape = pd.read_csv(scrape, delim_whitespace=True,names=column_names)
pf_gen = pd.read_csv(gen, delim_whitespace=True,names=column_names)
filtered_data_scrape = pf_scrape[pf_scrape["T"] <= 6000]
filtered_data_gen = pf_gen[pf_gen["T"] <= 6000]

axs[1, 0].plot(filtered_data_gen["T"], filtered_data_gen["Value"], marker='o',label="Pyexocross Generated PF")
axs[1, 0].plot(filtered_data_scrape["T"], filtered_data_scrape["Value"], marker='o',label="Kurucz Scraped PF")
axs[1, 0].set_title('Mg I')
axs[1, 0].legend()
axs[1, 0].set_xlabel('Temperature')
axs[1, 0].set_ylabel('Partition Function')



scrape = "Kurucz/Si/Si/Kurucz/Si__Kurucz.pf"
gen = "PyExoCross/TianyangXie/Si/28Si/Kurucz/28Si__Kurucz.pf"
column_names = ['T', 'Value']
pf_scrape = pd.read_csv(scrape, delim_whitespace=True,names=column_names)
pf_gen = pd.read_csv(gen, delim_whitespace=True,names=column_names)
filtered_data_scrape = pf_scrape[pf_scrape["T"] <= 6000]
filtered_data_gen = pf_gen[pf_gen["T"] <= 6000]

axs[1, 1].plot(filtered_data_gen["T"], filtered_data_gen["Value"], marker='o',label="Pyexocross Generated PF")
axs[1, 1].plot(filtered_data_scrape["T"], filtered_data_scrape["Value"], marker='o',label="Kurucz Scraped PF")
axs[1, 1].set_title('Si I')
axs[1, 1].legend()
axs[1, 1].set_xlabel('Temperature')
axs[1, 1].set_ylabel('Partition Function')



fig.suptitle('Pyexocross Generated PF vs Kurucz Scraped PF')
plt.tight_layout()
plt.legend()
plt.savefig("Plot/PF_scrape_gen_compare.png")
plt.show()



#The following code plots some of the scraped pf values in one figure.

# plt.figure(figsize=(10, 6))
#
#
# elements = ["Al","C","O","Cu","Ca","Mg"]
# for element in elements:
#     file_path = 'Kurucz/'+element+"/"+element+"/Kurucz/"+element+"__Kurucz.pf"
#     print(file_path)
#
#     column_names = ['T', 'Value']
#     data = pd.read_csv(file_path, delim_whitespace=True,names=column_names)
#
#     #filtered_data = data[data["T"] < 5000]
#     plt.plot(data["T"], data["Value"], marker='o',label=element)
#
# plt.xlabel('Temperature')
# plt.ylabel('Partition Function')
# plt.title('PF Plot')
# plt.legend()
# plt.savefig("Plot/Multi.png")
# #plt.savefig("Plot/Li_PF.png")
# plt.show()
#








