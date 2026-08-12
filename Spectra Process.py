import pandas as pd
import matplotlib.pyplot as plt

#This python file is for ploting the stick spectra calculated from the PyExoCross.
#By giving the stick file, we could plot the intensity.
#To map with the existing plots from other researches, we need to turn the PyExoCross generated wavenumber to wavelength.

filename = "Al_Spectra_Emission.stick"
column_names= ["Wavenumber","Intensity","1","2","3","4"]
df = pd.read_csv(filename, delim_whitespace=True, names=column_names)

df.drop(["1","2","3","4"],axis=1,inplace=True)

df["Wavelength"] = 1e7 / df["Wavenumber"]
df = df[['Wavelength',"Intensity"]]


filename = "Mg_Spectra_Emission.stick"
column_names= ["Wavenumber","Intensity","1","2","3","4"]
dfMg = pd.read_csv(filename, delim_whitespace=True, names=column_names)

dfMg.drop(["1","2","3","4"],axis=1,inplace=True)

dfMg["Wavelength"] = 1e7 / dfMg["Wavenumber"]
dfMg = dfMg[['Wavelength',"Intensity"]]



#The following code is to limit the range of the plot so that it could match the existing one in the same scale.
df["Intensity"] = df["Intensity"] * 6

df_filtered = df[(df['Wavelength'] >= 278) & (df['Wavelength'] <= 330)]
df_filtered_Mg = dfMg[(dfMg['Wavelength'] >= 278) & (dfMg['Wavelength'] <= 330)]


df_filtered_mg_278 = df_filtered_Mg[(df_filtered_Mg['Wavelength'] >= 278) & (df_filtered_Mg['Wavelength'] <= 280)]
peak_mg_wavelength = df_filtered_mg_278.loc[df_filtered_mg_278['Intensity'].idxmax(), 'Wavelength']
peak_al_wavelength = df_filtered.loc[df_filtered['Intensity'].idxmax(), 'Wavelength']

plt.figure(figsize=(10, 6))
plt.plot(df_filtered['Wavelength'], df_filtered['Intensity'])
plt.plot(df_filtered_Mg['Wavelength'], df_filtered_Mg['Intensity'],color="red")
plt.axvline(peak_mg_wavelength, color='green', linestyle='--', linewidth=2, label=f'Mg {peak_mg_wavelength:.2f} nm')
plt.axvline(peak_al_wavelength, color='orange', linestyle='--', linewidth=2, label=f'Al {peak_al_wavelength:.2f} nm')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity')
plt.title('Intensity vs. Wavelength')
plt.grid(True)
plt.legend()
plt.savefig("Plot/Match.png")
plt.show()