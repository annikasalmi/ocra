#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### plots_paper.py # contains functions to make plots available in the published paper

# Import ocra
from ocra import *
from ph import PH


## Ocean pH figures

### DIV decides the carbonate system: Ca, Mg or Fe

# Plot ocean pH as a function of PCO2 for the Ca-system
PH(DIV='Ca').pH_PCO2()

# Plot ocean pH as a function of PCO2 for the Mg-system
PH(DIV='Mg').pH_PCO2()

# Plot ocean pH as a function of PCO2 for the Fe-system
PH(DIV='Fe').pH_PCO2()

# Plot numerical and analytical solutions of ocean pH
PH().pH_PCO2_an()

# Plot ocean pH as a function of P
PH().pH_P()

# Plot ocean pH as a function of T
PH().pH_T()


## CCD figures

### nSiO2 = 1 includes silica and 0 excludes silica
### numQ1 is the number of steps in x-axis
### numQ2 is the number of steps in y-axis
### totnum is the number of steps in z-axis

# Plot Ca-CCD as a function of PCO2 and T with no silicates
CaCCD_PCO2_T(nSiO2 = 0, totnum = 20, numQ1 = 100, numQ2 = 100)

# Plot Ca-CCD as a function of PCO2 and T with silicates
CaCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)

# Plot Mg-CCD as a function of PCO2 and T with no silicates
MgCCD_PCO2_T(nSiO2 = 0, totnum = 20, numQ1 = 100, numQ2 = 100)

# Plot Mg-CCD as a function of PCO2 and T with silicates
MgCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)

# Plot Fe-CCD as a function of PCO2 and T with no silicates
FeCCD_PCO2_T(nSiO2 = 0, totnum = 20, numQ1 = 100, numQ2 = 100)

# Plot Fe-CCD as a function of PCO2 and T with silicates
FeCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)


## Partitioning figures

### DIV decides the carbonate system: Ca, Mg or Fe
### nSiO2 = 1 includes silica and 0 excludes silica
### Temp is temperature in kelvin

# Plot stable phases as a function of PCO2 for the Ca-system with no silicates
phases_PCO2(DIV='Ca', nSiO2 = 0, Temp = 310)

# Plot stable phases as a function of PCO2 for the Ca-system with silicates
phases_PCO2(DIV='Ca', nSiO2 = 1, Temp = 310)

# Plot stable phases as a function of PCO2 for the Mg-system with no silicates
phases_PCO2(DIV='Mg', nSiO2 = 0, Temp = 310)

# Plot stable phases as a function of PCO2 for the Mg-system with silicates
phases_PCO2(DIV='Mg', nSiO2 = 1, Temp = 310)

# Plot stable phases as a function of PCO2 for the Fe-system with no silicates
phases_PCO2(DIV='Fe', nSiO2 = 0, Temp = 310)

# Plot stable phases as a function of PCO2 for the Fe-system with silicates
phases_PCO2(DIV='Fe', nSiO2 = 1, Temp = 310)
