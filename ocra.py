#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### okra.py # contains functions to calculate ocean pH, CCD and stable phases

# Import libraries

import numpy as np
import pandas as pd

from store import *
from solve import *

# Calculate Ca-CCD as a function of PCO2 and T

def CaCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns Ca-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems3 = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    system, specs, solver = setup_Ca()
                
    for k, Temp in enumerate(Temps):

        for i, PCO2 in enumerate(PCO2s):
            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot
            j = 0
            while j < totnum:
                totP = totPs[j]
                state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                chems3 = save_chems3_Ca(state, PCO2, chems3, i, j, k)
                j = j + 1

    for k, _ in enumerate(Temps):
        for i in range(len(PCO2s)):
            nCarb_surf = chems3['Calcite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems3['Calcite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    output_CaCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV,
                        plot_flag = plot_flag, table_flag = table_flag)
        
    return


# Calculate Mg-CCD as a function of PCO2 and T

def MgCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns Mg-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems3 = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    system, specs, solver = setup_Mg()
                
    for k, Temp in enumerate(Temps):
        for i, PCO2 in enumerate(PCO2s):
            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot
            j = 0
            while j < totnum:
                totP = totPs[j]
                state = solve_Mg(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)  
                chems3 = save_chems3_Mg(state, PCO2, chems3, i, j, k)
                j = j + 1

    for k in range(len(Temps)):
        for i in range(len(PCO2s)):
            nCarb_surf = chems3['Magnesite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems3['Magnesite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    output_MgCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV,
                        plot_flag = plot_flag, table_flag = table_flag)
        
    return

# Calculate Fe-CCD as a function of PCO2 and T

def FeCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns Fe-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems3 = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    system, specs, solver = setup_Fe()
                
    for k, Temp in enumerate(Temps):
        for i, PCO2 in enumerate(PCO2s):
            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot
            j = 0
            while j < totnum:
                totP = totPs[j]
                state = solve_Fe(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                chems3 = save_chems3_Fe(state, PCO2, chems3, i, j, k)
                j = j + 1

    for k, _ in enumerate(Temps):
        for i, _ in enumerate(PCO2s):
            nCarb_surf = chems3['Siderite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems3['Siderite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    output_FeCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV,
                        plot_flag = plot_flag, table_flag = table_flag)
        
    return


# Calculate stable phases as a function of PCO2

def phases_PCO2 (DIV = 'Ca', Temp = 298, totP = 1, beta = 0.3, nDIV = 1, nSiO2 = 1, totnum = 100,
                table_flag = True, plot_flag = True):
    '''
    Returns stable phases as a function of PCO2 [bar]
    '''  
    if DIV != 'Ca' and DIV != 'Mg' and DIV != 'Fe':
        print('Error: Enter DIV = "Ca" or "Mg" or "Fe"')
    
    PCO2s = np.logspace(-8, -0.5, num=totnum) # bar
    chems1 = chem_dict1(totnum)
    
    if DIV == 'Ca':
        
        system, specs, solver = setup_Ca()

        j = 0
        while j < totnum:
            PCO2 = PCO2s[j]
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot
            state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)  
            chems1 = save_chems1_Ca(state, PCO2, chems1, j)
            j = j + 1
    
        df = pd.DataFrame({
            'Ca++': chems1['Ca+2'],
            'Calcite': chems1['Calcite'],
            'Silicates': chems1['Wollastonite'],
        }, index=PCO2s)
        
        output_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2,
                           table_flag = table_flag, plot_flag = plot_flag)
        
    
    elif DIV == 'Mg':
        
        system, specs, solver = setup_Mg()

        j = 0
        while j < totnum:
            PCO2 = PCO2s[j]
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot
            state = solve_Mg(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
            chems1 = save_chems1_Mg(state, PCO2, chems1, j)
            j = j + 1

        df = pd.DataFrame({
            'Mg++': chems1['Mg+2'],
            'Magnesite': chems1['Magnesite'],
            'Silicates': 2*chems1['Clino-Enstatite'],
        }, index=PCO2s)
        
        output_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2,
                           table_flag = table_flag, plot_flag = plot_flag)
    
    elif DIV == 'Fe':
        
        system, specs, solver = setup_Fe()

        j = 0
        while j < totnum:
            PCO2 = PCO2s[j]
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot
            state = solve_Fe(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP) 
            chems1 = save_chems1_Fe(state, PCO2, chems1, j)
            j = j + 1

        df = pd.DataFrame({
            'Fe++': chems1['Fe+2'],
            'Siderite': chems1['Siderite'],
            'Silicates': 2*chems1['Fayalite'],
        }, index=PCO2s)
        
        output_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2,
                           table_flag = table_flag, plot_flag = plot_flag)
        
    return

