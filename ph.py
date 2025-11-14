#!/usr/bin/env python
# coding: utf-8

import numpy as np
from pandas.core.frame import AnyArrayLike

from inputs import PH_DEFAULTS, GRID_DEFAULTS
from store import *
from solve import *
from output import *


class PH:
    def __init__(
        self,
        DIV = PH_DEFAULTS["DIV"],
        totP = PH_DEFAULTS["totP"],
        Temp = PH_DEFAULTS["Temp"],
        nDIV = PH_DEFAULTS["nDIV"],
        nSiO2 = PH_DEFAULTS["nSiO2"],
        totnum = PH_DEFAULTS["totnum"],
        plot_flag = PH_DEFAULTS["plot_flag"],
        table_flag = PH_DEFAULTS["table_flag"],
        analytical_flag = PH_DEFAULTS["analytical_flag"],
        comparison = PH_DEFAULTS["comparison"],
    ):
                 
        self.DIV = DIV
        if self.DIV != 'Ca' and self.DIV != 'Mg' and self.DIV != 'Fe':
            print('Error: Enter DIV = "Ca" or "Mg" or "Fe"')
        self.totP = totP
        self.Temp = Temp
        self.nDIV = nDIV
        self.nSiO2 = nSiO2
        self.totnum = totnum
        self.plot_flag = plot_flag
        self.table_flag = table_flag
        self.analytical_flag = analytical_flag
        self.comparison = comparison

        if self.comparison == 'PCO2':
            if self.analytical_flag == True:
                self.pH_PCO2_an()
            else:
                self.pH_PCO2()

        elif self.comparison == 'P':
            self.pH_P()
        
        elif self.comparison == 'T':
            self.pH_T()


    def _run(self, PCO2s, setup_fn, solve_fn, save_fn):
            numQ = 3
            betas = np.array([-1, 0, 0.3])
            chems2 = chem_dict2(numQ, self.totnum)
            system, specs, solver = setup_fn()
            for i in range(numQ):
                j = 0
                beta = betas[i]
                while j < self.totnum:
                        PCO2 = PCO2s[j]
                        if beta == -1:
                            addDIVtot = 0
                            addSiO2   = 0
                        elif beta == 0:
                            addDIVtot = 1e-2
                            addSiO2   = 0
                        else:
                            addDIVtot = self.nDIV * weath_scaling(PCO2, self.Temp, beta=beta) / numden
                            addSiO2   = self.nSiO2 * addDIVtot
                        state = solve_fn(system, specs, solver, addDIVtot, addSiO2, PCO2, self.Temp, self.totP)
                        chems2 = save_fn(state, PCO2, chems2, i, j)
                        j = j + 1
            return chems2


    def pH_PCO2(self):
        '''
        Returns ocean pH as a function of PCO2 [bar] for Ca, Mg or Fe carbonate systems
        '''
        PCO2s = GRID_DEFAULTS["pco2s"](self.totnum) # bar

        if self.DIV == 'Ca':
            chems2 = self._run(PCO2s, setup_Ca, solve_Ca, save_chems2_Ca)
        elif self.DIV == 'Mg':
            chems2 = self._run(PCO2s, setup_Mg, solve_Mg, save_chems2_Mg)
        elif self.DIV == 'Fe':
            chems2 = self._run(PCO2s, setup_Fe, solve_Fe, save_chems2_Fe)
        else:
            return

        output_pH_PCO2(PCO2s, chems2, DIV = self.DIV, plot_flag = self.plot_flag, table_flag = self.table_flag)

        return

    def pH_PCO2_an(self, nDIV_fixed = 1):
        '''
        Returns analytical/numerical solutions of ocean pH as a function of PCO2 [bar]
        '''    
        if self.DIV != 'Ca':
            print('Error: Enter DIV = "Ca"')

        PCO2s = GRID_DEFAULTS["pco2s"](self.totnum) # bar
        
        if self.DIV == 'Ca': # Almost the same as self._run(), but slight differences for analytical

            numQ = 2
            addDIVtots = numden*np.linspace(0, 1e-2, num=numQ) # in units of 1000 * mol/dm3 = 1 mol/m3
            chems2 = chem_dict2(numQ,self.totnum)
            chems2_san = chem_dict2(numQ,self.totnum)
            chems2_an = chem_dict2(numQ,self.totnum)

            system, specs, solver = setup_Ca()
            
            logK3, logK9, logK16 = setup_an_Ca(self.Temp, self.totP)
            
            for i in range(numQ):
                j = 0
                addDIVtot = addDIVtots[i] / numden
                addSiO2   = 0

                while j < self.totnum:
                    PCO2 = PCO2s[j]
                    state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, self.Temp, self.totP)
                    chems2 = save_chems2_Ca(state, PCO2, chems2, i, j)
                    chems2_an, chems2_san = save_chems2_an_Ca(PCO2, logK3, logK9, logK16, nDIV_fixed,
                                                              chems2['Ca+2'][i][j], chems2_an, chems2_san, i, j)
                    j = j + 1

            output_pH_PCO2_an(PCO2s, chems2, chems2_an, chems2_san, DIV = self.DIV, nDIV_fixed = nDIV_fixed,
                              plot_flag = self.plot_flag, table_flag = self.table_flag)
        
        return

    def pH_P(self, PCO2 = 0.3e-3):
        '''
        Returns ocean pH as a function of P [bar]
        '''  
        if self.DIV != 'Ca':
            print('Error: Enter DIV = "Ca"')

        totPs = GRID_DEFAULTS["totps"](self.totnum) # bar

        ####################################
        
        if self.DIV == 'Ca': #Almost the same as self._run(), but slight differences for P

            numQ = 3
            betas = np.array([-1, 0, 0.3]) 
            chems2 = chem_dict2(numQ, self.totnum)

            system, specs, solver = setup_Ca()

            for i in range(numQ):
                j = 0
                beta = betas[i]
                while j < self.totnum:

                    totP = totPs[j]
                    if beta == -1:
                        addDIVtot = 0
                        addSiO2   = 0
                    elif beta == 0:
                        addDIVtot = 1e-2
                        addSiO2 = 0
                    else:
                        addDIVtot = self.nDIV * weath_scaling(PCO2, self.Temp, beta=beta) / numden
                        addSiO2   = self.nSiO2 * addDIVtot

                    state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, self.Temp, totP)
                    chems2 = save_chems2_Ca(state, PCO2, chems2, i, j)
                    j = j + 1
                    
        output_pH_P(totPs, chems2, DIV = self.DIV, plot_flag = self.plot_flag, table_flag = self.table_flag)

        return

    def pH_T(self, PCO2 = 0.3e-3):
        '''
        Returns ocean pH as a function of Temp [K]
        '''  
        if self.DIV != 'Ca':
            print('Error: Enter DIV = "Ca"')

        Temps = GRID_DEFAULTS["temps"](self.totnum) # bar

        ####################################
        
        if self.DIV == 'Ca':

            numQ = 3
            betas = np.array([-1, 0, 0.3]) 

            chems2 = chem_dict2(numQ, self.totnum)

            system, specs, solver = setup_Ca()

            for i in range(numQ): # Almost the same as self._run(), but slight differences for T
                j = 0

                beta = betas[i]
                while j < self.totnum:

                    Temp = Temps[j]
                    if beta == -1:
                        addDIVtot = 0
                        addSiO2   = 0
                    elif beta == 0:
                        addDIVtot = 1e-2
                        addSiO2   = 0
                    else:
                        addDIVtot = self.nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
                        addSiO2   = self.nSiO2 * addDIVtot

                    state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, self.totP)
                    chems2 = save_chems2_Ca(state, PCO2, chems2, i, j)
                    j = j + 1

            output_pH_T(Temps, chems2, DIV = self.DIV, plot_flag = self.plot_flag, table_flag = self.table_flag)
        
        return
