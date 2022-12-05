from asyncio.windows_events import NULL
from tokenize import Double
import pythoncom
pythoncom.CoInitialize()

import os
import clr
import sys
import sqlite3
import pandas as pd
import numpy as np
import json
from gekko import GEKKO
import pandas as pd
from scipy.optimize import minimize, NonlinearConstraint
import System
from System import *
from System.IO import Directory, Path, File
from System import String, Environment, Array

mod = sys.modules[__name__]

dwsim_install_path = "C:/Users/lovpa/AppData/Local/DWSIM8"
sys.path.append(dwsim_install_path)
sys.path.append(os.getcwd() + "\\")
new_path = os.getcwd() + "\\"

data = pd.read_excel('data.xlsx')

conn = sqlite3.connect(new_path + 'Blank_db.db')
cur = conn.cursor()

clr.AddReference(new_path + "System.Buffers.dll")
clr.AddReference(dwsim_install_path + "/DWSIM.Automation.dll")
clr.AddReference(dwsim_install_path + "/DWSIM.Interfaces.dll")
clr.AddReference('DWSIM')

from DWSIM.Automation import Automation, Automation2

Directory.SetCurrentDirectory(dwsim_install_path)
interf = Automation2()

flowsheet = interf.LoadFlowsheet(new_path + "WORK1.dwxmz") 

STR_01 = flowsheet.GetFlowsheetSimulationObject("STR_01")
FURN_01 = flowsheet.GetFlowsheetSimulationObject("FURN_01")
LOSS_01 = flowsheet.GetFlowsheetSimulationObject("LOSS_01")
CRACK_01 = flowsheet.GetFlowsheetSimulationObject("CRACK_01")
STACK_01 = flowsheet.GetFlowsheetSimulationObject("STACK_01")
BFW_01 = flowsheet.GetFlowsheetSimulationObject("BFW_01")
HE_01 = flowsheet.GetFlowsheetSimulationObject("HE_01")

class CAL_GAS():

    def MAKE_ARR(self, k):
        with open(new_path + 'gas_comp.json', 'r', encoding="utf-8") as f:
            sample = json.load(f)
        arr = np.array(sample)
        arr = arr.flatten('F').reshape(arr.shape[1], arr.shape[0])
        self.arr_name = arr[0]
        self.arr_MWT = arr[1].astype(np.float16)
        self.arr_R1 = arr[2].astype(np.float16)
        self.arr_R2 = arr[3].astype(np.float16)
        self.arr_R3 = arr[4].astype(np.float16)
        self.arr_R4 = arr[5].astype(np.float16)
        self.arr_LHV = arr[6].astype(np.float16)
        self.arr_length = arr.shape[1]
        self.GTG_comp = [24.2331, 0.1261, 75.4859, 0.0007, 0.1367, 0.0176]                                                                 # GTG 연료 조성 [H2, CO, CH4, C2H6, C2H4, C2H2]   
        self.EXCESS_O2 = [0, 0, 0, 0, 0, data['F6_EXO2'][k]/100, 0, 0]                                          # 1호기 ~ 8호기
        
        # self.EXCESS_O2 = [data['F1_EXO2'][k]/100, data['F2_EXO2'][k]/100, data['F4_EXO2'][k]/100, data['F4_EXO2'][k]/100, data['F5_EXO2'][k]/100, data['F6_EXO2'][k]/100, data['F7_EXO2'][k]/100, data['F7_EXO2'][k]/100]                                          # 1호기 ~ 8호기
        # self.FUEL_MASS = [data['F1_FUEL_FLOW'][k], data['F2_FUEL_FLOW'][k], data['F3_FUEL_FLOW'][k], data['F4_FUEL_FLOW'][k], data['F5_FUEL_FLOW'][k], data['F6_FUEL_FLOW'][k], data['F7_FUEL_FLOW'][k], data['F8_FUEL_FLOW'][k]]
        self.FUEL_MASS = [0, 0, 0, 0, 0, data['F6_FUEL_FLOW'][k], 0, 0]
        
        self.AIRFLOW_FURN = [60.13174, 57.99801, 85.21120, 70.58341, 86.18465, data['F6_AIR_FLOW'][k], 88.16903, 79.14638]                               # 데이터
        self.STACK_T = [data["F1_STACK_T"][k],data["F2_STACK_T"][k],data["F3_STACK_T"][k],data["F4_STACK_T"][k],data["F5_STACK_T"][k],data["F6_STACK_T"][k],data["F7_STACK_T"][k],data["F8_STACK_T"][k]]
        self.GTG_FUEL = data['GTG_FUEL'][k]    
        self.GTG_O2_ANAL = data['GTG_Stack_O2'][k]/100
        var_t0 = data['F1_INLET_T'][k]
        var_t = data['F8_INLET_T'][k] - data['F1_INLET_T'][k]
        self.INLET_TEMP = [var_t0 + var_t*0/7, var_t0 + var_t*1/7,var_t0 + var_t*2/7,var_t0 + var_t*3/7,var_t0 + var_t*4/7,var_t0 + var_t*5/7,var_t0 + var_t*6/7,var_t0 + var_t*7/7]
        # self.BFW_FLOW = [data['F1_BFW_FLOW'][k],data['F2_BFW_FLOW'][k],data['F3_BFW_FLOW'][k],data['F4_BFW_FLOW'][k],data['F5_BFW_FLOW'][k],data['F6_BFW_FLOW'][k],data['F7_BFW_FLOW'][k],data['F8_BFW_FLOW'][k]]
        # self.STM_TEMP = [data['F1_STEAM_TEMP'][k],data['F2_STEAM_TEMP'][k],data['F3_STEAM_TEMP'][k],data['F4_STEAM_TEMP'][k],data['F5_STEAM_TEMP'][k],data['F6_STEAM_TEMP'][k],data['F7_STEAM_TEMP'][k],data['F8_STEAM_TEMP'][k]]
        self.BFW_FLOW = [0,0,0,0,0,data['F6_BFW_FLOW'][k],0,0]
        self.STM_TEMP = [0,0,0,0,0,data['F6_STEAM_TEMP'][k],0,0]
        
        self.OK = True
        self.OK2 = False
        
    def CAL_GTG(self):
        if self.GTG_FUEL < 100 :
            GTG_OUT_COMP = [0, 0, 0]
        else :
            GTG_MWT = (self.GTG_comp[0] * self.arr_MWT[0] + self.GTG_comp[1] * self.arr_MWT[23] + self.GTG_comp[2] * self.arr_MWT[1] 
            + self.GTG_comp[3] * self.arr_MWT[2] + self.GTG_comp[4] * self.arr_MWT[3] + self.GTG_comp[5] * self.arr_MWT[7]) / 100
            GTG_O2_GEN = (self.GTG_FUEL / GTG_MWT) * ((self.GTG_comp[0] / 100) * 0.5 + (self.GTG_comp[1] / 100) * 0.5 
            + (self.GTG_comp[2] / 100) * 2 + (self.GTG_comp[3] / 100) * 3.5 + (self.GTG_comp[4] / 100) * 3 + (self.GTG_comp[5] / 100) * 2.5)
            self.GTG_CO2_GEN = (self.GTG_FUEL / GTG_MWT) * ((self.GTG_comp[0] / 100) * 0 + (self.GTG_comp[1] / 100) * 1 
            + (self.GTG_comp[2] / 100) * 1 + (self.GTG_comp[3] / 100) * 2 + (self.GTG_comp[4] / 100) * 2 + (self.GTG_comp[5] / 100) * 2)
            self.GTG_H2O_GEN = (self.GTG_FUEL / GTG_MWT) * ((self.GTG_comp[0] / 100) * 1 + (self.GTG_comp[1] / 100) * 0 
            + (self.GTG_comp[2] / 100) * 2 + (self.GTG_comp[3] / 100) * 3 + (self.GTG_comp[4] / 100) * 2 + (self.GTG_comp[5] / 100) * 1)
            GTG_N2_GEN = GTG_O2_GEN * ( 0.79 / 0.21)
            
            self.GTG_N2 = ((self.GTG_O2_ANAL * (GTG_N2_GEN + self.GTG_CO2_GEN + self.GTG_H2O_GEN)) / (0.21 - self.GTG_O2_ANAL)) * 0.79 + GTG_N2_GEN
            self.GTG_O2 = ((self.GTG_O2_ANAL * (GTG_N2_GEN + self.GTG_CO2_GEN + self.GTG_H2O_GEN)) / (0.21 - self.GTG_O2_ANAL)) * 0.21
            
            GTG_OUT_CO2 = self.GTG_CO2_GEN / (self.GTG_O2 + self.GTG_CO2_GEN + self.GTG_H2O_GEN + self.GTG_N2)
            GTG_OUT_H2O = self.GTG_H2O_GEN / (self.GTG_O2 + self.GTG_CO2_GEN + self.GTG_H2O_GEN + self.GTG_N2)
            GTG_OUT_N2 = self.GTG_N2 / (self.GTG_O2 + self.GTG_CO2_GEN + self.GTG_H2O_GEN + self.GTG_N2)
            GTG_OUT_COMP = [GTG_OUT_CO2, GTG_OUT_H2O, GTG_OUT_N2]

        return GTG_OUT_COMP

    def CAL_FUEL_MWT(self):
        arr_COMP1 = np.array([22.0841, 67.8167, 5.1084, 1.0459, 1.4997, 0, 0.0345, 0.0138])                         # COMP. [H2, CH4, C2H6, C2H4, C3H8, △-C3H6, C3H6, C2H2]
        arr_COMP2 = np.array([0.0408, 0, 0.1303, 0.1897, 0.1044, 0.1532, 0.4364])                                   # COMP. [i-C4*, PD, n-C4*, T-2-C4-, 1-C4-, i-C4-, C-2-C4-]
        arr_COMP3 = np.array([0.0003, 0.0613, 0.0679, 0.0010, 0.2205, 0.2606, 0.0203, 0.1508])                      # COMP. [22-DM-PROP, i-C5*, 12-BUTAD, MA, 13-BUTAD, 3M-1-C4-, VA, EA]
        arr_COMP4 = np.array([0.1185, 0.0022, 0.6421, 0.0768])                                                      # COMP. [CO, CO2, N2, n-C5*]
        arr_COMP5 = np.array([0.0084, 0.0240, 0.0003, 0.0012, 0.0001, 0.0001])                                      # COMP. [HEXANE, BENZENE, HEPTANE, TOLUENE, EB, XYLENE]

        self.arr_comp = np.block([arr_COMP1,arr_COMP2,arr_COMP3,arr_COMP4,arr_COMP5])
        
        return (np.dot(self.arr_comp, self.arr_MWT))/100

    def CAL_LHV(self, x):
        FUEL_MWT = self.CAL_FUEL_MWT()
        LHV = 0
        for i in range(self.arr_length):
            LHV_ = self.arr_comp[i] * self.arr_LHV[i] * self.arr_MWT[i]
            LHV = LHV + LHV_
        FUEL_LHV = self.FUEL_MASS[x-1] * LHV / (100* FUEL_MWT)
        
        return FUEL_LHV

    def CAL_GEN(self, x):
        FUEL_MOLE = self.FUEL_MASS[x-1] / self.CAL_FUEL_MWT()
        self.GEN_O2 = 0
        self.GEN_CO2 = 0
        self.GEN_H2O = 0

        for i in range(self.arr_length):
            GEN_O2_ = FUEL_MOLE * (self.arr_comp[i]/100) * self.arr_R2[i]
            self.GEN_O2 = self.GEN_O2 + GEN_O2_
        for i in range(self.arr_length):
            GEN_CO2_ = FUEL_MOLE * (self.arr_comp[i]/100) * self.arr_R3[i]
            self.GEN_CO2 = self.GEN_CO2 + GEN_CO2_
        for i in range(self.arr_length):
            GEN_H2O_ = FUEL_MOLE * (self.arr_comp[i]/100) * self.arr_R4[i]
            self.GEN_H2O = self.GEN_H2O + GEN_H2O_
        self.GEN_N2 = FUEL_MOLE * (self.arr_comp[25]/100)
        self.GEN_X = (self.GEN_CO2 + self.GEN_N2 + self.GEN_H2O)/self.GEN_O2
       

    def CAL_AIR_RATIO(self, x):
        self.CAL_GEN(x)
        self.AIR_REQ = self.GEN_O2 * 22.4 / self.HAIR_comp[3]    
        self.CO2_REQ = self.GEN_CO2 * 22.4 + self.AIR_REQ * self.HAIR_comp[0] 
        self.N2_REQ = self.AIR_REQ * self.HAIR_comp[2] + self.GEN_N2 * 22.4
        self.H2O_REQ = self.GEN_H2O * 22.4 + self.AIR_REQ * self.HAIR_comp[1] 
        
        m = GEKKO()
        k = m.Var()
        m.Equations([k / (k + self.CO2_REQ + self.N2_REQ + self.H2O_REQ) == self.EXCESS_O2[x-1]])
        m.solve(disp=False)
        self.O2_REQ = k.value[0]
        
        # EX_AIR_REQ = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (self.EXCESS_O2[x-1]) / (self.HAIR_comp[3] - self.EXCESS_O2[x-1])
        EX_AIR_REQ = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ + self.O2_REQ ) - self.AIR_REQ
        EX_AIR_RATIO = 100 * EX_AIR_REQ / self.AIR_REQ
        self.HAIR_REQ = self.AIR_REQ + EX_AIR_REQ
        self.CO2_FRAC = self.CO2_REQ / self.HAIR_REQ
        self.H2O_FRAC = self.H2O_REQ / self.HAIR_REQ
        self.O2_FRAC = self.O2_REQ / self.HAIR_REQ
        self.N2_FRAC = self.N2_REQ / self.HAIR_REQ
        

        return EX_AIR_RATIO

    def CAL_H_HAIR(self, x):
        self.CAL_AIR_RATIO(x)
        # HAIR_CO2_W = (self.HAIR_REQ * self.HAIR_comp[0] * 44.010) / 22.4
        # HAIR_H2O_W = (self.HAIR_REQ * self.HAIR_comp[1] * 18.015) / 22.4
        # HAIR_N2_W = (self.HAIR_REQ * self.HAIR_comp[2] * 28.013) / 22.4
        # self.HAIR_O2_W = (self.HAIR_REQ * self.HAIR_comp[3] * 31.999) / 22.4

        HAIR_CO2_W = (self.CO2_REQ * 44.010) / 22.4
        HAIR_H2O_W = (self.H2O_REQ * 18.015) / 22.4
        HAIR_N2_W = (self.N2_REQ * 28.013) / 22.4
        self.HAIR_O2_W = (self.O2_REQ * 31.999) / 22.4

        self.HAIR_TOT = (HAIR_CO2_W + HAIR_H2O_W + HAIR_N2_W + self.HAIR_O2_W)
        HAIR_CO2_WT_FRAC = HAIR_CO2_W / self.HAIR_TOT
        HAIR_H2O_WT_FRAC = HAIR_H2O_W / self.HAIR_TOT
        HAIR_N2_WT_FRAC = HAIR_N2_W / self.HAIR_TOT
        HAIR_O2_WT_FRAC = self.HAIR_O2_W / self.HAIR_TOT

        T_1 = 520                            # ['R] 
        T_2 = (self.INLET_TEMP[x-1]) * 1.8 + 32 + 460        # ['R]

        H_CO2 = (0.114433*(T_2 - T_1)) + (1.011E-04 *(T_2**2 - T_1**2)) + (-2.649E-08 * (T_2**3 - T_1**3)) + (3.471E-12 *(T_2**4 - T_1**4)) + (-1.314E-16 *(T_2**5 - T_1**5))
        H_H2O = (0.457392*(T_2 - T_1)) + (-5.251E-05 *(T_2**2 - T_1**2)) + (6.459E-08 * (T_2**3 - T_1**3)) + (-2.028E-11 *(T_2**4 - T_1**4)) + (2.361E-15 *(T_2**5 - T_1**5))
        H_N2 = (0.255204*(T_2 - T_1)) + (-1.779E-05 *(T_2**2 - T_1**2)) + (1.589E-08 * (T_2**3 - T_1**3)) + (-3.220E-12 *(T_2**4 - T_1**4)) + (1.589E-16 *(T_2**5 - T_1**5))
        H_O2 = (0.227486*(T_2 - T_1)) + (-3.731E-05 *(T_2**2 - T_1**2)) + (4.830E-08 * (T_2**3 - T_1**3)) + (-1.852E-11 *(T_2**4 - T_1**4)) + (2.475E-15 *(T_2**5 - T_1**5))

        WEI_AVG_H = (H_CO2 * HAIR_CO2_WT_FRAC + H_H2O * HAIR_H2O_WT_FRAC + H_N2 * HAIR_N2_WT_FRAC + H_O2 * HAIR_O2_WT_FRAC) * (0.252164 / 0.453592)     # [kcal/kg]

        return WEI_AVG_H * self.HAIR_TOT


    def CAL_FURN_GAS(self, x):

        self.CAL_H_HAIR(x)
        FURN_GAS_O2 =  self.HAIR_O2_W*0.15 - self.GEN_O2 * 32
        FURN_GAS_N2 = (self.GEN_N2 + self.GEN_O2*31.9998*self.HAIR_comp[2] / (self.HAIR_comp[3]*28.0134)) 
        FURN_GAS_H2O = (self.GEN_H2O + self.GEN_O2*31.9998*self.HAIR_comp[1]/(self.HAIR_comp[3]*18.01534)) 
        FURN_GAS_CO2 = (self.GEN_CO2 + self.GEN_O2*31.9998*self.HAIR_comp[0]/(self.HAIR_comp[3]*44.00995)) 
        
        FURN_GAS_MASS = FURN_GAS_O2*31.9998 + FURN_GAS_N2*28.0134 + FURN_GAS_H2O*18.01534 + FURN_GAS_CO2*44.00995

        FURN_GAS_O2_FRAC = FURN_GAS_O2 / (FURN_GAS_O2 + FURN_GAS_N2 + FURN_GAS_H2O + FURN_GAS_CO2)
        FURN_GAS_N2_FRAC = FURN_GAS_N2 / (FURN_GAS_O2 + FURN_GAS_N2 + FURN_GAS_H2O + FURN_GAS_CO2)
        FURN_GAS_CO2_FRAC = FURN_GAS_CO2 / (FURN_GAS_O2 + FURN_GAS_N2 + FURN_GAS_H2O + FURN_GAS_CO2)
        FURN_GAS_H2O_FRAC = FURN_GAS_H2O / (FURN_GAS_O2 + FURN_GAS_N2 + FURN_GAS_H2O + FURN_GAS_CO2)
        FURN_GAS = np.array([FURN_GAS_MASS ,FURN_GAS_O2_FRAC, FURN_GAS_N2_FRAC, FURN_GAS_H2O_FRAC, FURN_GAS_CO2_FRAC])
        
        return FURN_GAS

    def SOLVE_COMP(self):
        
        # m = GEKKO(remote=False) 
        # m.options.SOLVER = 3        # SOLVER = 1 -> 'APOPT'(solve MINLP) // SOLVER = 3 -> 'IPOPT'(solve NLP)

        arr_GTG = self.CAL_GTG()

        # #-----MV-----

        # # GTG_RATIO = m.Var(lb=0, ub=1)
        # EX_O2_var1 = m.Var()
        # EX_O2_var2 = m.Var()
        # EX_O2_var3 = m.Var()
        # EX_O2_var4 = m.Var()
        # EX_O2_var5 = m.Var()
        # EX_O2_var6 = m.Var()
        # EX_O2_var7 = m.Var()
        # EX_O2_var8 = m.Var()
        # #-------------
        GTG_RATIO = 0.9
        self.HAIR_comp = [arr_GTG[0]*GTG_RATIO, arr_GTG[1]*GTG_RATIO, GTG_RATIO*(arr_GTG[2]-0.79) + 0.79, GTG_RATIO*(self.GTG_O2_ANAL - 0.21) + 0.21]   #CO2, H2O, N2, O2

        # X_O2 = []
        # X_HAIR = []
        # for i in range(8):
        #     self.CAL_H_HAIR(i+1)
        #     try:
        #         HAIR_O2_WX = self.GEN_O2 * 31.999 * (1/self.EXCESS_O2[i] - 1/self.HAIR_comp[3] + self.GEN_X + self.HAIR_comp[0]/self.HAIR_comp[3] + self.HAIR_comp[1]/self.HAIR_comp[3] + self.HAIR_comp[2]/self.HAIR_comp[3])/(1/self.EXCESS_O2[i] - 1/self.HAIR_comp[3])
        #     except ZeroDivisionError:
        #         self.OK = False
        #         print("Error : ZeroDivisionError!")
        #     X_O2.append(HAIR_O2_WX)
        #     HAIR_TOTX = (self.HAIR_comp[0] * 44.010 + self.HAIR_comp[1] * 18.015 + self.HAIR_comp[2] * 28.013 + self.HAIR_comp[3] * 31.999)*HAIR_O2_WX/(self.HAIR_comp[3]*31.999)
        #     X_HAIR.append(HAIR_TOTX)
        # X_O2 = np.sum(X_O2)
        # Y_AIR = np.sum(X_HAIR)
        # Y_O2 = self.GTG_O2 * 32
        # EX_GTG_AIR = (self.GTG_CO2_GEN * 44.01 + self.GTG_H2O_GEN * 18.01528 + self.GTG_O2 * 32 + self.GTG_N2 * 28.014)
        # X_AIR = (self.GTG_CO2_GEN * 44.01 + self.GTG_H2O_GEN * 18.01528 + self.GTG_O2 * 32 + self.GTG_N2 * 28.014) + (X_O2 - Y_O2) * 100 / 23

        # self.CAL_AIR_RATIO(1)
        # EX_AIR_REQ_1 = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (EX_O2_var1/100) / (self.HAIR_comp[3] - EX_O2_var1/100)
        # HOT_AIR_TOT_V1 = EX_AIR_REQ_1 + self.AIR_REQ
        # HOT_AIR_CO2_1 = (HOT_AIR_TOT_V1 * self.HAIR_comp[0]) * 44.00995 / 22.4
        # HOT_AIR_H2O_1 = (HOT_AIR_TOT_V1 * self.HAIR_comp[1]) * 18.01534 / 22.4
        # HOT_AIR_N2_1 = (HOT_AIR_TOT_V1 * self.HAIR_comp[2]) * 28.0134 / 22.4
        # HOT_AIR_O2_1 = (HOT_AIR_TOT_V1 * self.HAIR_comp[3]) * 31.9988 / 22.4
        # HOT_AIR_TOT_W1 = (HOT_AIR_CO2_1 + HOT_AIR_H2O_1 + HOT_AIR_N2_1 + HOT_AIR_O2_1)

        # self.CAL_AIR_RATIO(2)
        # EX_AIR_REQ_2 = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (EX_O2_var2/100) / (self.HAIR_comp[3] - EX_O2_var2/100)
        # HOT_AIR_TOT_V2 = EX_AIR_REQ_2 + self.AIR_REQ
        # HOT_AIR_CO2_2 = (HOT_AIR_TOT_V2 * self.HAIR_comp[0]) * 44.00995 / 22.4
        # HOT_AIR_H2O_2 = (HOT_AIR_TOT_V2 * self.HAIR_comp[1]) * 18.01534 / 22.4
        # HOT_AIR_N2_2 = (HOT_AIR_TOT_V2 * self.HAIR_comp[2]) * 28.0134 / 22.4
        # HOT_AIR_O2_2 = (HOT_AIR_TOT_V2 * self.HAIR_comp[3]) * 31.9988 / 22.4
        # HOT_AIR_TOT_W2 = (HOT_AIR_CO2_2 + HOT_AIR_H2O_2 + HOT_AIR_N2_2 + HOT_AIR_O2_2)

        # self.CAL_AIR_RATIO(3)
        # EX_AIR_REQ_3 = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (EX_O2_var3/100) / (self.HAIR_comp[3] - EX_O2_var3/100)
        # HOT_AIR_TOT_V3 = EX_AIR_REQ_3 + self.AIR_REQ
        # HOT_AIR_CO2_3 = (HOT_AIR_TOT_V3 * self.HAIR_comp[0]) * 44.00995 / 22.4
        # HOT_AIR_H2O_3 = (HOT_AIR_TOT_V3 * self.HAIR_comp[1]) * 18.01534 / 22.4
        # HOT_AIR_N2_3 = (HOT_AIR_TOT_V3 * self.HAIR_comp[2]) * 28.0134 / 22.4
        # HOT_AIR_O2_3 = (HOT_AIR_TOT_V3 * self.HAIR_comp[3]) * 31.9988 / 22.4
        # HOT_AIR_TOT_W3 = (HOT_AIR_CO2_3 + HOT_AIR_H2O_3 + HOT_AIR_N2_3 + HOT_AIR_O2_3)

        # self.CAL_AIR_RATIO(4)
        # EX_AIR_REQ_4 = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (EX_O2_var4/100) / (self.HAIR_comp[3] - EX_O2_var4/100)
        # HOT_AIR_TOT_V4 = EX_AIR_REQ_4 + self.AIR_REQ
        # HOT_AIR_CO2_4 = (HOT_AIR_TOT_V4 * self.HAIR_comp[0]) * 44.00995 / 22.4
        # HOT_AIR_H2O_4 = (HOT_AIR_TOT_V4 * self.HAIR_comp[1]) * 18.01534 / 22.4
        # HOT_AIR_N2_4 = (HOT_AIR_TOT_V4 * self.HAIR_comp[2]) * 28.0134 / 22.4
        # HOT_AIR_O2_4 = (HOT_AIR_TOT_V4 * self.HAIR_comp[3]) * 31.9988 / 22.4
        # HOT_AIR_TOT_W4 = (HOT_AIR_CO2_4 + HOT_AIR_H2O_4 + HOT_AIR_N2_4 + HOT_AIR_O2_4)

        # self.CAL_AIR_RATIO(5)
        # EX_AIR_REQ_5 = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (EX_O2_var5/100) / (self.HAIR_comp[3] - EX_O2_var5/100)
        # HOT_AIR_TOT_V5 = EX_AIR_REQ_5 + self.AIR_REQ
        # HOT_AIR_CO2_5 = (HOT_AIR_TOT_V5 * self.HAIR_comp[0]) * 44.00995 / 22.4
        # HOT_AIR_H2O_5 = (HOT_AIR_TOT_V5 * self.HAIR_comp[1]) * 18.01534 / 22.4
        # HOT_AIR_N2_5 = (HOT_AIR_TOT_V5 * self.HAIR_comp[2]) * 28.0134 / 22.4
        # HOT_AIR_O2_5 = (HOT_AIR_TOT_V5 * self.HAIR_comp[3]) * 31.9988 / 22.4
        # HOT_AIR_TOT_W5 = (HOT_AIR_CO2_5 + HOT_AIR_H2O_5 + HOT_AIR_N2_5 + HOT_AIR_O2_5)

        # self.CAL_AIR_RATIO(6)
        # EX_AIR_REQ_6 = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (EX_O2_var6/100) / (self.HAIR_comp[3] - EX_O2_var6/100)
        # HOT_AIR_TOT_V6 = EX_AIR_REQ_6 + self.AIR_REQ
        # HOT_AIR_CO2_6 = (HOT_AIR_TOT_V6 * self.HAIR_comp[0]) * 44.00995 / 22.4
        # HOT_AIR_H2O_6 = (HOT_AIR_TOT_V6 * self.HAIR_comp[1]) * 18.01534 / 22.4
        # HOT_AIR_N2_6 = (HOT_AIR_TOT_V6 * self.HAIR_comp[2]) * 28.0134 / 22.4
        # HOT_AIR_O2_6 = (HOT_AIR_TOT_V6 * self.HAIR_comp[3]) * 31.9988 / 22.4
        # HOT_AIR_TOT_W6 = (HOT_AIR_CO2_6 + HOT_AIR_H2O_6 + HOT_AIR_N2_6 + HOT_AIR_O2_6)

        # self.CAL_AIR_RATIO(7)
        # EX_AIR_REQ_7 = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (EX_O2_var7/100) / (self.HAIR_comp[3] - EX_O2_var7/100)
        # HOT_AIR_TOT_V7 = EX_AIR_REQ_7 + self.AIR_REQ
        # HOT_AIR_CO2_7 = (HOT_AIR_TOT_V7 * self.HAIR_comp[0]) * 44.00995 / 22.4
        # HOT_AIR_H2O_7 = (HOT_AIR_TOT_V7 * self.HAIR_comp[1]) * 18.01534 / 22.4
        # HOT_AIR_N2_7 = (HOT_AIR_TOT_V7 * self.HAIR_comp[2]) * 28.0134 / 22.4
        # HOT_AIR_O2_7 = (HOT_AIR_TOT_V7 * self.HAIR_comp[3]) * 31.9988 / 22.4
        # HOT_AIR_TOT_W7 = (HOT_AIR_CO2_7 + HOT_AIR_H2O_7 + HOT_AIR_N2_7 + HOT_AIR_O2_7)

        # self.CAL_AIR_RATIO(8)
        # EX_AIR_REQ_8 = (self.CO2_REQ + self.N2_REQ + self.H2O_REQ) * (EX_O2_var8/100) / (self.HAIR_comp[3] - EX_O2_var8/100)
        # HOT_AIR_TOT_V8 = EX_AIR_REQ_8 + self.AIR_REQ
        # HOT_AIR_CO2_8 = (HOT_AIR_TOT_V8 * self.HAIR_comp[0]) * 44.00995 / 22.4
        # HOT_AIR_H2O_8 = (HOT_AIR_TOT_V8 * self.HAIR_comp[1]) * 18.01534 / 22.4
        # HOT_AIR_N2_8 = (HOT_AIR_TOT_V8 * self.HAIR_comp[2]) * 28.0134 / 22.4
        # HOT_AIR_O2_8 = (HOT_AIR_TOT_V8 * self.HAIR_comp[3]) * 31.9988 / 22.4
        # HOT_AIR_TOT_W8 = (HOT_AIR_CO2_8 + HOT_AIR_H2O_8 + HOT_AIR_N2_8 + HOT_AIR_O2_8)
        
        # # obj = EX_GTG_AIR / X_AIR
        # obj = (EX_O2_var6 - self.EXCESS_O2[5])**2 
        # #-----Contraints-----
        # m.Equation(X_AIR == Y_AIR)
        # m.Equation(HOT_AIR_TOT_W1 == (self.AIRFLOW_FURN[0] * 1000))
        # m.Equation(HOT_AIR_TOT_W2 == (self.AIRFLOW_FURN[1] * 1000))
        # m.Equation(HOT_AIR_TOT_W3 == (self.AIRFLOW_FURN[2] * 1000))
        # m.Equation(HOT_AIR_TOT_W4 == (self.AIRFLOW_FURN[3] * 1000))
        # m.Equation(HOT_AIR_TOT_W5 == (self.AIRFLOW_FURN[4] * 1000))
        # m.Equation(HOT_AIR_TOT_W6 == (self.AIRFLOW_FURN[5] * 1000))
        # m.Equation(HOT_AIR_TOT_W7 == (self.AIRFLOW_FURN[6] * 1000))
        # m.Equation(HOT_AIR_TOT_W8 == (self.AIRFLOW_FURN[7] * 1000))
        # m.Equation(obj > 0)
        # # m.Equation(EX_O2_var1 == self.EXCESS_O2[0])
        # # m.Equation(EX_O2_var2 == self.EXCESS_O2[1])
        # # m.Equation(EX_O2_var3 == self.EXCESS_O2[2])
        # # m.Equation(EX_O2_var4 == self.EXCESS_O2[3])
        # # m.Equation(EX_O2_var5 == self.EXCESS_O2[4])
        # # m.Equation(EX_O2_var6 == self.EXCESS_O2[5])
        # # m.Equation(EX_O2_var7 == self.EXCESS_O2[6])
        # # m.Equation(EX_O2_var8 == self.EXCESS_O2[7])
        # m.Equation(EX_O2_var1 > 0)
        # m.Equation(EX_O2_var2 > 0)
        # m.Equation(EX_O2_var3 > 0)
        # m.Equation(EX_O2_var4 > 0)
        # m.Equation(EX_O2_var5 > 0)
        # m.Equation(EX_O2_var6 > 0)
        # m.Equation(EX_O2_var7 > 0)
        # m.Equation(EX_O2_var8 > 0)
        # #--------------------
    
        # m.Minimize(obj)
 
        # try:
        #     m.solve(disp=True)
        # except:
        #     # self.SOLVE_COMP_RE()
        #     self.OK = False
        #     pass
        # if (self.OK == True) :
        #     print("Complete optimization!")
        #     # print("GTG/(GTG+fresh air) RATIO : " + str(GTG_RATIO[0]))
        #     # print("------계산 분해로 O2값------")
        #     # print("1호기 O2 계산값 : " + str(EX_O2_var1[0]))
        #     # print("2호기 O2 계산값 : " + str(EX_O2_var2[0]))
        #     # print("3호기 O2 계산값 : " + str(EX_O2_var3[0]))
        #     # print("4호기 O2 계산값 : " + str(EX_O2_var4[0]))
        #     # print("5호기 O2 계산값 : " + str(EX_O2_var5[0]))
        #     # print("6호기 O2 계산값 : " + str(EX_O2_var6[0]))
        #     # print("7호기 O2 계산값 : " + str(EX_O2_var7[0]))
        #     # print("8호기 O2 계산값 : " + str(EX_O2_var8[0]))
        #     # print("---------------------------")
        #     GTG_RATIO = GTG_RATIO[0]
        #     self.HAIR_comp = [arr_GTG[0]*GTG_RATIO, arr_GTG[1]*GTG_RATIO, GTG_RATIO*(arr_GTG[2]-0.79) + 0.79, GTG_RATIO*(self.GTG_O2_ANAL - 0.21) + 0.21]
        # else : 
        #     pass

    def SOLVE_COMP_RE(self):
        
        m = GEKKO(remote=False) 
        m.options.SOLVER = 3        # SOLVER = 1 -> 'APOPT'(solve MINLP) // SOLVER = 3 -> 'IPOPT'(solve NLP)

        arr_GTG = self.CAL_GTG()

        #-----MV-----
        GTG_RATIO = m.Var(lb=0, ub=1)
        #-------------

        self.HAIR_comp = [arr_GTG[0]*GTG_RATIO, arr_GTG[1]*GTG_RATIO, GTG_RATIO*(arr_GTG[2]-0.79) + 0.79, GTG_RATIO*(self.GTG_O2_ANAL - 0.21) + 0.21]
        
        X_O2 = []
        X_HAIR = []
        for i in range(8):
            self.CAL_H_HAIR(i+1)
            try:
                HAIR_O2_WX = self.GEN_O2 * 31.999 * (1/self.EXCESS_O2[i] - 1/self.HAIR_comp[3] + self.GEN_X + self.HAIR_comp[0]/self.HAIR_comp[3] + self.HAIR_comp[1]/self.HAIR_comp[3] + self.HAIR_comp[2]/self.HAIR_comp[3])/(1/self.EXCESS_O2[i] - 1/self.HAIR_comp[3])
            except ZeroDivisionError:
                self.OK = False
                print("Error : ZeroDivisionError!")
            X_O2.append(HAIR_O2_WX)
            HAIR_TOTX = (self.HAIR_comp[0] * 44.010 + self.HAIR_comp[1] * 18.015 + self.HAIR_comp[2] * 28.013 + self.HAIR_comp[3] * 31.999)*HAIR_O2_WX/(self.HAIR_comp[3]*31.999)
            X_HAIR.append(HAIR_TOTX)
        X_O2 = np.sum(X_O2)
        Y_AIR = np.sum(X_HAIR)
        Y_O2 = self.GTG_O2 * 32
        EX_GTG_AIR = (self.GTG_CO2_GEN * 44.01 + self.GTG_H2O_GEN * 18.01528 + self.GTG_O2 * 32 + self.GTG_N2 * 28.014)
        X_AIR = (self.GTG_CO2_GEN * 44.01 + self.GTG_H2O_GEN * 18.01528 + self.GTG_O2 * 32 + self.GTG_N2 * 28.014) + (X_O2 - Y_O2) * 100 / 23

        
        obj = EX_GTG_AIR / X_AIR
        #-----Contraints-----
        m.Equation(X_AIR == Y_AIR)
        m.Equation(obj > 0)
        #--------------------
    
        m.Minimize(obj)
 
        try:
            m.solve(disp=False)
            self.OK2 = True
        except:
            self.OK = False

            pass
        if (self.OK2 == True) :
            GTG_RATIO = GTG_RATIO[0]
            self.HAIR_comp = [arr_GTG[0]*GTG_RATIO, arr_GTG[1]*GTG_RATIO, GTG_RATIO*(arr_GTG[2]-0.79) + 0.79, GTG_RATIO*(self.GTG_O2_ANAL - 0.21) + 0.21]
        else : 
            pass


    def SOLVE(self, x):
        for i in range(2900, 2901):
            self.MAKE_ARR(i+1)
            if self.GTG_FUEL < 100:
                self.HAIR_comp = [0, 0, 0.79, 0.21]
                self.CAL_H_HAIR(x)
                BFW_P = data["BFW_PRESS"][i+1]
                BFW_T = data["BFW_TEMP"][i+1] + 273
                STR_01.SetMassFlow(f"{self.HAIR_TOT} kg/h")
                STR_01.SetOverallComposition(Array[float]((self.CO2_FRAC, self.N2_FRAC, self.O2_FRAC, self.H2O_FRAC)))
                BFW_01.SetPressure(f"{BFW_P} kgf/cm2")
                BFW_01.SetTemperature(BFW_T)
                BFW_01.SetMassFlow(f"{self.BFW_FLOW[x-1] * 1000} kg/h")
                HE_01.set_ColdSideOutletTemperature(self.STM_TEMP[x-1]+ 273.15)

                LHV = self.CAL_LHV(x) * 4.1868 / 3600
                set_q = LHV + (self.CAL_H_HAIR(x) * 4.1868 / 3600)      # [kW]
                FURN_01.set_DeltaQ(set_q)
                CRACK_01.set_OutletTemperature(self.STACK_T[x-1] + 273)
                interf.CalculateFlowsheet2(flowsheet)
                
                R_FLUE_T2 = STACK_01.GetTemperature() - 273.15
                R_FLUE_O2 = STACK_01.GetOverallComposition()[2] * 100

                LOSS = LOSS_01.get_DeltaQ() - STACK_01.GetOverallComposition()[3] * STACK_01.GetMolarFlow() * 3.6 * 18.0152 * 2260 / 3600
                RAD_LOSS = 0.025
                EFF_HL = 1 - (LOSS/LHV + RAD_LOSS)
                self.OK = True

                arr = [R_FLUE_T2, R_FLUE_O2, LHV, LOSS, EFF_HL, self.OK]
                tup_arr = tuple(arr)
                sql = "INSERT INTO FURN_6_2 VALUES (?,?,?,?,?,?)"
                cur.execute(sql, tup_arr)
                conn.commit()
                print(str(i+1) + " -> success!" )
                interf.SaveFlowsheet(flowsheet, new_path + "dump.dwxmz", True)

            else:
                self.SOLVE_COMP()
                if (self.OK == True) or (self.OK2 == True) :
                    Q1 = self.CAL_H_HAIR(x)
                    BFW_P = data["BFW_PRESS"][i+1]
                    BFW_T = data["BFW_TEMP"][i+1] + 273
                    # STR_01.SetMassFlow(f"{FURN_GAS[0]} kg/h")
                    # STR_01.SetOverallComposition(Array[float]((FURN_GAS[4], FURN_GAS[2], FURN_GAS[1], FURN_GAS[3])))
                    STR_01.SetMassFlow(f"{self.HAIR_TOT} kg/h")
                    STR_01.SetOverallComposition(Array[float]((self.CO2_FRAC, self.N2_FRAC, self.O2_FRAC, self.H2O_FRAC)))
                    BFW_01.SetPressure(f"{BFW_P} kgf/cm2")
                    BFW_01.SetTemperature(BFW_T)
                    BFW_01.SetMassFlow(f"{self.BFW_FLOW[x-1] * 1000} kg/h")
                    HE_01.set_ColdSideOutletTemperature(self.STM_TEMP[x-1]+ 273.15)

                    LHV = self.CAL_LHV(x) * 4.1868 / 3600
                    set_q = LHV + (Q1 * 4.1868 / 3600)      # [kW]
                    FURN_01.set_DeltaQ(set_q)
                    CRACK_01.set_OutletTemperature(self.STACK_T[x-1] + 273)
                    interf.CalculateFlowsheet2(flowsheet)
                    
                    R_FLUE_T2 = STACK_01.GetTemperature() - 273.15
                    R_FLUE_O2 = STACK_01.GetOverallComposition()[2] * 100

                    LOSS = LOSS_01.get_DeltaQ() - STACK_01.GetOverallComposition()[3] * STACK_01.GetMolarFlow() * 3.6 * 18.0152 * 584
                    RAD_LOSS = 0.025
                    EFF_HL = 1 - (LOSS/LHV + RAD_LOSS)

                    arr = [R_FLUE_T2, R_FLUE_O2, LHV, LOSS, EFF_HL, self.OK, self.HAIR_TOT]
                    tup_arr = tuple(arr)
                    sql = "INSERT INTO FURN_6_2 VALUES (?,?,?,?,?,?,?)"
                    cur.execute(sql, tup_arr)
                    conn.commit()
                    print(str(i+1) + " -> success!" )
                    interf.SaveFlowsheet(flowsheet, new_path + "dump.dwxmz", True)
                else :
                    arr = [NULL, NULL, NULL, NULL, NULL, self.OK]
                    tup_arr = tuple(arr)
                    sql = "INSERT INTO FURN_6_2 VALUES (?,?,?,?,?,?)"
                    cur.execute(sql, tup_arr)
                    conn.commit()
                    print(str(i+1) + " -> fail!" )

#################################################################################
#################################################################################
class BLR:

    def callback(self, xk):
        if (self.obj_1 < 0.0001) and (self.obj_2 < 0.003):
            print("Terminating optimization")
            raise RuntimeError("stop!")
        else:
            return False
    
    def callback2(self, xk):
        if (self.obj_0 < 1E-07):
            print("Terminating optimization")
            raise RuntimeError("stop!")
        else:
            return False

    def BLR_DWSIM(self):
        num = 1



    def opt_o2(self, x, *args):
        init0 = args[0]
        # AIR_R.SetMassFlow(f"{x[0]} kg/h")
        
        # interf.CalculateFlowsheet2(flowsheet2)

        # obj_a = float(init0)
        # obj_b = FLUE_GAS_OUT_R.GetOverallComposition()[9]
        # self.obj_0 = (obj_b - obj_a)**2
        # print("Air Flow : " + str(round(x[0],4)))
        # print("obj : " + str(round(obj_a * 100,5)) + "<->" + str(round(obj_b * 100,5)))

        return self.obj_0

    def opt_gas(self, x, *args):
        init1 = args[0]
        init2 = args[1]
        # GAS_IN.SetVolumetricFlow(f"{x[0]} m3/h")
        # AIR_IN.SetMassFlow(f"{x[1]} kg/h")

        # interf.CalculateFlowsheet2(flowsheet)

        # print(datetime.datetime.today())
        # obj_a = float(init1)
        # obj_b = FLUE_GAS_OUT2.GetTemperature() - 273.15
        # obj_c = float(init2) 
        # obj_d = FLUE_GAS_OUT2.GetOverallComposition()[9] 
        # self.obj_1 = ((obj_b - obj_a)/10)**2
        # self.obj_2 = (((obj_d - obj_c)/0.01)**2)
        # self.obj = self.obj_1 + self.obj_2
        # print(x[0], x[1])
        # print("temp : "+str(round(init1,2))+"<->"+str(round(FLUE_GAS_OUT2.GetTemperature() - 273.15,2)))
        # print("o2 : "+str(round(obj_c,4))+"<->"+str(round(obj_d,4)))
        # print("OBJ1 : " + str(self.obj_1))
        # print("OBJ2 : " + str(self.obj_2))
        # print("obj : " + str(self.obj))

        return self.obj

    def run(self):
        self.BLR_DWSIM()  



if __name__ == "__main__":
    CAL_GAS().SOLVE(6)




