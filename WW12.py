import pythoncom
pythoncom.CoInitialize()

import os
import clr
import sys
import sqlite3
import pandas as pd
import numpy as np
import json
import pandas as pd
from scipy.optimize import minimize
from System.IO import Directory
from System import Array


# Set & Load DWSIM Flowsheet
mod = sys.modules[__name__]

dwsim_install_path = "C:/Users/lovpa/AppData/Local/DWSIM8"
sys.path.append(dwsim_install_path)
sys.path.append(os.getcwd() + "\\")
new_path = os.getcwd() + "\\"

clr.AddReference(new_path + "System.Buffers.dll")
clr.AddReference(dwsim_install_path + "/DWSIM.Automation.dll")
clr.AddReference(dwsim_install_path + "/DWSIM.Interfaces.dll")
clr.AddReference('DWSIM')

from DWSIM.Automation import Automation, Automation2

Directory.SetCurrentDirectory(dwsim_install_path)
interf = Automation2()

flowsheet = interf.LoadFlowsheet(new_path + "main_furn.dwxmz")
flowsheet2 = interf.LoadFlowsheet(new_path + "fuel_ch.dwxmz")

Fuel_in = flowsheet2.GetFlowsheetSimulationObject("Fuel_in")
Fuel_out = flowsheet2.GetFlowsheetSimulationObject("Fuel_out")
Fuel_Heat = flowsheet2.GetFlowsheetSimulationObject("Fuel_Heat")

Inlet_1 = flowsheet.GetFlowsheetSimulationObject("INLET_1")
Inlet_2 = flowsheet.GetFlowsheetSimulationObject("INLET_2")
Inlet = flowsheet.GetFlowsheetSimulationObject("INLET")

STR_01 = flowsheet.GetFlowsheetSimulationObject("STR_01")
FURN_01 = flowsheet.GetFlowsheetSimulationObject("FURN_01")
BFW_01 = flowsheet.GetFlowsheetSimulationObject("BFW_01")
HE_01 = flowsheet.GetFlowsheetSimulationObject("HE_01")
CRACK_01 = flowsheet.GetFlowsheetSimulationObject("CRACK_01")
LOSS_01 = flowsheet.GetFlowsheetSimulationObject("LOSS_01")
STACK_01 = flowsheet.GetFlowsheetSimulationObject("STACK_01")
STACK_OUT_01 = flowsheet.GetFlowsheetSimulationObject("STACK_OUT_01")
#========================================================================
conn = sqlite3.connect(new_path + 'Blank_db.db')
cur = conn.cursor()


class CAL_FURN_EFF():
    def reduce_mem_usage(self, df):
        start_mem = df.memory_usage().sum() / 1024**2
        print('Memory usage of dataframe is {:.2f} MB'.format(start_mem))
        
        for col in df.columns:
            col_type = df[col].dtype

            if col_type != object:
                c_min = df[col].min()
                c_max = df[col].max()
                if str(col_type)[:3] == 'int':
                    if c_min > np.iinfo(np.int8).min and c_max < np.iinfo(np.int8).max:
                        df[col] = df[col].astype(np.int8)
                    elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                        df[col] = df[col].astype(np.int16)
                    elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                        df[col] = df[col].astype(np.int32)
                    elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                        df[col] = df[col].astype(np.int64)

                elif str(col_type)[:4] == 'date':   
                    pass

                else:
                    if c_min > np.finfo(np.float16).min and c_max < np.finfo(np.float16).max:
                        df[col] = df[col].astype(np.float16)
                    elif c_min > np.finfo(np.float32).min and c_max < np.finfo(np.float32).max:
                        df[col] = df[col].astype(np.float32)
                    else:
                        df[col] = df[col].astype(np.float64)
            else:
                pass

        end_mem = df.memory_usage().sum() / 1024**2
        print('Memory usage after optimization is: {:.2f} MB'.format(end_mem))
        print('Decreased by {:.1f}%'.format(100 * (start_mem - end_mem) / start_mem))
        
        return df

    def LOAD_DATA(self):
        df = pd.read_excel(new_path + 'data2.xlsx')
        self.data = self.reduce_mem_usage(df)
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
        self.arr_HHV = arr[7].astype(np.float16)
        self.arr_length = arr.shape[1]

    def LOAD_FURN_DATA (self, k):
        arr_COMP1 = np.array([22.0841, 67.8167, 5.1084, 1.0459, 1.4997, 0, 0.0345, 0.0138])                         # COMP. [H2, CH4, C2H6, C2H4, C3H8, △-C3H6, C3H6, C2H2]
        arr_COMP2 = np.array([0.0408, 0, 0.1303, 0.1897, 0.1044, 0.1532, 0.4364])                                   # COMP. [i-C4*, PD, n-C4*, T-2-C4-, 1-C4-, i-C4-, C-2-C4-]
        arr_COMP3 = np.array([0.0003, 0.0613, 0.0679, 0.0010, 0.2205, 0.2606, 0.0203, 0.1508])                      # COMP. [22-DM-PROP, i-C5*, 12-BUTAD, MA, 13-BUTAD, 3M-1-C4-, VA, EA]
        arr_COMP4 = np.array([0.1185, 0.0022, 0.6421, 0.0768])                                                      # COMP. [CO, CO2, N2, n-C5*]
        arr_COMP5 = np.array([0.0084, 0.0240, 0.0003, 0.0012, 0.0001, 0.0001])                                      # COMP. [HEXANE, BENZENE, HEPTANE, TOLUENE, EB, XYLENE]

        self.arr_comp = np.block([arr_COMP1,arr_COMP2,arr_COMP3,arr_COMP4,arr_COMP5])

        var_t0 = self.data['F1_INLET_T'][k]
        var_t = self.data['F8_INLET_T'][k] - self.data['F1_INLET_T'][k]
        if self.flag == 1:
            self.Fuel_Gas_Mass = self.data['F1_Fuel_Gas_M'][k]
            self.Fuel_Gas_Temp = self.data['F1_Fuel_Gas_T'][k]
            self.Air_Mass = self.data['F1_Air_M'][k]
            self.Air_CO2 = self.data['F1_Air_CO2'][k]
            self.Air_H2O = self.data['F1_Air_H2O'][k]
            self.Air_N2 = self.data['F1_Air_N2'][k]
            self.Air_O2 = self.data['F1_Air_O2'][k]
            self.Flue_T = self.data['F1_Flue_T'][k]
            self.BFW_P = self.data['F1_BFW_PRESS'][k]
            self.BFW_T = self.data['F1_BFW_TEMP'][k]
            self.BFW_F = self.data['F1_BFW_FLOW'][k]
            self.STM_T = self.data['F1_STEAM_TEMP'][k]
            self.Steam_P = self.data['F1_BFW_PRESS'][k] - self.data['F1_STEAM_PRESS'][k]
            self.init = self.data['F1_EXO2'][k]
            self.INLET_TEMP = var_t0 + var_t*0/7
        
        elif self.flag == 2:
            self.Fuel_Gas_Mass = self.data['F2_Fuel_Gas_M'][k]
            self.Fuel_Gas_Temp = self.data['F2_Fuel_Gas_T'][k]
            self.Air_Mass = self.data['F2_Air_M'][k]
            self.Air_CO2 = self.data['F2_Air_CO2'][k]
            self.Air_H2O = self.data['F2_Air_H2O'][k]
            self.Air_N2 = self.data['F2_Air_N2'][k]
            self.Air_O2 = self.data['F2_Air_O2'][k]
            self.Flue_T = self.data['F2_Flue_T'][k]
            self.BFW_P = self.data['F2_BFW_PRESS'][k]
            self.BFW_T = self.data['F2_BFW_TEMP'][k]
            self.BFW_F = self.data['F2_BFW_FLOW'][k]
            self.STM_T = self.data['F2_STEAM_TEMP'][k]
            self.Steam_P = self.data['F2_BFW_PRESS'][k] - self.data['F2_STEAM_PRESS'][k]
            self.init = self.data['F2_EXO2'][k]
            self.INLET_TEMP = var_t0 + var_t*1/7

        elif self.flag == 3:
            self.Fuel_Gas_Mass = self.data['F3_Fuel_Gas_M'][k]
            self.Fuel_Gas_Temp = self.data['F3_Fuel_Gas_T'][k]
            self.Air_Mass = self.data['F3_Air_M'][k]
            self.Air_CO2 = self.data['F3_Air_CO2'][k]
            self.Air_H2O = self.data['F3_Air_H2O'][k]
            self.Air_N2 = self.data['F3_Air_N2'][k]
            self.Air_O2 = self.data['F3_Air_O2'][k]
            self.Flue_T = self.data['F3_Flue_T'][k]
            self.BFW_P = self.data['F3_BFW_PRESS'][k]
            self.BFW_T = self.data['F3_BFW_TEMP'][k]
            self.BFW_F = self.data['F3_BFW_FLOW'][k]
            self.STM_T = self.data['F3_STEAM_TEMP'][k]
            self.Steam_P = self.data['F3_BFW_PRESS'][k] - self.data['F3_STEAM_PRESS'][k]
            self.init = self.data['F3_EXO2'][k]
            self.INLET_TEMP = var_t0 + var_t*2/7

        elif self.flag == 4:
            self.Fuel_Gas_Mass = self.data['F4_Fuel_Gas_M'][k]
            self.Fuel_Gas_Temp = self.data['F4_Fuel_Gas_T'][k]
            self.Air_Mass = self.data['F4_Air_M'][k]
            self.Air_CO2 = self.data['F4_Air_CO2'][k]
            self.Air_H2O = self.data['F4_Air_H2O'][k]
            self.Air_N2 = self.data['F4_Air_N2'][k]
            self.Air_O2 = self.data['F4_Air_O2'][k]
            self.Flue_T = self.data['F4_Flue_T'][k]
            self.BFW_P = self.data['F4_BFW_PRESS'][k]
            self.BFW_T = self.data['F4_BFW_TEMP'][k]
            self.BFW_F = self.data['F4_BFW_FLOW'][k]
            self.STM_T = self.data['F4_STEAM_TEMP'][k]
            self.Steam_P = self.data['F4_BFW_PRESS'][k] - self.data['F4_STEAM_PRESS'][k]
            self.init = self.data['F4_EXO2'][k]
            self.INLET_TEMP = var_t0 + var_t*3/7

        elif self.flag == 5:
            self.Fuel_Gas_Mass = self.data['F5_Fuel_Gas_M'][k]
            self.Fuel_Gas_Temp = self.data['F5_Fuel_Gas_T'][k]
            self.Air_Mass = self.data['F5_Air_M'][k]
            self.Air_CO2 = self.data['F5_Air_CO2'][k]
            self.Air_H2O = self.data['F5_Air_H2O'][k]
            self.Air_N2 = self.data['F5_Air_N2'][k]
            self.Air_O2 = self.data['F5_Air_O2'][k]
            self.Flue_T = self.data['F5_Flue_T'][k]
            self.BFW_P = self.data['F5_BFW_PRESS'][k]
            self.BFW_T = self.data['F5_BFW_TEMP'][k]
            self.BFW_F = self.data['F5_BFW_FLOW'][k]
            self.STM_T = self.data['F5_STEAM_TEMP'][k]
            self.Steam_P = self.data['F5_BFW_PRESS'][k] - self.data['F5_STEAM_PRESS'][k]
            self.init = self.data['F5_EXO2'][k]
            self.INLET_TEMP = var_t0 + var_t*4/7

        elif self.flag == 6:
            self.Fuel_Gas_Mass = self.data['F6_Fuel_Gas_M'][k]
            self.Fuel_Gas_Temp = self.data['F6_Fuel_Gas_T'][k]
            self.Air_Mass = self.data['F6_Air_M'][k]
            self.Air_CO2 = self.data['F6_Air_CO2'][k]
            self.Air_H2O = self.data['F6_Air_H2O'][k]
            self.Air_N2 = self.data['F6_Air_N2'][k]
            self.Air_O2 = self.data['F6_Air_O2'][k]
            self.Flue_T = self.data['F6_Flue_T'][k]
            self.BFW_P = self.data['F6_BFW_PRESS'][k]
            self.BFW_T = self.data['F6_BFW_TEMP'][k]
            self.BFW_F = self.data['F6_BFW_FLOW'][k]
            self.STM_T = self.data['F6_STEAM_TEMP'][k]
            self.Steam_P = self.data['F6_BFW_PRESS'][k] - self.data['F6_STEAM_PRESS'][k]
            self.init = self.data['F6_EXO2'][k]
            self.INLET_TEMP = var_t0 + var_t*5/7

        elif self.flag == 7:
            self.Fuel_Gas_Mass = self.data['F7_Fuel_Gas_M'][k]
            self.Fuel_Gas_Temp = self.data['F7_Fuel_Gas_T'][k]
            self.Air_Mass = self.data['F7_Air_M'][k]
            self.Air_CO2 = self.data['F7_Air_CO2'][k]
            self.Air_H2O = self.data['F7_Air_H2O'][k]
            self.Air_N2 = self.data['F7_Air_N2'][k]
            self.Air_O2 = self.data['F7_Air_O2'][k]
            self.Flue_T = self.data['F7_Flue_T'][k]
            self.BFW_P = self.data['F7_BFW_PRESS'][k]
            self.BFW_T = self.data['F7_BFW_TEMP'][k]
            self.BFW_F = self.data['F7_BFW_FLOW'][k]
            self.STM_T = self.data['F7_STEAM_TEMP'][k]
            self.Steam_P = self.data['F7_BFW_PRESS'][k] - self.data['F7_STEAM_PRESS'][k]
            self.INLET_TEMP = var_t0 + var_t*6/7

        elif self.flag == 8:
            self.Fuel_Gas_Mass = self.data['F8_Fuel_Gas_M'][k]
            self.Fuel_Gas_Temp = self.data['F8_Fuel_Gas_T'][k]
            self.Air_Mass = self.data['F8_Air_M'][k]
            self.Air_CO2 = self.data['F8_Air_CO2'][k]
            self.Air_H2O = self.data['F8_Air_H2O'][k]
            self.Air_N2 = self.data['F8_Air_N2'][k]
            self.Air_O2 = self.data['F8_Air_O2'][k]
            self.Flue_T = self.data['F8_Flue_T'][k]
            self.BFW_P = self.data['F8_BFW_PRESS'][k]
            self.BFW_T = self.data['F8_BFW_TEMP'][k]
            self.BFW_F = self.data['F8_BFW_FLOW'][k]
            self.STM_T = self.data['F8_STEAM_TEMP'][k]
            self.Steam_P = self.data['F8_BFW_PRESS'][k] - self.data['F8_STEAM_PRESS'][k]
            self.INLET_TEMP = var_t0 + var_t*7/7   

    def CAL_GEN(self):
        FUEL_MWT = (np.dot(self.arr_comp, self.arr_MWT))/100
        FUEL_MOLE = self.Fuel_Gas_Mass / FUEL_MWT
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

    def CAL_LHV(self):
        FUEL_MWT = (np.dot(self.arr_comp, self.arr_MWT))/100
        LHV = 0
        for i in range(self.arr_length):
            LHV_ = self.arr_comp[i] * self.arr_LHV[i] * self.arr_MWT[i]
            LHV = LHV + LHV_
        FUEL_LHV = 4.184 * LHV / (100* FUEL_MWT)
        
        return FUEL_LHV

    def CAL_HHV(self):
        FUEL_MWT = (np.dot(self.arr_comp, self.arr_MWT))/100
        HHV = 0
        for i in range(self.arr_length):
            HHV_ = self.arr_comp[i] * self.arr_HHV[i] * self.arr_MWT[i]
            HHV = HHV + HHV_
        FUEL_HHV = 4.184 * HHV / (100* FUEL_MWT)
        
        return FUEL_HHV      

    def CAL_FURN(self):
        LHV = self.CAL_LHV()
        HHV = self.CAL_HHV()
        self.CAL_GEN()
        Fuel_Heat.set_OutletTemperature(self.Fuel_Gas_Temp + 273.15)
        Fuel_in.SetMassFlow(f"{self.Fuel_Gas_Mass} kg/h")
        interf.CalculateFlowsheet2(flowsheet2)
        
        Q1 = Fuel_Heat.get_DeltaQ()                                             # [kW]
        FUEL_LHV = self.Fuel_Gas_Mass * LHV/3600                                # [kW]
        FUEL_HHV = self.Fuel_Gas_Mass * HHV/3600                                # [kW]

        Flue_gas_mass = self.Fuel_Gas_Mass + self.Air_Mass
        air_low_bound = 0.4 * Flue_gas_mass
        air_high_bound = 100.0 * Flue_gas_mass
        bnds = [(air_low_bound, air_high_bound)]
        Inlet.set_OutletTemperature(self.INLET_TEMP + 273.15)
        print("=========" + str(self.flag) + "호기 최적화=========" )
        res = minimize(self.OPT_O2, self.Air_Mass, args=(self.init), method='Nelder-Mead', bounds=bnds, options={'disp':True, 'maxiter': 100, 'fatol' : 0.05, 'xatol' : 0.1})
        air_mass = Inlet_1.GetMassFlow() * 3600

        Q2_TOT = Inlet.get_DeltaQ()
        try:
            Q2_vap_heat = Inlet_1.Phases[3].Properties.get_massflow() * 2443.456
        except:
            Q2_vap_heat = 0
        Q2 = Q2_TOT - Q2_vap_heat

        Flue_CO2 = self.x1 / (self.x1 + self.x2 + self.x3 + self.x4)
        Flue_H2O = self.x2 / (self.x1 + self.x2 + self.x3 + self.x4)
        Flue_N2 = self.x3 / (self.x1 + self.x2 + self.x3 + self.x4)
        Flue_O2 = self.x4 / (self.x1 + self.x2 + self.x3 + self.x4)

        STR_01.SetOverallComposition(Array[float]((Flue_CO2, Flue_N2, Flue_O2, Flue_H2O)))
        Flue_gas_mass = self.Fuel_Gas_Mass + air_mass
        STR_01.SetMassFlow(f"{Flue_gas_mass} kg/h")
        Input_heat_HHV = (FUEL_HHV + Q1 + Q2_TOT)                               # [kW]
        Input_heat_LHV = (FUEL_LHV + Q1 + Q2)                                   # [kW]
        LOSS_RAD_HHV = Input_heat_HHV * 0.01
        LOSS_RAD_LHV = Input_heat_LHV * 0.01
        FURN_01.set_DeltaQ(Input_heat_HHV)
        BFW_01.SetPressure(f"{self.BFW_P} kgf/cm2")
        # BFW_01.SetPressure(f"{5} kgf/cm2")
        BFW_01.SetTemperature(self.BFW_T + 273.15)
        BFW_01.SetMassFlow(f"{self.BFW_F * 1000} kg/h")
        HE_01.set_ColdSideOutletTemperature(self.STM_T + 273.15)
        HE_01.set_ColdSidePressureDrop(self.Steam_P*98070*2)                    # [Pa] ##??
        CRACK_01.set_OutletTemperature(self.Flue_T + 273.15)
        interf.CalculateFlowsheet2(flowsheet)

        Steam_heat = HE_01.get_MaxHeatExchange() * (HE_01.get_ThermalEfficiency()/100) 
        Crack_heat_HHV = CRACK_01.get_DeltaQ()
        Crack_heat_HHV = Crack_heat_HHV - LOSS_RAD_HHV                    

        Q3_TOT = LOSS_01.get_DeltaQ()
        Q3_vap_heat = STACK_OUT_01.Phases[3].Properties.get_massflow() * 2443.456 
        Q3 = Q3_TOT - Q3_vap_heat
        Crack_heat_LHV = Input_heat_LHV - LOSS_RAD_LHV - Steam_heat - Q3
        Eff_HHV = 1 - Q3_TOT/Input_heat_HHV - 0.01
        Eff_LHV = 1 - (Q3 / Input_heat_LHV) - 0.01
        print("=============================")
        print("Complete "+ str(self.flag) +"호기 Furnace" )
        print("INPUT HEAT : " + str(Input_heat_HHV))
        print("STEAM HEAT : " + str(Steam_heat))
        print("CRACK HEAT : " + str(Crack_heat_HHV))
        print("STACK LOSS : " + str(Q3))
        print("효율(HHV) : " + str(Eff_HHV))
        print("효율(LHV) : " + str(Eff_LHV))
        print("=============================")
        interf.SaveFlowsheet(flowsheet, new_path + "dump.dwxmz", True)

        stack_t = STACK_01.GetTemperature() - 273.15                        # ['C]
        stack_o2 = STACK_OUT_01.GetOverallComposition()[2] * 100            # [(%)]

        arr = [air_mass, stack_o2, stack_t, FUEL_HHV, FUEL_LHV ,Q1, Q2, Q2_vap_heat, Input_heat_HHV, Input_heat_LHV, Crack_heat_HHV, Crack_heat_LHV, Steam_heat, Q3_vap_heat, Q3, Eff_HHV, Eff_LHV]
        self.tup_arr = tuple(arr)

        
    def OPT_O2(self, x, *args):
        init0 = args[0]
        Inlet_1.SetOverallComposition(Array[float]((self.Air_CO2, self.Air_N2, self.Air_O2, self.Air_H2O)))
        Inlet_1.SetMassFlow(f"{x[0]} kg/h")
        interf.CalculateFlowsheet2(flowsheet)

        inlet_mf = Inlet_1.GetMolarFlow() * 3.6
        CO2_mass = inlet_mf * self.Air_CO2 * 44.01
        H2O_mass = inlet_mf * self.Air_H2O * 18.015
        N2_mass = inlet_mf * self.Air_N2 * 28.013
        O2_mass = inlet_mf * self.Air_O2 * 31.999

        CO2_frac = CO2_mass / (CO2_mass + H2O_mass + N2_mass + O2_mass)
        H2O_frac = H2O_mass / (CO2_mass + H2O_mass + N2_mass + O2_mass)
        N2_frac = N2_mass / (CO2_mass + H2O_mass + N2_mass + O2_mass)
        O2_frac = O2_mass / (CO2_mass + H2O_mass + N2_mass + O2_mass)

        self.x1 = (x[0] * CO2_frac + self.GEN_CO2 * 44.01) / 44.01
        self.x2 = (x[0] * H2O_frac + self.GEN_H2O * 18.015) / 18.015
        self.x3 = (x[0] * N2_frac + self.GEN_N2 * 28.013) / 28.013
        self.x4 = (x[0] * O2_frac - self.GEN_O2 * 31.999) / 31.999
        
        Flue_O2 = self.x4 / (self.x1 + self.x2 + self.x3 + self.x4)
        print(f"Iteration --> O2% = {Flue_O2*100} // Hot_air = {x[0]}")
        self.obj_0 = (Flue_O2*100 - init0)**2

        return self.obj_0  

    def RUN(self):
        self.LOAD_DATA()
        for i in range(0, 1):
            print("Processing .. Data No. " + str(i))
            # # 1호기 furnace 
            # self.flag = 1                              
            # self.LOAD_FURN_DATA(i)
            # self.CAL_FURN()
            # sql = "INSERT INTO FURN_1 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
            # cur.execute(sql, self.tup_arr)
            # conn.commit()
            # # 2호기 furnace 
            # self.flag = 2                              
            # self.LOAD_FURN_DATA(i)
            # self.CAL_FURN()
            # sql = "INSERT INTO FURN_2 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
            # cur.execute(sql, self.tup_arr)
            # conn.commit()
            # # 3호기 furnace 
            # self.flag = 3                              
            # self.LOAD_FURN_DATA(i)
            # self.CAL_FURN()
            # sql = "INSERT INTO FURN_3 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
            # cur.execute(sql, self.tup_arr)
            # conn.commit()
            # # 4호기 furnace 
            # self.flag = 4                              
            # self.LOAD_FURN_DATA(i)
            # self.CAL_FURN()
            # sql = "INSERT INTO FURN_4 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
            # cur.execute(sql, self.tup_arr)
            # conn.commit()
            # # 5호기 furnace 
            # self.flag = 5                            
            # self.LOAD_FURN_DATA(i)
            # self.CAL_FURN()
            # sql = "INSERT INTO FURN_5 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
            # cur.execute(sql, self.tup_arr)
            # conn.commit()
            # 6호기 furnace 
            self.flag = 6                               
            self.LOAD_FURN_DATA(i)
            self.CAL_FURN()
            sql = "INSERT INTO FURN_6 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
            cur.execute(sql, self.tup_arr)
            conn.commit()
            # # 7호기 furnace 
            # self.flag = 7                               
            # self.LOAD_FURN_DATA(i)
            # self.CAL_FURN()
            # sql = "INSERT INTO FURN_7 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
            # cur.execute(sql, self.tup_arr)
            # conn.commit()
            # # 8호기 furnace 
            # self.flag = 8                               
            # self.LOAD_FURN_DATA(i)
            # self.CAL_FURN()
            # sql = "INSERT INTO FURN_6 VALUES (?,?,?,?,?,?,?,?,?,?,?,?)"
            # cur.execute(sql, self.tup_arr)
            # conn.commit()

    def debug(self):
        for i in range(100, 1300):
            self.MAKE_ARR(i)
            self.CAL_GEN()
            Inlet_1.SetOverallComposition(Array[float]((self.Air_CO2, self.Air_N2, self.Air_O2, self.Air_H2O)))
            Inlet_1.SetMassFlow(f"{150000} kg/h")
            interf.CalculateFlowsheet2(flowsheet)

            inlet_mf = Inlet_1.GetMolarFlow() * 3.6
            CO2_mass = inlet_mf * self.Air_CO2 * 44.01
            H2O_mass = inlet_mf * self.Air_H2O * 18.015
            N2_mass = inlet_mf * self.Air_N2 * 28.013
            O2_mass = inlet_mf * self.Air_O2 * 31.999

            CO2_frac = CO2_mass / (CO2_mass + H2O_mass + N2_mass + O2_mass)
            H2O_frac = H2O_mass / (CO2_mass + H2O_mass + N2_mass + O2_mass)
            N2_frac = N2_mass / (CO2_mass + H2O_mass + N2_mass + O2_mass)
            O2_frac = O2_mass / (CO2_mass + H2O_mass + N2_mass + O2_mass)

            self.x1 = (150000 * CO2_frac + self.GEN_CO2 * 44.01) / 44.01
            self.x2 = (150000 * H2O_frac + self.GEN_H2O * 18.015) / 18.015
            self.x3 = (150000 * N2_frac + self.GEN_N2 * 28.013) / 28.013
            self.x4 = (150000 * O2_frac - self.GEN_O2 * 31.999) / 31.999
            
            Flue_O2 = self.x4 / (self.x1 + self.x2 + self.x3 + self.x4)
            print(Flue_O2)
#################################################################################
#################################################################################


if __name__ == "__main__":
    CAL_FURN_EFF().RUN()
