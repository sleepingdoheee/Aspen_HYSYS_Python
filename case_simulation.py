import pandas as pd
import numpy as np
import json
import sqlite3


# from utils.insert_data import _make_insert_sql

class CASE_SIMULATION:
    def load_data(self):
        """시뮬레이션을 위한 필요 데이터 불러오기"""
        # Connect TechDas DB(sqlite)
        # conn = sqlite3.connect("") ## db 이름 필요
        # c = conn.cursor()
        # # Load Input Data from Table
        # res = c.execute("select fuel1, fuel2, fuel3, fuel4, fuel5, fuel6, fuel7, fuel8 from case_simul_input;")
        # fgas_flow = list(res.fetchall()[0])

        # res = c.execute("select air1, air2, air3, air4, air5, air6, air7, air8 from case_simul_input;")
        # air_flow = list(res.fetchall()[0])
        
        # res = c.execute("select f_bias1, f_bias2, f_bias3, f_bias4, f_bias5, f_bias6, f_bias7, f_bias8 from case_simul_input;")
        # f_bias = list(res.fetchall()[0])
        
        # res = c.execute("select a_bias1, a_bias2, a_bias3, a_bias4, a_bias5, a_bias6, a_bias7, a_bias8 from case_simul_input;")
        # a_bias = list(res.fetchall()[0])

        # res = c.execute("select gtg_fuel, gtg_flue_o2 from case_simul_input;")
        # gtg_data = list(res.fetchall()[0])
        
        # # Load Values for gtg_composition
        # res = c.execute("select * from gtg_composition;")
        # gtg_composition = res.fetchall()
        
        # # Load Values for decomp_composition
        # res = c.execute("select * from decomp_composition;")
        # decomposer_composition = res.fetchall()
        # conn.close()

        # self.furn_fgas_flow = [round(i,3) for i in fgas_flow]
        # self.air_list = [round(i,3) for i in air_flow]
        # self.fuel_bias_calc = [round(i,3) for i in f_bias]
        # self.air_bias_calc = [round(i,3) for i in a_bias]
        # self.gtg_data = [round(i,3) for i in gtg_data]

        # self.gtgfuel_comp = [float(gtg_composition[i][-1]) for i in range(len(gtg_composition))] 
        # self.fuel_comp = [float(decomposer_composition[i][-1]) for i in range(len(decomposer_composition))]

        arr_COMP1 = np.array([22.0841, 67.8167, 5.1084, 1.0459, 1.4997, 0, 0.0345, 0.0138])                         # COMP. [H2, CH4, C2H6, C2H4, C3H8, △-C3H6, C3H6, C2H2]
        arr_COMP2 = np.array([0.0408, 0, 0.1303, 0.1897, 0.1044, 0.1532, 0.4364])                                   # COMP. [i-C4*, PD, n-C4*, T-2-C4-, 1-C4-, i-C4-, C-2-C4-]
        arr_COMP3 = np.array([0.0003, 0.0613, 0.0679, 0.0010, 0.2205, 0.2606, 0.0203, 0.1508])                      # COMP. [22-DM-PROP, i-C5*, 12-BUTAD, MA, 13-BUTAD, 3M-1-C4-, VA, EA]
        arr_COMP4 = np.array([0.1185, 0.0022, 0.6421, 0.0768])                                                      # COMP. [CO, CO2, N2, n-C5*]
        arr_COMP5 = np.array([0.0084, 0.0240, 0.0003, 0.0012, 0.0001, 0.0001])                                      # COMP. [HEXANE, BENZENE, HEPTANE, TOLUENE, EB, XYLENE]
        self.fuel_comp = np.block([arr_COMP1,arr_COMP2,arr_COMP3,arr_COMP4,arr_COMP5])

        # Load Fuel Gas Coefficient
        with open('gas_comp.json', 'r', encoding="utf-8") as f:
            sample = json.load(f)
        arr = np.array(sample)
        arr = arr.flatten('F').reshape(arr.shape[1], arr.shape[0])
        self.arr_name = arr[0]
        self.arr_MWT = arr[1].astype(np.float16)
        self.arr_R1 = arr[2].astype(np.float16)
        self.arr_R2 = arr[3].astype(np.float16)
        self.arr_R3 = arr[4].astype(np.float16)
        self.arr_R4 = arr[5].astype(np.float16)
        self.arr_length = arr.shape[1]

    def gtg_calc(self):
        ''' GTG 배기가스 유량 및 조성 계산 '''
        # INPUT-------------
        # self.gtg_fgas_flow = self.gtg_data[0] * 0.9                     # kg/hr, GTG 연료 유량, NC1FIT600
        # self.gtg_stack_o2 = self.gtg_data[1]                            # mol%, GTG 배기가스 O2%, NC1AI1010
        # GTG_comp = self.gtgfuel_comp                                    # COMP. [H2, CO, CH4, C2H6, C2H4, C2H2]
        self.gtg_fgas_flow = 1 * 0.9
        self.gtg_stack_o2 = 20.6
        GTG_comp = [24.2622, 0.1219, 75.4055, 0.0011, 0.1791, 0.0197]
        #-------------------
        GTG_Fuel_Flow = self.gtg_fgas_flow  # kg/hr
        GTG_Stack_O2 = self.gtg_stack_o2  # mol%

        ## GTG Fuel 몰질량 계산(kg/kmol)
        GTG_Fuel_MW = (GTG_comp[0] * 2.0160 + GTG_comp[1] * 28.01 + GTG_comp[2] * 16.0429 + GTG_comp[3] * 30.0699 + GTG_comp[4] * 28.0538 + GTG_comp[5] * 26.0380) / 100 # kg/kmol

        ## GTG 배기가스 유량 계산(m3/hr)
        GTG_O2 = GTG_Fuel_Flow / GTG_Fuel_MW * (GTG_comp[0] / 2 + GTG_comp[1] / 2 + GTG_comp[2] * 2 + GTG_comp[3] * 3.5 + GTG_comp[4] * 3 + GTG_comp[5] * 2.5) / 100 * 22.4
        GTG_CO2 = GTG_Fuel_Flow / GTG_Fuel_MW * (GTG_comp[0] * 0 + GTG_comp[1] * 1 + GTG_comp[2] * 1 + GTG_comp[3] * 2 + GTG_comp[4] * 2 + GTG_comp[5] * 2) / 100 * 22.4
        GTG_H2O = GTG_Fuel_Flow / GTG_Fuel_MW * (GTG_comp[0] * 1 + GTG_comp[1] * 0 + GTG_comp[2] * 2 + GTG_comp[3] * 3 + GTG_comp[4] * 2 + GTG_comp[5] * 1) / 100 * 22.4
        GTG_N2 = GTG_O2 / 0.21 * 0.79
        GTG_N2_2 = (GTG_Stack_O2 / 100 * (GTG_H2O + GTG_CO2 + GTG_N2)) / (0.21 - GTG_Stack_O2 / 100) * 0.79 + GTG_N2
        GTG_O2_2 = (GTG_Stack_O2 / 100 * (GTG_N2 + GTG_CO2 + GTG_H2O)) / (0.21 - GTG_Stack_O2 / 100) * 0.21
        GTG_EXHAUST_AIR = GTG_O2_2 + GTG_N2_2 + GTG_CO2 + GTG_H2O

        ## GTG 배기가스 조성 계산(mol%)
        self.GTG_O2_Frac = GTG_O2_2 / GTG_EXHAUST_AIR
        self.GTG_N2_Frac = GTG_N2_2 / GTG_EXHAUST_AIR
        self.GTG_H2O_Frac = GTG_H2O / GTG_EXHAUST_AIR
        self.GTG_CO2_Frac = GTG_CO2 / GTG_EXHAUST_AIR

        ## GTG 배기가스 유량 계산(kg/hr)
        GTG_CO2_W = GTG_CO2 * 44 / 22.4
        GTG_H2O_W = GTG_H2O * 18 / 22.4
        GTG_N2_W = GTG_N2_2 * 28 / 22.4
        GTG_O2_W = GTG_O2_2 * 32 / 22.4
        self.EX_GTG_W = GTG_CO2_W + GTG_H2O_W + GTG_N2_W + GTG_O2_W
        self.GTG_O2_W = GTG_O2_W


    def hotair_calc(self, a):
        ''' HOT AIR 조성 계산'''
        Air_Ratio = a
        self.HOTAIR_O2_Frac = self.GTG_O2_Frac * Air_Ratio + 0.21 * (1 - Air_Ratio)
        self.HOTAIR_N2_Frac = self.GTG_N2_Frac * Air_Ratio + 0.79 * (1 - Air_Ratio)
        self.HOTAIR_CO2_Frac = self.GTG_CO2_Frac * Air_Ratio
        self.HOTAIR_H2O_Frac = self.GTG_H2O_Frac * Air_Ratio


    def furn_calc(self, a, b):
        ''' 분해로 공기 유량 계산 '''
        FURN_Fuel_Flow = self.furn_fgas_flow[a]     # kg/hr, Input
        FURN_EXCESS_O2 = float(b)                   # mol%, Input
        
        # FLUE GAS의 조성별 이론 에어량 및 EXCESS 에어량 계산 (m3/hr)
        FUEL_MWT = (np.dot(self.fuel_comp, self.arr_MWT))/100
        FUEL_MOLE =  FURN_Fuel_Flow / FUEL_MWT
        FURN_THEO_O2 = sum([FUEL_MOLE * (self.fuel_comp[i]/100) * self.arr_R2[i] * 22.4 for i in range(self.arr_length)])
        FURN_THEO_AIR = FURN_THEO_O2 / self.HOTAIR_O2_Frac
        FURN_THEO_CO2 = sum([FUEL_MOLE * (self.fuel_comp[i]/100) * self.arr_R3[i] * 22.4 for i in range(self.arr_length)]) + FURN_THEO_AIR * self.HOTAIR_CO2_Frac 
        FURN_THEO_H2O = sum([FUEL_MOLE * (self.fuel_comp[i]/100) * self.arr_R4[i] * 22.4 for i in range(self.arr_length)]) + FURN_THEO_AIR * self.HOTAIR_H2O_Frac
        FURN_THEO_N2 = FURN_THEO_AIR * self.HOTAIR_N2_Frac + (FUEL_MOLE * self.fuel_comp[25] / 100 * 22.4) 
        FURN_EXCESS_AIR = (FURN_THEO_CO2 + FURN_THEO_H2O + FURN_THEO_N2) * FURN_EXCESS_O2 / 100 / (self.HOTAIR_O2_Frac - FURN_EXCESS_O2 / 100)  

        # HOT AIR의 질량 유속 및 O2 조성 계산
        HOTAIR_V = FURN_THEO_AIR + FURN_EXCESS_AIR                  # m3/hr
        HOTAIR_CO2 = HOTAIR_V * self.HOTAIR_CO2_Frac / 22.4 * 44    # kg/hr
        HOTAIR_H2O = HOTAIR_V * self.HOTAIR_H2O_Frac / 22.4 * 18    # kg/hr
        HOTAIR_N2 = HOTAIR_V * self.HOTAIR_N2_Frac / 22.4 * 28      # kg/hr
        HOTAIR_O2 = HOTAIR_V * self.HOTAIR_O2_Frac / 22.4 * 32      # kg/hr
        HOTAIR_W = HOTAIR_CO2 + HOTAIR_H2O + HOTAIR_N2 + HOTAIR_O2  # kg/hr

        return HOTAIR_W, HOTAIR_O2

    def sim_obj(self):
        ''' What-if Simulation 최적화 문제 설정 '''
        # INPUT (test) ------------------------------
        self.furn_fgas_flow = [2038.4, 2267.4, 0, 0, 2690.8, 2677.1, 3431.2, 2765.5]
        self.air_list = [41.1, 33.4, 0.0, 0.0, 51.4, 24.2, 54.1, 44.0]
        self.fuel_bias_calc = [-90.6, -272.6, 0, 0, -336.2, -320.9, 457.2, -329.5]
        self.air_bias_calc = [-32.95, -27.44, 0, 0, -16.3, -45.15, -15.67, -39.4]
        # self.fuel_bias_calc = [-17.38, -23.87, -12.66, -24.45, -14.01, -21.77, -9.96, -23.99]
        # self.air_bias_calc = [-33.3, -27.8, 0.0, 0.0, 3.2, -45.6, -16.1, -40]
        furn_o2_arr = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
        # -------------------------------------

        # Bias 적용
        self.sim_air_list_bias = []     
        self.sim_fuel_list_bias = []    
        for i in range(8):
            if self.furn_fgas_flow == 0:                                                            # f/g = 0 일시, bias 미적용
                self.sim_air_list_bias.append(0)
            else:
                self.sim_air_list_bias.append(self.air_list[i] - (self.air_bias_calc[i]))      # Air 계산값 + Bias
        for i in range(8):
            if self.air_list == 0:                                                                  # air = 0 일시, bias 미적용
                self.sim_fuel_list_bias.append(0)     
            else:
                self.sim_fuel_list_bias.append(self.furn_fgas_flow[i] - self.fuel_bias_calc[i])     # Fuel 계산값 + Bias

        # 최적화
        init_point = [0.9, *furn_o2_arr]
        bnds = [(0,1), (-10,21), (-10,21), (-10,21), (-10,21), (-10,21), (-10,21), (-10,21), (-10,21)]
        from scipy.optimize import minimize
        res = minimize(self.sim_opt, init_point, method='Nelder-Mead', bounds=bnds, options={'fatol':0.0001, 'disp': False})
        
        # try:
        #     res = minimize(self.sim_opt, init_point, method='Nelder-Mead', bounds=bnds, options={'fatol':0.0001, 'disp': False})
        # except:
        #     pass
        
        self.FRESH_AIR = self.obj_HOT_AIR - self.EX_GTG_W
        air_d = sum(self.air_list) / (self.obj_HOT_AIR / 1000)

        result_tmp = [float(self.obj_HOT_AIR) / 1000 * air_d, float(self.EX_GTG_W) / 1000 * air_d, float(self.FRESH_AIR) / 1000 * air_d]

        for i in range(len(result_tmp)):
            self.result_list.append(result_tmp[i])

    def sim_opt(self, x):
        ''' What-if Simulation 최적화 실시 '''
        self.furn_air_list = []
        self.furn_o2_list = []

        mv1 = x[0]              
        self.hotair_calc(mv1)   # Cal. Hot Air composition
        self.mv_list = [x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]]
        # print(self.mv_list)
        # 8개 분해로 계산
        for i in range(0, 8):
            chk = self.furn_calc(i, self.mv_list[i])
            self.furn_air_list.append(chk[0])
            self.furn_o2_list.append(chk[1])
        furn_air_sum = sum(self.furn_air_list)
        furn_o2_sum = sum(self.furn_o2_list)

        fresh_o2_sum = furn_o2_sum - self.GTG_O2_W
        fresh_air_sum = fresh_o2_sum * 100 / 23

        obj_HOT_AIR = self.EX_GTG_W + fresh_air_sum     # Raw
        obj_FURN_AIR = furn_air_sum                     # Calc

        obj1_1 = (self.sim_air_list_bias[0] * 1000 - self.furn_air_list[0]) ** 2     # 목적함수, Air 입력값 vs Air 계산값
        obj1_2 = (self.sim_air_list_bias[1] * 1000 - self.furn_air_list[1]) ** 2
        obj1_3 = (self.sim_air_list_bias[2] * 1000 - self.furn_air_list[2]) ** 2
        obj1_4 = (self.sim_air_list_bias[3] * 1000 - self.furn_air_list[3]) ** 2
        obj1_5 = (self.sim_air_list_bias[4] * 1000 - self.furn_air_list[4]) ** 2
        obj1_6 = (self.sim_air_list_bias[5] * 1000 - self.furn_air_list[5]) ** 2
        obj1_7 = (self.sim_air_list_bias[6] * 1000 - self.furn_air_list[6]) ** 2
        obj1_8 = (self.sim_air_list_bias[7] * 1000 - self.furn_air_list[7]) ** 2

        obj2 = (obj_HOT_AIR - obj_FURN_AIR) ** 2                                    # 목적함수2, Hot Air vs 8기 분해로 Air
        print(obj1_1, obj1_2, obj1_3, obj1_4, obj1_5, obj1_6, obj1_7, obj1_8)
        obj1 = obj1_1 + obj1_2 + obj1_3 + obj1_4 + obj1_5 + obj1_6 + obj1_7 + obj1_8 
        obj = obj1 + obj2
        # print(self.sim_air_list_bias[0] * 1000, self.furn_air_list[0])

        
        self.obj_HOT_AIR = obj_HOT_AIR
        self.obj_FURN_AIR = obj_FURN_AIR

        return obj

    def insert_data_to_DB(self, insert_tup):

        conn = sqlite3.connect(self.db_name)
        cur = conn.cursor()

        sql = _make_insert_sql(table_name='case_simul_output', values_num=12)
        cur.execute(sql, insert_tup)
        conn.commit()
    
        conn.close()   

    def run(self):
        ''' What-if Case Simulation 실행 '''
        self.load_data()
        self.result_list = []
        self.gtg_calc()
        self.sim_obj()
        print(self.result_list)
        print(self.mv_list)
        # 결과값 저장
        # COMB_AIR = [self.result_list[0]]
        # COMB_AIR_O2 = [(self.HOTAIR_O2_Frac) * 100]
        # GTG_FLUE_GAS = [self.result_list[1]]
        # FRESH_AIR = [self.result_list[2]]
        # FURN_O2 = self.mv_list
        # insert_tup = tuple(COMB_AIR, COMB_AIR_O2, GTG_FLUE_GAS, FRESH_AIR, FURN_O2)
        # self.insert_data_to_DB(insert_tup = insert_tup)

if __name__ == "__main__":
    CASE_SIMULATION().run()
