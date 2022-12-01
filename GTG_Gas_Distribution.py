import pandas as pd
import numpy as np
import json
from scipy.optimize import minimize

df = pd.read_excel('raw_data.xlsx')

class AIR_BALANCE:
    def DATA_TUNING(self, i):
        ''' 밸런스 수립 및 Raw Data 튜닝 실시 '''
        ## GTG Data
        param_gtg_fgas = 0.9    # GTG 연료 유량 보정(향후 데이터분석 필요, 연료 생산(NC1FY63007) vs 연료 사용(NC1FIT600 + NC1G96DG_1) 차이 고려 파라미터값 설정함)
        # self.gtg_fgas_flow = df['GTGFGAS'][i] * param_gtg_fgas   # kg/hr, GTG 연료 유량, NC1FIT600
        # self.gtg_stack_o2 = df['GTGO2'][i] # mol%, GTG 배기가스 O2%, NC1AI1010

        self.gtg_fgas_flow = 150 * param_gtg_fgas   # kg/hr, GTG 연료 유량, NC1FIT600
        self.gtg_stack_o2 = 20.6 # mol%, GTG 배기가스 O2%, NC1AI1010


        ## 향후 GTG ON/OFF 시그널로 대체 필요
        if self.gtg_stack_o2 > 18 or self.gtg_fgas_flow <= 1000:        # 향후 GTG ON/OFF 시그널로 대체 필요
            self.gtg_fgas_flow = 0.00001    # 오류 방지(GTG 연료 유량)
            # self.gtg_stack_o2 = 21.00001    # 오류 방지(GTG 배기가스 O2%)
            self.gtg_onoff = 0              # GTG ON/OFF 확인 (1: ON, 0: OFF)
        else:
            self.gtg_onoff = 1

        ## 향후 모델 Coefficient 값은 외부에서 불러오도록 수정 필요함
        coeff_on_lpg = [0.149947477, 0.168400533, 0.192113600, 0.214467051, 0.226178363, 0, 0.236999310, 0]     # 분해로 별 Feed Model(Fuel Gas 계산, GTG ON & LPG)
        coeff_off_lpg = [0.158177880, 0.188621457, 0, 0.215316451, 0.222893098, 0, 0.272555958, 0]              # 분해로 별 Feed Model(Fuel Gas 계산, GTG OFF & LPG)
        coeff_on_naph = [0, 0, 0.165161442, 0.167588885, 0.167531435, 0.170951644, 0.167688808, 0.175201227]    # 분해로 별 Feed Model(Fuel Gas 계산, GTG ON & NAPH)
        coeff_off_naph = [0, 0, 0.191086941, 0.187579468, 0.196766548, 0.194852126, 0.193308698, 0.200981580]   # 분해로 별 Feed Model(Fuel Gas 계산, GTG OFF & NAPH)

        coeff_on_ethane = [0.149947477, 0.168400533, 0.192113600, 0.214467051, 0.226178363, 0, 0.236999310, 0]     # 분해로 별 Feed Model(Fuel Gas 계산, GTG ON & ETHANE)
        coeff_off_ethane = [0.158177880, 0.188621457, 0, 0.215316451, 0.222893098, 0, 0.272555958, 0]  # 분해로 별 Feed Model(Fuel Gas 계산, GTG OFF & ETHANE)

        ## LPG 유량 계산(KG/HR), 8호기 없음
        # 1호기: NC1FC1117 + NC1FC1119 + NC1FC1121 + NC1FC1123,     2호기: NC1FC1217 + NC1FC1219 + NC1FC1221 + NC1FC1223
        # 3호기: NC1FC1317 + NC1FC1319 + NC1FC1321 + NC1FC1323,     4호기: NC1FC1417 + NC1FC1419 + NC1FC1421 + NC1FC1423
        # 5호기: NC1FC61517 + NC1FC61519 + NC1FC61521 + NC1FC61523, 6호기: NC1FC61617 + NC1FC61619 + NC1FC61621 + NC1FC61623
        # 7호기: NC1FC61717 + NC1FC61719 + NC1FC61721 + NC1FC61723

        ## NAPHTHA 유량 계산(KG/HR), 1,2호기 없음
        # 3호기: NC1FC1301 + NC1FC1303 + NC1FC1305 + NC1FC1307,     4호기: NC1FC1401 + NC1FC1403 + NC1FC1405 + NC1FC1407
        # 5호기: NC1FC1501 + NC1FC1503 + NC1FC1505 + NC1FC1507,     6호기: NC1FC1601 + NC1FC1603 + NC1FC1605 + NC1FC1607
        # 7호기: NC1FC1701 + NC1FC1703 + NC1FC1705 + NC1FC1707,     8호기: NC1FC1801 + NC1FC1803 + NC1FC1805 + NC1FC1807

        Feed_LPG = [df['LPG_1'][i], df['LPG_2'][i], df['LPG_3'][i], df['LPG_4'][i], df['LPG_5'][i], df['LPG_6'][i], df['LPG_7'][i], 0]
        Feed_NAPH = [0, 0, df['NAPH_3'][i], df['NAPH_4'][i], df['NAPH_5'][i], df['NAPH_6'][i], df['NAPH_7'][i], df['NAPH_8'][i]]
        Feed_ETHANE = [0, 0, df['LPG_3'][i], df['LPG_4'][i], 0, 0, 0, 0]



        ## 분해로 Fuel Gas(KG/HR)
        # 1호기: NC1FC1127, 2호기: NC1FC1227, 3호기: NC1FC1327, 4호기: NC1FC1427
        # 5호기: NC1FC1527, 6호기: NC1FC1627, 7호기: NC1FC1727, 8호기: NC1FC1827
        Fuel_Raw = [df['FGAS_1'][i], df['FGAS_2'][i], df['FGAS_3'][i], df['FGAS_4'][i], df['FGAS_5'][i], df['FGAS_6'][i], df['FGAS_7'][i], df['FGAS_8'][i]]
        Fuel_LPG = []; Fuel_NAPH = []; Fuel_ETHANE = []; Fuel_Feed = [];

        # LPG or Naphtha 원료 사용 여부 -> 분해로 Fuel Gas 유량 계산(Feed 모델)
        # 향후 Feed 판단 로직 추가 필요(LPG? Naph? Ethane?)
        for j in range(8):
            if self.gtg_onoff == 1:
                Fuel_LPG.append(coeff_on_lpg[j] * Feed_LPG[j])
                Fuel_NAPH.append(coeff_on_naph[j] * Feed_NAPH[j])
                Fuel_ETHANE.append(coeff_on_ethane[j] * Feed_ETHANE[j])
                Fuel_Feed.append(max(Fuel_LPG[j], Fuel_NAPH[j], Fuel_ETHANE[j]))
                if Fuel_Feed[j] <= 500:
                    Fuel_Feed[j] = Fuel_Raw[j]
            elif self.gtg_onoff == 0:
                Fuel_LPG.append(coeff_off_lpg[j] * Feed_LPG[j])
                Fuel_NAPH.append(coeff_off_naph[j] * Feed_NAPH[j])
                Fuel_ETHANE.append(coeff_off_ethane[j] * Feed_ETHANE[j])
                Fuel_Feed.append(max(Fuel_LPG[j], Fuel_NAPH[j], Fuel_ETHANE[j]))
                if Fuel_Feed[j] <= 500:
                    Fuel_Feed[j] = Fuel_Raw[j]

        ## Furnace Data
        self.furn_fgas_flow = Fuel_Feed # kg/hr, Furnace F/Gas Flow(Feed Model)
        # self.furn_fgas_flow = Fuel_Raw
        self.furn_fgas_flow = [self.furn_fgas_flow[j] * df['RATIO'][i] for j in range(8)]   # Ratio, 분해로 Fuel Gas Density 보정, 1.97(HYSYS, 고정 or 계산) / NC1AI1002

        fuel_bias = [raw-cal for raw,cal in zip(Fuel_Raw, Fuel_Feed)]
        print(fuel_bias)
        ## 분해로 Air 유량(delta 계산용)
        self.Air_Raw_list = [df['AIR_1'][i], df['AIR_2'][i], df['AIR_3'][i], df['AIR_4'][i], df['AIR_5'][i], df['AIR_6'][i], df['AIR_7'][i], df['AIR_8'][i]]
        self.Air_Raw_Sum = df['AIR_1'][i] + df['AIR_2'][i] + df['AIR_3'][i] + df['AIR_4'][i] + df['AIR_5'][i] + df['AIR_6'][i] + df['AIR_7'][i] + df['AIR_8'][i]

        ## 분해로 Excess O2(mol%)
        # 1호기: NC1AI1103A, 2호기: NC1AI1203A, 3호기: NC1AI1303A, 4호기: NC1AI1403A
        # 5호기: NC1AI1503A, 6호기: NC1AI1603A, 7호기: NC1AI1703A, 8호기: NC1AI1803A
        self.furn_o2 = [df['O2_1'][i], df['O2_2'][i], df['O2_3'][i], df['O2_4'][i], df['O2_5'][i], df['O2_6'][i], df['O2_7'][i], df['O2_8'][i]] # Excess O2 mol%, Furnace O2%(가상센서)
        self.furn_o2_target = [df['O2_1'][i], df['O2_2'][i], df['O2_3'][i], df['O2_4'][i], df['O2_5'][i], df['O2_6'][i], df['O2_7'][i], df['O2_8'][i]]  # Excess O2 mol%, Furnace O2% Target
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


    def GTG_CALC(self):
        ''' GTG 배기가스 유량 및 조성 계산(PIS Tag 필요) '''
        ## GTG Fuel Flow, Comp.(Mol%)
        self.gtg_fgas_flow = 0 * 0.9
        self.gtg_stack_o2 = 20.6
        GTG_Fuel_Flow = self.gtg_fgas_flow  # kg/hr
        GTG_Fuel_Hydrogen = 24.2622         # H2
        GTG_Fuel_CO = 0.1219                # CO
        GTG_Fuel_Methane = 75.4055          # CH4
        GTG_Fuel_Ethane = 0.0011            # C2H6
        GTG_Fuel_Ethylene = 0.1791          # C2H4
        GTG_Fuel_Acetylene = 0.0197         # C2H2

        GTG_Stack_O2 = self.gtg_stack_o2  # mol%
        GTG_Fuel_MW = (GTG_Fuel_Hydrogen * 2.0160 + GTG_Fuel_CO * 28.01 + GTG_Fuel_Methane * 16.0429 + GTG_Fuel_Ethane * 30.0699 + GTG_Fuel_Ethylene * 28.0538 + GTG_Fuel_Acetylene * 26.0380) / 100 # kg/kmol

        ## GTG 배기가스 유량 계산(m3/hr)
        GTG_O2 = GTG_Fuel_Flow / GTG_Fuel_MW * (GTG_Fuel_Hydrogen / 2 + GTG_Fuel_CO / 2 + GTG_Fuel_Methane * 2 + GTG_Fuel_Ethane * 3.5 + GTG_Fuel_Ethylene * 3 + GTG_Fuel_Acetylene * 2.5) / 100 * 22.4
        GTG_CO2 = GTG_Fuel_Flow / GTG_Fuel_MW * (GTG_Fuel_Hydrogen * 0 + GTG_Fuel_CO * 1 + GTG_Fuel_Methane * 1 + GTG_Fuel_Ethane * 2 + GTG_Fuel_Ethylene * 2 + GTG_Fuel_Acetylene * 2) / 100 * 22.4
        GTG_H2O = GTG_Fuel_Flow / GTG_Fuel_MW * (GTG_Fuel_Hydrogen * 1 + GTG_Fuel_CO * 0 + GTG_Fuel_Methane * 2 + GTG_Fuel_Ethane * 3 + GTG_Fuel_Ethylene * 2 + GTG_Fuel_Acetylene * 1) / 100 * 22.4
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

    def HOTAIR_CALC(self, a):
        ''' HOT AIR 조성 계산'''
        ## Air Ratio Calc(EXHAUST FROM GTG + FRESH AIR)
        Air_Ratio = a
        self.HOTAIR_O2_Frac = self.GTG_O2_Frac * Air_Ratio + 0.21 * (1 - Air_Ratio)
        self.HOTAIR_N2_Frac = self.GTG_N2_Frac * Air_Ratio + 0.79 * (1 - Air_Ratio)
        self.HOTAIR_CO2_Frac = self.GTG_CO2_Frac * Air_Ratio
        self.HOTAIR_H2O_Frac = self.GTG_H2O_Frac * Air_Ratio

    def FURNACE_CALC(self, a):
        ''' 분해로 공기 유량 계산 '''
        ## Furnace Fuel Comp. (PIS Tag 필요)


        FURN_Fuel_Flow = self.furn_fgas_flow[a]  # kg/hr, Input
        FURN_EXCESS_O2 = self.furn_o2[a] # mol%, Input
        print(FURN_Fuel_Flow)
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

    def obj(self):
        '''Hot Air와 분해로 Air의 밸런스 수립을 위한 최적화 문제 설정 '''
        init_point = [0.9]
        bnds = [(0,1)]
        res = minimize(self.opt, init_point, method='Nelder-Mead', bounds=bnds, options={'fatol':0.0000001, 'disp': False})
        self.FRESH_AIR = self.obj_HOT_AIR - self.EX_GTG_W

        air_d = self.Air_Raw_Sum / (self.obj_HOT_AIR / 1000)    # 8기 분해로 측정값 vs 8기 이론 공기량

        result_tmp = [self.obj_HOT_AIR / 1000 * air_d, self.EX_GTG_W / 1000 * air_d, self.FRESH_AIR / 1000 * air_d ,
                      *self.furn_air_list, *self.HOTAIR_CO2_Frac * 100, *self.HOTAIR_H2O_Frac * 100, *self.HOTAIR_N2_Frac * 100, *self.HOTAIR_O2_Frac * 100]

        for i in range(len(result_tmp)):
            self.result_list.append(result_tmp[i])
        for i in range(len(self.furn_fgas_flow)):
            self.result_list.append(self.furn_fgas_flow[i])

    def opt(self, x):
        ''' Hot Air와 분해로 Air의 밸런스 수립을 위한 최적화 실시 '''
        self.furn_air_list = []
        self.furn_o2_list = []

        mv1 = x # Hot Air Ratio
        self.HOTAIR_CALC(x)

        # 8개 분해로 계산
        for i in range(0, 8):
            chk = self.FURNACE_CALC(i)
            self.furn_air_list.append(chk[0][0])
            self.furn_o2_list.append(chk[1][0])
        furn_air_sum = sum(self.furn_air_list)
        furn_o2_sum = sum(self.furn_o2_list)

        fresh_o2_sum = furn_o2_sum - self.GTG_O2_W
        fresh_air_sum = fresh_o2_sum * 100 / 23

        obj_HOT_AIR = self.EX_GTG_W + fresh_air_sum
        obj_FURN_AIR = furn_air_sum

        obj = (obj_HOT_AIR - obj_FURN_AIR) ** 2
        print(obj_HOT_AIR, obj_FURN_AIR)
        print(obj_FURN_AIR)
        self.obj_HOT_AIR = obj_HOT_AIR
        self.obj_FURN_AIR = obj_FURN_AIR

        return obj

    def run(self):
        result_table = pd.DataFrame()
        for i in range(0, 1): 
            self.result_list = []
            self.DATA_TUNING(i)
            self.GTG_CALC()
            self.obj()

            a = pd.DataFrame(self.result_list)
            print(self.result_list)
            a = a.transpose()
            result_table = pd.concat([result_table, a])

        result_table.columns = ["C_HOTAIR_F", "C_GTGAIR_F", "C_FRESHAIR_F", "C_FURNAIR_1", "C_FURNAIR_2", "C_FURNAIR_3",
                                "C_FURNAIR_4", "C_FURNAIR_5", "C_FURNAIR_6", "C_FURNAIR_7", "C_FURNAIR_8",
                                "C_HOTAIRCO2", "C_HOTAIRH2O", "C_HOTAIRN2", "C_HOTAIRO2",
                                "C_FURNFUEL_1", "C_FURNFUEL_2", "C_FURNFUEL_3", "C_FURNFUEL_4", "C_FURNFUEL_5",
                                "C_FURNFUEL_6", "C_FURNFUEL_7", "C_FURNFUEL_8"]

        result_table.to_csv('result.csv')

if __name__ == "__main__":
    AIR_BALANCE().run()


