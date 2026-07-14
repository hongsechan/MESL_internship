from thermo_properties import ThermoProperties

class Compression():
    def __init__(self, stream=None, p_out=None, efficiency=None, mass_flow=None):
        if stream is None:
            self.stream = ThermoProperties()
        else:
            self.stream = stream

        if p_out is None:
            self.p_out = float(input("Enter the outlet pressure (atm): "))
        else:
            self.p_out = p_out

        if efficiency is None:
            self.efficiency = float(input("Enter the isentropic efficiency (0-1): "))
        else:
            self.efficiency = efficiency        

        if mass_flow is None:
            self.mass_flow = float(input("Enter the mass flow rate (kg/s): "))
        else:
            self.mass_flow = mass_flow

        self.M_out = self.stream.M_mix
        self.outlet_mole_fraction = self.stream.mole_fraction

        if self.liquid_water():
            self.caculate_pump()
        else:
            self.calculate_compressor()
    
    def liquid_water(self):
            T_C = self.stream.t - 273.15 
            return self.stream.comp_name == ["H2O"] and len(self.stream.comp_name) == 1 and T_C < 100 
    #-------------------------------------------------------------- liquid water 인지 확인
    def caculate_pump(self):
        self.rho = 997
        self.v = 1/self.rho
        self.cp_liquid = 4.18
    
        w_p_s = self.v * (self.p_out - self.stream.p) * 101.325
        self.w_p= w_p_s / self.efficiency

        self.h_out = self.stream.h + self.w_p

        delta_T = self.w_p / self.cp_liquid
        self.T_out = self.stream.t + delta_T
        self.T_out_C = self.T_out - 273.15

        self.W_dot = self.mass_flow * self.w_p

        self.Ts_out = None
        self.Ts_out_C = None
        self.hs_out = None
    #-------------------------------------------------------------- pump 계산    

    def calculate_compressor(self): 
        
        self.Ts_out = self.calculate_Ts_out()
        self.Ts_out_C = self.Ts_out - 273.15
        #--------------------------------------------------------------- 출구 등엔트로피 온도 매서드    
        self.hs_out = self.stream.h_at(self.Ts_out)
        self.h_out = self.stream.h + (self.hs_out - self.stream.h) / self.efficiency

        #-------------------------------------------------------------- 출구 엔탈피 매서드    
        self.T_out = self.calculate_T_at_h(self.h_out)
        self.T_out_C = self.T_out - 273.15
        #-------------------------------------------------------------- 출구 온도 매서드
        self.W_dot = self.mass_flow * (self.h_out - self.stream.h)

        #-------------------------------------------------------------- 압축기 전력 매서드

    def entropy_balance(self, T):
        return self.stream.s_at(T, self.p_out) - self.stream.s
    
    def calculate_Ts_out(self, tol = 1e-6):
        T_low = self.stream.t
        T_high = self.stream.t*2

        f_low = self.entropy_balance(T_low)
        f_high = self.entropy_balance(T_high)


        while f_low * f_high > 0:
            T_high *= 2
            f_high = self.entropy_balance(T_high)

        for _ in range(100):
            T_mid = (T_low + T_high) / 2
            f_mid = self.entropy_balance(T_mid)

            if abs(f_mid) < tol:
                return T_mid
            elif f_low * f_mid < 0:
                T_high = T_mid
                f_high = f_mid
            else:
                T_low = T_mid
                f_low = f_mid

        return T_mid
    #-------------------------------------------------------------- 출구 등엔트로피 온도 계산
    def calculate_T_at_h(self, h, tol=1e-6):
        T_low = self.Ts_out
        T_high = self.Ts_out * 2

        f_low = self.stream.h_at(T_low) - h
        f_high = self.stream.h_at(T_high) - h

        if abs(f_low) < tol:
            return T_low

        while f_low * f_high > 0:
            T_high *= 2
            f_high = self.stream.h_at(T_high) - h

        for _ in range(100):
            T_mid = (T_low + T_high) / 2
            f_mid = self.stream.h_at(T_mid) - h

            if abs(f_mid) < tol:
                return T_mid
            elif f_low * f_mid < 0:
                T_high = T_mid
                f_high = f_mid
            else:
                T_low = T_mid
                f_low = f_mid

        return T_mid
    #-------------------------------------------------------------- 출구 엔탈피로 실제 온도 계산



if __name__ == '__main__':
    a= Compression()
    
    if a.Ts_out_C is not None:
        print(f"Outlet temperature (isentropic): {a.Ts_out_C:.2f} °C")
    else: pass
    
    print(f"Outlet temperature (actual): {a.T_out_C:.2f} °C")
    
    if a.hs_out is not None:
        print(f"Outlet enthalpy (isentropic): {a.hs_out:.2f} kJ/kg")
    else: pass

    print(f"Outlet enthalpy (actual): {a.h_out:.2f} kJ/kg")
    print(f"Power input: {a.W_dot:.2f} kW")
    print(f"Outlet mole fraction: {a.outlet_mole_fraction}")
    print(f"Molar mass of the outlet stream: {a.M_out:.4f} kg/kmol")