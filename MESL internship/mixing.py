from thermo_properties import ThermoProperties  
import numpy as np

class Mixing():
    def __init__(self, stream1=None, mass_flow1=None, stream2=None, mass_flow2=None):
        if stream1 is None:
            self.stream1 = ThermoProperties()
        else:
            self.stream1 = stream1

        if mass_flow1 is None:
            self.mass_flow1 = float(input("Enter mass flow rate of stream 1 (kg/s): "))
        else:
            self.mass_flow1 = mass_flow1

        if stream2 is None:
            self.stream2 = ThermoProperties()
        else:
            self.stream2 = stream2

        if mass_flow2 is None:
            self.mass_flow2 = float(input("Enter mass flow rate of stream 2 (kg/s): "))
        else:
            self.mass_flow2 = mass_flow2

        self.mole_flow1 = self.mass_flow1 / self.stream1.M_mix
        self.mole_flow2 = self.mass_flow2 / self.stream2.M_mix

        self.outlet_mole_flow_each = self.calculate_outlet_mole_flow_each()
        self.outlet_mole_fraction = self.calculate_outlet_mole_fraction()
        self.outlet_idx_comp = np.where(self.outlet_mole_fraction > 0)[0]
        self.outlet_mole_fraction_selected = self.outlet_mole_fraction[self.outlet_idx_comp]
    #----------------------------------------------------------------- 출구 혼합물의 몰분율 매서드
        self.M_out = self.calculate_M_out()
    #----------------------------------------------------------------- 출구 혼합물의 몰질량 매서드
        self.t1 = self.stream1.t
        self.t2 = self.stream2.t
        T_out = self.find_T_out()
        if T_out is None:
            self.T_out = None
        else:
            self.T_out = T_out -273.15
    #----------------------------------------------------------------- 출구 혼합물의 온도 매서드
        self.p_out = min(self.stream1.p, self.stream2.p)
    #----------------------------------------------------------------- 출구 혼합물의 압력 매서드

    def calculate_outlet_mole_flow_each(self):
        outlet_mole_flow_each = np.zeros(len(self.stream1.element))
        outlet_mole_flow_each[self.stream1.idx_comp] += self.mole_flow1 * self.stream1.mole_fraction
        outlet_mole_flow_each[self.stream2.idx_comp] += self.mole_flow2 * self.stream2.mole_fraction
        return outlet_mole_flow_each
    
    def calculate_outlet_mole_fraction(self):
        return self.outlet_mole_flow_each / np.sum(self.outlet_mole_flow_each) 


    #----------------------------------------------------------------- 출구 혼합물의 몰유량 계산
    def calculate_M_out(self):
        return np.sum(self.outlet_mole_fraction * self.stream1.M)
    #----------------------------------------------------------------- 출구 혼합물의 몰질량 계산

    def h_out_at(self, T):
        if T >= 1000:
            coeff = self.stream1.coeff_data[self.outlet_idx_comp, 0, :]
        else:
            coeff = self.stream1.coeff_data[self.outlet_idx_comp, 1, :]

        
        t_ref = self.stream1.t_ref

        T_delta_h_calc = np.array([T-t_ref,
                                    (T*T-t_ref*t_ref)/2,
                                    (T*T*T-t_ref*t_ref*t_ref)/3,
                                    (T*T*T*T-t_ref*t_ref*t_ref*t_ref)/4,
                                    (T*T*T*T*T-t_ref*t_ref*t_ref*t_ref*t_ref)/5])   
        delta_h = 8.314*(coeff[:, :5]*T_delta_h_calc).sum(axis=1)

        h_f = self.stream1.h_f
        idx_comp = self.outlet_idx_comp
        h_each_mole = h_f[idx_comp]+delta_h
        mole_fraction = self.outlet_mole_fraction_selected
        M_mix = self.M_out
        return np.sum(h_each_mole*mole_fraction)/M_mix
    #----------------------------------------------------------------- 출구 혼합물의 엔탈피 계산

    def energy_balance(self, T_out):
        inlet_energy = self.mass_flow1 * self.stream1.h + self.mass_flow2 * self.stream2.h
        outlet_energy = (self.mass_flow1 + self.mass_flow2) * self.h_out_at(T_out)
        return inlet_energy - outlet_energy
    

    def find_T_out(self, tol = 1e-6 ):
        T_low = min(self.t1, self.t2)
        T_high = max(self.t1, self.t2)

        f_low = self.energy_balance(T_low)
        f_high = self.energy_balance(T_high)

        if T_low == T_high:
            return T_low

        elif f_low * f_high >0 :
            print ("Error")
            return None
        
        else:
            for _ in range(100):
                T_mid = (T_low + T_high)/2
                f_mid = self.energy_balance(T_mid)

                if abs(f_mid) < tol:
                    return T_mid
                elif abs(f_high) < tol:
                    return T_high
                elif abs(f_low) < tol:
                    return T_low
                elif f_low*f_mid < 0:
                    T_high = T_mid
                    f_high = f_mid
                else: 
                    T_low = T_mid
                    f_low = f_mid

            return T_mid
        
    #----------------------------------------------------------------- 출구 혼합물의 온도 계산


if __name__ == '__main__':
    mixing = Mixing()

    print(f"Mass flow rate of stream 1: {mixing.mass_flow1} kg/s")
    print(f"Mass flow rate of stream 2: {mixing.mass_flow2} kg/s")
    print(f"Molar mass of the mixed stream: {mixing.M_out:.4f} kg/kmol")
    for i in mixing.outlet_idx_comp:
        print(f"Mole fraction of {mixing.stream1.element[i]} in the mixed stream: {mixing.outlet_mole_fraction[i]:.4f}")

    if mixing.T_out is None:
        print("T_out could not be calculated.")
    else:
        print(f"Temperature of the mixed stream: {mixing.T_out:.2f} °C")
    print(f"Pressure of the mixed stream: {mixing.p_out:.2f} atm")        
        
        





