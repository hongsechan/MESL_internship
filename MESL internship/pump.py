from thermo_properties import ThermoProperties

class Pump():
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

        

      


if __name__ == '__main__':
    p = Pump()

    print(f"Outlet temperature: {p.T_out_C:.2f} °C")
    print(f"Outlet enthalpy: {p.h_out:.2f} kJ/kg")
    print(f"Power input: {p.w_p:.2f} kJ/kg")
    print(f"Power input: {p.W_dot:.2f} kW")
