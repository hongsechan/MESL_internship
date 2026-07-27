from thermo_properties import ThermoProperties
import constants as const

class HeatExchanger():
    def __init__(self, hot_fluid = None, cold_fluid = None, hot_fluid_massflow = None, cold_fluid_massflow = None, cold_flow_target_temp = None):
        if hot_fluid is None:
            self.hot_fluid = ThermoProperties()
        else:
            self.hot_fluid = hot_fluid

        if cold_fluid is None:
            self.cold_fluid = ThermoProperties()
        else:
            self.cold_fluid = cold_fluid

        if hot_fluid_massflow is None:
            self.hot_fluid_massflow = float(input("Enter the mass flow rate of hot fluid (kg/s): "))
        else:
            self.hot_fluid_massflow = hot_fluid_massflow

        if cold_fluid_massflow is None:
            self.cold_fluid_massflow = float(input("Enter the mass flow rate of cold fluid (kg/s): "))
        else:
            self.cold_fluid_massflow = cold_fluid_massflow

        if cold_flow_target_temp is None:
            self.cold_flow_target_temp = float(input("Enter the target temperature of cold fluid (C): ")) + const.kelvin_offset
        else:
            self.cold_flow_target_temp = cold_flow_target_temp + const.kelvin_offset

        self.h_out_cold = self.cold_fluid.h_at(self.cold_flow_target_temp)
        self.h_out_hot = self.calculate_h_out_hot()
        self.heat_duty = self.cold_fluid_massflow * (self.h_out_cold - self.cold_fluid.h)

        self.hot_flow_out_temp = self.calculate_T_at_h()
        self.hot_flow_out_temp_C =self.hot_flow_out_temp - const.kelvin_offset


    def calculate_h_out_hot(self):
        h_out_hot = self.hot_fluid.h + (self.cold_fluid_massflow * ( self.cold_fluid.h - self.h_out_cold)) / self.hot_fluid_massflow
        return h_out_hot

    
    def calculate_T_at_h(self, tol =1e-6):
        T_low = self.cold_flow_target_temp
        T_high = self.hot_fluid.t

        f_low = self.hot_fluid.h_at(T_low) - self.h_out_hot
        f_high = self.hot_fluid.h_at(T_high) - self.h_out_hot

        for _ in range(100):
            T_mid = (T_low + T_high)/2
            f_mid = self.hot_fluid.h_at(T_mid) - self.h_out_hot

            if abs(f_mid) < tol:
                return T_mid

            if f_low*f_mid < 0:
                T_high = T_mid
                f_high = f_mid

            else:
                T_low = T_mid
                f_low = f_mid
        return T_mid
            


    

if __name__ == "__main__":
    heat_exchanger = HeatExchanger()
    print(f"Heat duty: {heat_exchanger.heat_duty:.2f} kW")
    print(f"Hot fluid outlet enthalpy: {heat_exchanger.h_out_hot:.2f} kJ/kg")
    print(f"Cold fluid outlet enthalpy: {heat_exchanger.h_out_cold:.2f} kJ/kg")
    print(f"Hot fluid outlet temperature: {heat_exchanger.hot_flow_out_temp_C:.2f} C")




        