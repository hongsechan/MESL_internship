from thermo_properties import ThermoProperties


class Splitter():
    def __init__(self, stream=None, mass_flow=None, split_ratio=None):
        if stream is None:
            self.stream = ThermoProperties()
        else:
            self.stream = stream

        if mass_flow is None:
            self.mass_flow = float(input("Enter mass flow rate of the inlet stream (kg/s): "))
        else:
            self.mass_flow = mass_flow

        if split_ratio is None:
            self.split_ratio = float(input("Enter split ratio (0-1): "))
        else:
            self.split_ratio = split_ratio
       
        self.outlet_mole_fraction_stream1 = self.stream.mole_fraction 
        self.outlet_mole_fraction_stream2 = self.stream.mole_fraction

        self.outlet_mass_flow1, self.outlet_mass_flow2 = self.calculate_outlet_mass_flows()
  
    def calculate_outlet_mass_flows(self): 
        outlet_mass_flow1 = self.mass_flow * self.split_ratio 
        outlet_mass_flow2 = self.mass_flow * (1 - self.split_ratio)
        return outlet_mass_flow1, outlet_mass_flow2
