from mixing import Mixing
from mass_splitter import Splitter
from thermo_properties import ThermoProperties

class Model2():
    def __init__(self, stream1 = None, mass_flow1 = None, split_ratio = None):
        if stream1 is None:
            self.stream1 = ThermoProperties()
        else: 
            self.stream1 = stream1

        if mass_flow1 is None:
            self.mass_flow1 = float(input("Enter mass flow rate of stream 1 (kg/s): "))

        else:
            self.mass_flow1 = mass_flow1    

        if split_ratio is None:
            self.split_ratio = float(input("Enter split ratio (0-1): "))
        else:
            self.split_ratio = split_ratio

        self.tolerance = 1e-6
        self.max_iterations = 100
        self.initial_mass_flow3 = self.mass_flow1 * self.split_ratio

        self.solve()

    def solve(self):

            stream3_guess = self.stream1    
            mass_flow3_guess = self.initial_mass_flow3

            for _ in range(self.max_iterations):
                mixing = Mixing(stream1=self.stream1, mass_flow1=self.mass_flow1, stream2= stream3_guess, mass_flow2= mass_flow3_guess)
                splitter = Splitter(stream=mixing.stream_out, mass_flow=mixing.mass_flow_out, split_ratio=self.split_ratio)

                mass_flow3_calculated = splitter.outlet_mass_flow1
                mass_flow4_calculated = splitter.outlet_mass_flow2
                residual = mass_flow3_calculated - mass_flow3_guess

                if abs(residual) < self.tolerance:

                    self.stream2 = mixing.stream_out
                    self.mass_flow2 = mixing.mass_flow_out
                    self.stream3 = splitter.stream
                    self.mass_flow3 = mass_flow3_calculated 
                    self.stream4 = splitter.stream
                    self.mass_flow4 = mass_flow4_calculated
                    return    
                
                stream3_guess = splitter.stream
                mass_flow3_guess = mass_flow3_guess + residual * 0.5

            raise ValueError("Model2 did not converge within the maximum number of iterations.")
