
import numpy as np

class MassSolver():
    def __init__(self, stream_number = None):
        if stream_number is None:
            self.stream_numbers = int(input("Enter the number of streams: "))
        else:
            self.stream_numbers = stream_number

        self.right_hand_side = np.zeros(self.stream_numbers)
        self.stream_name = np.array(["Stream " + str(i+1) for i in range(self.stream_numbers)])

        self.equation = np.zeros((self.stream_numbers, self.stream_numbers))

        self.current_equation_count = 0


    def add_mixing_eq(self, inlet_stream_numbers = None, outlet_stream_number = None):
        if inlet_stream_numbers is None:
            self.inlet_stream_numbers = np.array([int(x) for x in input("Enter the inlet stream numbers (comma-separated): ").split(",")])
        else:
            self.inlet_stream_numbers = np.array(inlet_stream_numbers) 

        if outlet_stream_number is None:
            self.outlet_stream_number = int(input("Enter the outlet stream number: ")) 
        else:
            self.outlet_stream_number = outlet_stream_number 

        inlet_stream_idxs = self.inlet_stream_numbers - 1
        outlet_stream_idx = self.outlet_stream_number - 1

        coeff = np.zeros(self.stream_numbers)

         # -m_in1 - m_in2 + m_out = 0
        coeff[inlet_stream_idxs] = +1.0
        coeff[outlet_stream_idx] = -1.0

        row = self.current_equation_count
        self.equation[row, :] = coeff
        self.current_equation_count += 1

        return self.equation

    def add_split_eq(self, inlet_stream_number = None, outlet_stream_number1 = None, outlet_stream_number2 = None, split_ratio = None):
        if inlet_stream_number is None:
            self.inlet_stream_number = int(input("Enter the inlet stream number: "))
        else:
            self.inlet_stream_number = inlet_stream_number

        if outlet_stream_number1 is None:
            self.outlet_stream_number1 = int(input("Enter the first outlet stream number: "))
        else:
            self.outlet_stream_number1 = outlet_stream_number1

        if outlet_stream_number2 is None:
            self.outlet_stream_number2 = int(input("Enter the second outlet stream number: "))
        else:
            self.outlet_stream_number2 = outlet_stream_number2

        if split_ratio is None:
            self.split_ratio =  float(input("Enter the split ratio for the outlet stream: "))
        else:
            self.split_ratio = split_ratio


        inlet_stream_idx = self.inlet_stream_number - 1
        outlet_stream_idx1 = self.outlet_stream_number1 - 1
        outlet_stream_idx2 = self.outlet_stream_number2 - 1


        # m_out1 - r*m_in = 0
        coeff1 = np.zeros(self.stream_numbers)
        coeff1[inlet_stream_idx] = +1.0 * self.split_ratio
        coeff1[outlet_stream_idx1] = -1.0
        row = self.current_equation_count
        self.current_equation_count += 1
        self.equation[row, :] = coeff1

        # m_out2 - (1-r)*m_in = 0
        coeff2 = np.zeros(self.stream_numbers)
        coeff2[inlet_stream_idx] = +1.0 * (1 - self.split_ratio)
        coeff2[outlet_stream_idx2] = -1.0
        row = self.current_equation_count
        self.current_equation_count += 1
        self.equation[row, :] = coeff2

        return self.equation

    def add_known_mass_flow(self, stream_number = None, mass_flow = None):
        if stream_number is None:
            self.stream_number = int(input("Enter the stream number with known mass flow: "))
        else:
            self.stream_number = stream_number

        if mass_flow is None:
            self.mass_flow = float(input("Enter the known mass flow (kg/s): "))
        else:
            self.mass_flow = mass_flow

        stream_idx = self.stream_number - 1
        coeff = np.zeros(self.stream_numbers)
        coeff[stream_idx] = 1.0
        self.right_hand_side[stream_idx] = self.mass_flow

        row = self.current_equation_count
        self.equation[row, :] = coeff
        self.current_equation_count += 1

        return self.equation, self.right_hand_side   