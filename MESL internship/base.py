import numpy as np

class BasicInfo:
    def __init__(self):
        self.comp_name = input("Enter composition(H2, H2O, CO, CO2, O2, N2, CH4): ").split()
        self.p = float(input("Enter pressure (atm): "))
        self.t = float(input("Enter temperature (C): "))
        self.t = self.t + 273.15
        self.element = np.array(["H2", "H2O", "CO", "CO2", "O2", "N2", "CH4"])
        self.comp = np.array(self.comp_name)


if __name__ == "__main__":
    b= BasicInfo()
    print(b.comp)
    print(b.p)
    print(b.t)
