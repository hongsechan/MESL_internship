import numpy as np

class BasicInfo:
    def __init__(self, comp_name=None, p=None, t_C=None):
        if comp_name is None:
            comp_name = input("Enter composition(H2, H2O, CO, CO2, O2, N2, CH4): ").split()
        
        if p is None:
           p = float(input("Enter pressure (atm): "))
      
        if t_C is None:
            t_C = float(input("Enter temperature (C): "))
    
        self.comp_name = comp_name
        self.p = p
        self.t = t_C + 273.15
        self.element = np.array(["H2", "H2O", "CO", "CO2", "O2", "N2", "CH4"])
        self.comp = np.array(comp_name)


if __name__ == "__main__":
    b= BasicInfo()
    print(b.comp)
    print(b.p)
    print(b.t)
    