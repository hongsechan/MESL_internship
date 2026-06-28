from unittest.mock import patch 
from thermo_properties import ThermoProperties




def run_case(name, inputs):
    with patch("builtins.input", side_effect=inputs):
        a = ThermoProperties()

    print(f"Case: {name}")
    print(f"Mole Fraction: {a.mole_fraction}")
    print(f"Mixing Molecular Weight: {a.M_mix} (kg/kmol)")
    print(f"Partial Pressure: {a.partial_p} (atm)")
    print(f"Specific Heat at Constant Pressure: {a.Cp_mix:.6f} (kJ/kmol·K)")
    print(f"Enthalpy: {a.h:.6f} (kJ/kg)")
    print(f"Entropy: {a.s:.6f} (kJ/kg·K)")
    print(f"Gibbs Free Energy: {a.g:.6f} (kJ/kg)")
    print("\n" + "-"*50 + "\n")

dry_air_inputs = [
    "O2 N2",  # composition
    1, # pressure (atm)
    500, # temperature (C)
    21, # O2 mole fraction (%)
    79 # N2 mole fraction (%)
]

fuel_inputs = [
    "H2 CO CO2 CH4",  # composition
    1, # pressure (atm)
    600, # temperature (C)
    20, # H2 mole fraction (%)
    10, # CO mole fraction (%)
    15, # CO2 mole fraction (%)
    55 # CH4 mole fraction (%)
]



if __name__ == "__main__":
    run_case("Dry Air, O2:N2 = 21:79, 1 atm, 500 C", dry_air_inputs)
    run_case("Fuel, H2:CO:CO2:CH4 = 20:10:15:55, 1 atm, 600 C", fuel_inputs)
