from mass_solver_model1 import Model1
from thermo_properties import ThermoProperties
import constants as const

air = ThermoProperties(
    comp_name=["O2", "N2"],
    p=1,
    t_C=500,
    mole_fraction_percentage=[21, 79]
)

steam = ThermoProperties(
    comp_name=["H2O"],
    p=1,
    t_C=300,
    mole_fraction_percentage=[100]
)

model1 = Model1(stream1=air, stream2=steam, mass_flow1=10, mass_flow2=5, split_ratio=0.6)

def print_model1_result(model1):
    print("Stream 3 (Mixer Outlet):")
    print(f"Mass flow rate: {model1.mass_flow3:.2f} kg/s")
    print(f"Temperature: {model1.stream3.t - const.kelvin_offset:.2f} °C")
    print(f"Pressure: {model1.stream3.p:.2f} atm")
    print(f"Enthalpy: {model1.stream3.h:.2f} kJ/kg")

    print("\nStream 4 (Splitter Outlet 1):")
    print(f"Mass flow rate: {model1.mass_flow4:.2f} kg/s")
    print(f"Temperature: {model1.stream4.t - const.kelvin_offset:.2f} °C")
    print(f"Pressure: {model1.stream4.p:.2f} atm")
    print(f"Enthalpy: {model1.stream4.h:.2f} kJ/kg")

    print("\nStream 5 (Splitter Outlet 2):")
    print(f"Mass flow rate: {model1.mass_flow5:.2f} kg/s")
    print(f"Temperature: {model1.stream5.t - const.kelvin_offset:.2f} °C")
    print(f"Pressure: {model1.stream5.p:.2f} atm")
    print(f"Enthalpy: {model1.stream5.h:.2f} kJ/kg")


print_model1_result(model1)