from mass_solver_model2 import Model2
from thermo_properties import ThermoProperties
import constants as const

air = ThermoProperties(
    comp_name=["O2", "N2"],
    p=1,
    t_C=500,
    mole_fraction_percentage=[21, 79]
)

model2 = Model2(stream1=air, mass_flow1=10, split_ratio=0.6)

def print_model2_result(model2):
    print("Stream 2 (Mixer Outlet):")
    print(f"Mass flow rate: {model2.mass_flow2:.2f} kg/s")
    print(f"Temperature: {model2.stream2.t - const.kelvin_offset:.2f} °C")
    print(f"Pressure: {model2.stream2.p:.2f} atm")
    print(f"Enthalpy: {model2.stream2.h:.2f} kJ/kg")

    print("\nStream 3 (Splitter Outlet 1):")
    print(f"Mass flow rate: {model2.mass_flow3:.2f} kg/s")
    print(f"Temperature: {model2.stream3.t - const.kelvin_offset:.2f} °C")
    print(f"Pressure: {model2.stream3.p:.2f} atm")
    print(f"Enthalpy: {model2.stream3.h:.2f} kJ/kg")

    print("\nStream 4 (Splitter Outlet 2):")
    print(f"Mass flow rate: {model2.mass_flow4:.2f} kg/s")
    print(f"Temperature: {model2.stream4.t - const.kelvin_offset:.2f} °C")
    print(f"Pressure: {model2.stream4.p:.2f} atm")
    print(f"Enthalpy: {model2.stream4.h:.2f} kJ/kg")


print_model2_result(model2)