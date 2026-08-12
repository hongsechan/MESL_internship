from thermo_properties import ThermoProperties
from heatexchanger import HeatExchanger


def print_result(case):
    print("=" * 60)
    print(f"Hot fluid outlet temperature: {case.hot_flow_out_temp_C:.2f} °C")
    print(f"Hot fluid outlet enthalpy: {case.h_out_hot:.2f} kJ/kg")
    print(f"Cold fluid outlet enthalpy: {case.h_out_cold:.2f} kJ/kg")
    print(f"Heat duty: {case.heat_duty:.2f} kW")


hot_air = ThermoProperties(
    comp_name=["O2", "N2"],
    p=1,
    t_C=900,
    mole_fraction_percentage=[21, 79]
)

cold_fuel = ThermoProperties(
    comp_name=["CH4", "H2O"],
    p=1,
    t_C=25,
    mole_fraction_percentage=[25, 75]
)


case = HeatExchanger(
    hot_fluid=hot_air,
    cold_fluid=cold_fuel,
    hot_fluid_massflow=100,
    cold_fluid_massflow=20,
    cold_flow_target_temp=500
)

print_result(case)