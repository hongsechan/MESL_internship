from thermo_properties import ThermoProperties
from mixing import Mixing


def print_mixing_result(title, mixing):
    print("=" * 60)
    print(title)
    print(f"Mass flow rate of stream 1: {mixing.mass_flow1:.2f} kg/s")
    print(f"Mass flow rate of stream 2: {mixing.mass_flow2:.2f} kg/s")
    print(f"Molar mass of mixed stream: {mixing.M_out:.4f} kg/kmol")

    print("Mole fraction of mixed stream:")
    for name, y in zip(mixing.stream1.element[mixing.outlet_idx_comp],
                       mixing.outlet_mole_fraction_selected):
        print(f"  {name}: {y:.6f}")

    if mixing.T_out is None:
        print("T_out could not be calculated.")
    else:
        print(f"Temperature of mixed stream: {mixing.T_out:.2f} °C")

    print(f"Pressure of mixed stream: {mixing.p_out:.2f} atm")


air = ThermoProperties(
    comp_name=["O2", "N2"],
    p=1,
    t_C=500,
    mole_fraction_percentage=[21, 79]
)

steam = ThermoProperties(
    comp_name=["H2O"],
    p=1,
    t_C=500,
    mole_fraction_percentage=[100]
)

case1 = Mixing(
    stream1=air,
    mass_flow1=10,
    stream2=steam,
    mass_flow2=5
)

print_mixing_result("Case 1: Air + Steam", case1)


methane = ThermoProperties(
    comp_name=["CH4"],
    p=1,
    t_C=500,
    mole_fraction_percentage =[100]
)

syngas = ThermoProperties(
    comp_name=["H2", "CO", "CO2"],
    p=1,
    t_C=100,
    mole_fraction_percentage   =[40, 30, 30]
)

case2 = Mixing(
    stream1=methane,
    mass_flow1=4,
    stream2=syngas,
    mass_flow2=6
)

print_mixing_result("Case 2: Methane + Syngas", case2)