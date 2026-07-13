from thermo_properties import ThermoProperties
from compression import Compression
from pump import Pump

def print_result(title, case):
    print("=" * 60)
    print(title)
    print(f"Mass flow rate: {case.mass_flow:.2f} kg/s")
    print(f"Molar mass of the outlet stream: {case .M_out:.4f} kg/kmol")
    print(f"Outlet mole fraction: ")
    for name, y in zip(case.stream.element[case.stream.idx_comp], case.outlet_mole_fraction):
        print(f"  {name}: {y:.6f}")
    print(f"Outlet pressure: {case.p_out:.2f} atm")

    if hasattr(case, "Ts_out"):
        print(f"Outlet temperature (isentropic): {case.Ts_out_C:.2f} °C")
    
    print(f"Outlet temperature (actual): {case.T_out_C:.2f} °C")

    if hasattr(case, "hs_out"):
        print(f"Outlet enthalpy (isentropic): {case.hs_out:.2f} kJ/kg")

    print(f"Outlet enthalpy (actual): {case.h_out:.2f} kJ/kg")
    print(f"Power input: {case.W_dot:.2f} kW")


air = ThermoProperties(
    comp_name=["O2", "N2"],
    p=1,    
    t_C=25,
    mole_fraction_percentage=[21, 79]
)

case1 = Compression(
    stream=air,
    p_out=10,
    mass_flow=10,
    efficiency=1.0
)
   
print_result("Case 1: Compression of Air", case1)


water_liquid = ThermoProperties(
    comp_name=["H2O"],
    p=1,
    t_C=25,
    mole_fraction_percentage=[100]
)

case2 = Pump(
    stream=water_liquid,
    p_out=10,
    mass_flow=10,
    efficiency=0.85
)

print_result("Case 2: Pumping of Water", case2)


fuel = ThermoProperties(
    comp_name=["CH4", "H2O"],
    p=1,
    t_C=150,
    mole_fraction_percentage=[25, 75]
)

case3 = Compression(
    stream=fuel,
    p_out=2,
    mass_flow=10,  
    efficiency=0.85
)
print_result("Case 3: Compression of Fuel", case3)