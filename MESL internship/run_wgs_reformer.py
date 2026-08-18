from wgs_reformer import WGSReformer
from thermo_properties import ThermoProperties

def print_result(case):
    print("=" * 60)
    print(f"Equilibrium outlet temperature: {case.T_out_C:.2f} °C")
    print(f"Reaction progress (x): {case.x:.6f}")
    print(f"Outlet mole flows: ")
    for name, y in zip(case.stream_basis.comp_name,
                       case.outlet_mole_flows):
        print(f"  {name}: {y:.6f} kmol/s")    


    print("Mole fraction of outlet stream:")
    for name, y in zip(case.stream_basis.comp_name,
                       case.outlet_mole_fraction):
        print(f"  {name}: {y:.6f} ")

    print("Outlet mass flows:")
    for name, y in zip(case.stream_basis.comp_name,
                           case.outlet_mass_flows):
        print(f"  {name}: {y:.6f} kg/s")

    print(f"Outlet total mass flow: {case.outlet_total_mass_flow:.4f} kg/s")



case1 = WGSReformer(
    mass_flow=10,
    stream=ThermoProperties(
        comp_name=["H2O", "CO"],
        p=1,
        t_C=300,
        mole_fraction_percentage=[30, 70]
    )
    )
print_result(case1)    