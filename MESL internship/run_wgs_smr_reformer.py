from wgs_smr_reformer import WGS_SMR_Reformer
from thermo_properties import ThermoProperties

def print_result(case):
    print("=" * 60)
    print(f"Equilibrium outlet temperature: {case.T_out_C:.2f} °C")
    print(f"Reaction progress of WGS: {case.x:.6f}")
    print(f"Reaction progress SMR: {case.y:.6f}")
    print(f"Outlet mole flows: ")
    for name, z in zip(case.stream_basis.comp_name,
                       case.outlet_mole_flows):
        print(f"  {name}: {z:.6f} kmol/s")    

    print("Mole fraction of outlet stream:")
    for name, z in zip(case.stream_basis.comp_name,
                       case.outlet_mole_fraction):
        print(f"  {name}: {z:.6f} ")

    print("Outlet mass flows:")
    for name, z in zip(case.stream_basis.comp_name,
                           case.outlet_mass_flows):
        print(f"  {name}: {z:.6f} kg/s")

    print(f"Outlet total mass flow: {case.outlet_total_mass_flow:.4f} kg/s")



case1 = WGS_SMR_Reformer(
    mass_flow_CO=10,
    mass_flow_H2O=5,
    mass_flow_CH4=5,
    stream_CO=ThermoProperties(
        comp_name=["CO"],
        p=1,
        t_C=300,
        mole_fraction_percentage=[100]
    ),
    stream_H2O=ThermoProperties(
        comp_name=["H2O"],
        p=1,
        t_C=300,
        mole_fraction_percentage=[100]
    ),
        stream_CH4=ThermoProperties(
        comp_name=["CH4"],
        p=1,
        t_C=300,
        mole_fraction_percentage=[100]
    ),
)
print_result(case1)    