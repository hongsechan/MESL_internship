from thermo_properties import ThermoProperties

a = ThermoProperties()

print(f"{a.mole_fraction}")
print(f"{a.M_mix} (kg/kmol)")
print(f"{a.partial_p} (atm)")
print(f"{a.Cp_mix:.6f} (kJ/kmol·K)")
print(f"{a.h:.6f} (kJ/kg)")
print(f"{a.s:.6f} (kJ/kg·K)")
print(f"{a.g:.6f} (kJ/kg)")