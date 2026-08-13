import numpy as np

#WGS
WGS_reaction_coefficients = np.array([+1.0, -1.0, -1.0, +1.0, 0.0, 0.0, 0.0]) # H2, H2O, CO, CO2, O2, N2, CH4
                                                                              # CO + H2O <-> H2 + CO2 

SMR_reaction_coefficients = np.array([+3.0, -1.0, +1.0, 0.0, 0.0, 0.0, -1.0]) # H2, H2O, CO, CO2, O2, N2, CH4
                                                                              # CH4 + H2O <-> 3H2 + CO