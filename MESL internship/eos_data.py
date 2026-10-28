import numpy as np


class EOSData:
    def __init__(self):
        # Columns: M [kg/kmol], Tc [K], Pc [kPa], omega [-].

        self.data = np.array([
            [2.016,  33.145,    1296.4, -0.2190],  # H2 (normal)
            [18.02, 647.096,   22064.0,  0.3443],  # H2O
            [28.01, 132.86,     3494.0,  0.0497],  # CO
            [44.01, 304.1282,   7377.3,  0.2239],  # CO2
            [32.00, 154.581,    5043.0,  0.0222],  # O2
            [28.01, 126.192,    3395.8,  0.0372],  # N2
            [16.04, 190.564,    4599.2,  0.0114],  # CH4: p.1092의 채택값 
        ], dtype=float)

        # Binary interaction parameters (피팅하지 않은 초기 가정값).
        # 첫 번째 축: vdW, SRK, PR
        # 두 번째·세 번째 축: H2, H2O, CO, CO2, O2, N2, CH4
        self.kij_data = np.zeros((3, 7, 7))
