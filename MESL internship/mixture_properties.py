import numpy as np
from base import BasicInfo

class MixtureProperties(BasicInfo):
    def __init__(self):
        super().__init__()
        self.mole_fraction = self.calc_mole_fraction()
#====================================================================== 몰분율 매서드
        self.M = np.array([
            2.016,  # H2 몰질량
            18.02,  # H2O 몰질량
            28.01,  # CO 몰질량
            44.01,  # CO2 몰질량
            32.00,  # O2 몰질량
            28.01,  # N2 몰질량
            16.04  # CH4 몰질량
        ])
        self.idx_comp = np.array([np.where(i == self.element)[0][0] for i in self.comp])
        self.M_mix = self.molar_mass()

#====================================================================== 몰질량 매서드

        self.partial_p = self.p * self.mole_fraction

#====================================================================== 분압 매서드

    def calc_mole_fraction(self):
        mole_fraction_list = []
        for i in self.comp:
            y = float(input(f"Enter {i} mole fraction(%): "))
            mole_fraction_list.append(y*0.01)

        if sum(mole_fraction_list) != 1:
            print("Error: The mole percentages should add up to 100%.")
            exit()

        return np.array(mole_fraction_list)

#==========================================================================몰분율 계산

    def molar_mass(self):
        M = self.M
        idx_comp = self.idx_comp
        mole_fraction = self.mole_fraction

        return np.sum(mole_fraction * M[idx_comp])

#==========================================================================몰질량 계산


if __name__ == '__main__':
    a = MixtureProperties()
    print(a.M_mix)