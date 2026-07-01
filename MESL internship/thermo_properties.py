import numpy as np
import math
from nasa_poly_coeff import NasaPolyCoeff

class ThermoProperties(NasaPolyCoeff):
    def __init__(self):
        super().__init__()
        self.Cp_each = self.calc_Cp_each()
        self.Cp_mix = self.calc_Cp_mix()

        #================================================== Cp 매서드
        self.t_ref = 298.15
        self.h_f = np.array([
            0,  # H2
            -241820,  # H2O
            -110530,  # CO
            -393520,  # CO2
            0,  # O2
            0,  # N2
            -74850  # CH4
        ])
        self.delta_h = self.calc_delta_h()
        self.h_each_mole = self.calc_h_each_mole()
        self.h = self.calc_h()
        #================================================== 엔탈피 매서드

        self.p_ref = 1
        self.s_abs = self.calc_s_abs()
        self.s_each_mole = self.calc_s_each_mole()
        self.s = self.calc_s()
        #================================================= 엔트로피 매서드

        self.mu = self.calc_mu()
        self.g = self.calc_g()

        #================================================= 깁스 에너지 매서드
    def calc_Cp_each(self):
        t = self.t
        t_Cp_calc = np.array([1, t, t*t, t*t*t, t*t*t*t])
        coeff = self.coeff
        return 8.314 * (coeff[:, :5]*t_Cp_calc).sum(axis=1)

    def calc_Cp_mix(self):
        mole_fraction = self.mole_fraction
        Cp_each = self.Cp_each
        M_mix = self.M_mix

        return np.sum(Cp_each*mole_fraction)/M_mix

#==========================================================Cp 계산
    def calc_delta_h(self):
        coeff = self.coeff
        t_ref = self.t_ref
        t = self.t
        t_delta_h_calc = np.array([t-t_ref, (t*t-t_ref*t_ref)/2, (t*t*t-t_ref*t_ref*t_ref)/3, (t*t*t*t-t_ref*t_ref*t_ref*t_ref)/4, (t*t*t*t*t-t_ref*t_ref*t_ref*t_ref*t_ref)/5])    
        return       8.314*(coeff[:, :5]*t_delta_h_calc).sum(axis=1)
                              

    def calc_h_each_mole(self):
        delta_h = self.delta_h
        h_f = self.h_f
        idx_comp = self.idx_comp

        return h_f[idx_comp]+delta_h

    def calc_h(self):
        h_each_mole = self.h_each_mole
        mole_fraction = self.mole_fraction
        M_mix = self.M_mix
        return np.sum(h_each_mole*mole_fraction)/M_mix
#======================================================================== 엔탈피 계산

    def calc_s_abs(self):
        coeff = self.coeff
        t= self.t
        t_s_calc = np.array([math.log(t), t, t*t/2, t*t*t/3, t*t*t*t/4])
        return 8.314*((coeff[:, :5]*t_s_calc).sum(axis=1) + coeff[:, 6])


    def calc_s_each_mole(self):
        partial_p = self.partial_p
        s_abs = self.s_abs
        p_ref = self.p_ref
        return s_abs - 8.314 * np.log(partial_p/p_ref)


    def calc_s(self):
        s_each_mole = self.s_each_mole
        mole_fraction = self.mole_fraction
        M_mix = self.M_mix
        return np.sum(s_each_mole*mole_fraction)/M_mix

#======================================================= 엔트로피 계산

    def calc_mu(self):
        return self.h_each_mole - self.t * self.s_each_mole

    def calc_g(self):
        return np.sum(self.mu*self.mole_fraction)/self.M_mix











if __name__ == '__main__':
    a = ThermoProperties()
    print(a.Cp_mix)
    print(a.h)
    print(a.s)
    print(a.mu)
    print(a.g)


