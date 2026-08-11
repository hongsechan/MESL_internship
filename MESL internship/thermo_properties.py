import numpy as np
import math
from nasa_poly_coeff import NasaPolyCoeff
import constants as const

class ThermoProperties(NasaPolyCoeff):
    def __init__(self, comp_name=None, p=None, t_C=None, mole_fraction_percentage=None):
        super().__init__(comp_name, p, t_C, mole_fraction_percentage)
        self.Cp_each = self.calc_Cp_each()
        self.Cp_mix = self.calc_Cp_mix()

        #================================================== Cp 매서드

        self.h_each_mole = self.calc_h_each_mole()
        self.h = self.calc_h()
        #================================================== 엔탈피 매서드

        self.p_ref = const.P_ref
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
        return const.R_universal * (coeff[:, :5]*t_Cp_calc).sum(axis=1)

    def calc_Cp_mix(self):
        mole_fraction = self.mole_fraction
        Cp_each = self.Cp_each
        M_mix = self.M_mix

        return np.sum(Cp_each*mole_fraction)/M_mix

#==========================================================Cp 계산
    def h_abs_at(self, T):
            coeff = self.coeff_at(T)
            T_h_calc = np.array([T, T*T/2, T*T*T/3, T*T*T*T/4, T*T*T*T*T/5, 1])  
            return const.R_universal*(coeff[:, :6]*T_h_calc).sum(axis=1)
                              

    def calc_h_each_mole(self):
        return self.h_abs_at(self.t)

    def calc_h(self):
        h_each_mole = self.h_each_mole
        mole_fraction = self.mole_fraction
        M_mix = self.M_mix
        return np.sum(h_each_mole*mole_fraction)/M_mix

    def h_at(self, T):
        h_each_mole = self.h_abs_at(T)    
        mole_fraction = self.mole_fraction
        M_mix = self.M_mix
        return np.sum(h_each_mole*mole_fraction)/M_mix
#======================================================================== 엔탈피 계산
    def s_abs_at(self, T):
        coeff = self.coeff_at(T)
        T_s_calc = np.array([math.log(T), T, T*T/2, T*T*T/3, T*T*T*T/4])
        return const.R_universal*((coeff[:, :5]*T_s_calc).sum(axis=1) + coeff[:, 6])
    
    def calc_s_abs(self):
        return self.s_abs_at(self.t)

    def s_each_at(self, T, p):
        s_abs = self.s_abs_at(T)
        partial_p = p * self.mole_fraction
        p_ref = self.p_ref
        return s_abs - const.R_universal * np.log(partial_p/p_ref)

    def calc_s_each_mole(self):
        return self.s_each_at(self.t, self.p)

    def s_at(self, T, p):
        s_each_mole = self.s_each_at(T, p)
        mole_fraction = self.mole_fraction
        M_mix = self.M_mix
        return np.sum(s_each_mole*mole_fraction)/M_mix

    def calc_s(self):
        return self.s_at(self.t, self.p)


#======================================================= 엔트로피 계산

    def calc_mu(self):
        return self.h_each_mole - self.t * self.s_each_mole

    def g_abs_at(self, T):
        return self.h_abs_at(T)  - T * self.s_abs_at(T)

    def calc_g(self):
        return np.sum(self.mu*self.mole_fraction)











if __name__ == '__main__':
    a = ThermoProperties()
    print(a.Cp_mix)
    print(a.h)
    print(a.h_at(1500))
    print(a.s)
    print(a.mu)
    print(a.g)


