import numpy as np
import constants as const

from eos import EOS

class EOSDeparture(EOS):
    def __init__(self, comp_name=None, p=None, t_C=None, eos_type=None, mole_fraction_percentage=None):
        super().__init__(comp_name=comp_name, p=p, t_C=t_C, eos_type=eos_type, mole_fraction_percentage=mole_fraction_percentage)

        self.enthalpy_departure = self.calc_enthalpy_departure()
        self.entropy_departure = self.calc_entropy_departure()
        self.cp_departure = self.calc_cp_departure()

        

    def calc_enthalpy_departure(self):
        R = const.R_universal    
        term1 = R * self.t * (self.z - 1.0)
        a = self.a
        v = self.v

        if self.eos_type == "vdW":
            return term1 - (a / v)

        elif self.eos_type in ["SRK", "PR"]:
            t = self.t       
            b = self.b
            B = self.eos_B
            z = self.z
            t = self.t
            v = self.v
            da_dT = self.da_dT
            u = self.eos_u[self.idx_eos]
            w = self.eos_w[self.idx_eos] 
            d = np.sqrt(u*u - 4*w)

            term2 = (t*da_dT - a) / (b*d) * np.log(
                (2*z + B*(u + d)) / (2*z + B*(u - d))
            )
            return term1 + term2

    def calc_entropy_departure(self):
        """Molar entropy departure [kJ/(kmol K)] at fixed T, P, composition."""
        R = const.R_universal
        term1 = R * np.log(self.z - self.eos_B)

        if self.eos_type == "vdW":
            return term1
        z = self.z
        b = self.b
        B = self.eos_B
        da_dT = self.da_dT
        u = self.eos_u[self.idx_eos]
        w = self.eos_w[self.idx_eos]
        d = np.sqrt(u*u - 4.0*w)
        log_ratio = np.log(
            (2.0*z + B*(u+d))
            / (2.0*z + B*(u-d))
        )

        return term1 + da_dT * log_ratio / (b*d)

    def calc_d2a_dT2(self):
        "Mixture d²a/dT²."
        if self.eos_type == "vdW":
            return 0.0

        T = self.t
        f = 1.0 + self.m * (1.0 - np.sqrt(T / self.Tc))
        df = -self.m / (2.0 * np.sqrt(T * self.Tc))
        d2f = self.m / (4.0 * np.sqrt(self.Tc) * T**1.5)
        d2a_pure = 2.0 * self.a_c * (df*df + f*d2f)

        relative_da = self.da_pure_dT / self.a_pure
        relative_d2a = d2a_pure / self.a_pure
        d2aij = self.aij * (
            0.5 * (relative_d2a[:, None] + relative_d2a[None, :])
            - 0.25 * (relative_da[:, None] - relative_da[None, :])**2
        )
        return self.mole_fraction @ d2aij @ self.mole_fraction

    def calc_dz_dT(self):
        """Implicit temperature derivative of the selected cubic Z root at fixed P, x."""
        Z, A, B = self.z, self.eos_A, self.eos_B
        u = self.eos_u[self.idx_eos]
        w = self.eos_w[self.idx_eos]
        dA = A * (self.da_dT / self.a - 2.0 / self.t)
        dB = -B / self.t
        K2 = (u - 1.0)*B - 1.0
        K3 = A - u*B - (u - w)*B*B
        dK2 = (u - 1.0)*dB
        dK3 = dA - u*dB - 2.0*(u - w)*B*dB
        dK4 = -(dA*B + A*dB + 2.0*w*B*dB + 3.0*w*B*B*dB)
        return -(dK2*Z*Z + dK3*Z + dK4) / (3.0*Z*Z + 2.0*K2*Z + K3)

    def calc_cp_departure(self):
        """Molar isobaric heat-capacity departure [kJ/(kmol K)]."""
        T, Z, R = self.t, self.z, const.R_universal
        dZ = self.calc_dz_dT()
        term1 = R * (Z - 1.0) + R*T*dZ

        if self.eos_type == "vdW":
            # hR = RT(Z-1) - a/v; a is constant for vdW.
            return term1 + (self.a / self.v) * (dZ / Z + 1.0 / T)

        u = self.eos_u[self.idx_eos]
        w = self.eos_w[self.idx_eos]
        d = np.sqrt(u*u - 4.0*w)
        B, dB = self.eos_B, -self.eos_B / T
        D_plus = 2.0*Z + B*(u+d)
        D_minus = 2.0*Z + B*(u-d)
        log_ratio = np.log(D_plus / D_minus)
        dlog_ratio = (
            (2.0*dZ + dB*(u+d)) / D_plus
            - (2.0*dZ + dB*(u-d)) / D_minus
        )
        return term1 + (
            T*self.calc_d2a_dT2()*log_ratio
            + (T*self.da_dT - self.a)*dlog_ratio
        ) / (self.b*d)


if __name__ == "__main__":
    for eos_type in ("vdW", "SRK", "PR"):
        state = EOSDeparture(
            comp_name=["CH4"],
            p= 10.0,                         # atm
            t_C= 600,                     # °C
            eos_type=eos_type,
            mole_fraction_percentage=[100.0],
        )

        print(
            f"{eos_type}: Z={state.z:.8f}, "
            f"hR={state.enthalpy_departure:.8f} kJ/kmol, "
            f"sR={state.entropy_departure:.8f} kJ/(kmol K), "
            f"CpR={state.cp_departure:.8f} kJ/(kmol K)"
        )
        print(
            f"    mass basis: hR={state.enthalpy_departure / state.M_mix:.8f} kJ/kg, "
            f"sR={state.entropy_departure / state.M_mix:.8f} kJ/(kg K), "
            f"CpR={state.cp_departure / state.M_mix:.8f} kJ/(kg K)"
        )
