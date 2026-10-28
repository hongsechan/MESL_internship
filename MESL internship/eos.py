import numpy as np
import constants as const
from mixture_properties import MixtureProperties
from eos_data import EOSData

class EOS(MixtureProperties):
    def __init__(
        self, comp_name=None, p=None, t_C=None,
        eos_type=None, mole_fraction_percentage=None,
    ):
        super().__init__(
            comp_name=comp_name,
            p=p,
            t_C=t_C,
            mole_fraction_percentage=mole_fraction_percentage,
        )

        if eos_type is None:
            eos_type = input("Enter EOS type (vdW, SRK, PR): ")

        self.eos_type = eos_type
        self.eos = np.array(["vdW", "SRK", "PR"])
        
        self.idx_eos = np.where(eos_type == self.eos)[0][0]
        eos_coeff_data = np.array([
            [0.421875, 0.125],  # vdW
            [0.42748, 0.08664],  # SRK
            [0.45724, 0.07780]   # PR
        ])
        self.eos_coeff = eos_coeff_data[self.idx_eos]

        eos_data = EOSData()
        self.data = eos_data.data[self.idx_comp]
        self.Tc = self.data[:, 1]
        self.Pc = self.data[:, 2]     # Pa
        self.omega = self.data[:, 3]
       
        self.alpha = self.calc_alpha()
        self.a_pure = self.calc_a_pure()  # Pa·m⁶/mol²
        self.b_pure = self.calc_b_pure()  # m³/mol

        self.kij = EOSData().kij_data[self.idx_eos][np.ix_(self.idx_comp, self.idx_comp)]

        self.aij = np.sqrt(self.a_pure[:, None] * self.a_pure[None, :]) * (1.0 - self.kij)

        self.a = self.mole_fraction @ self.aij @ self.mole_fraction
        self.b = self.mole_fraction @ self.b_pure
    
    def calc_alpha(self):
        if self.eos_type == "vdW":
            return np.ones_like(self.Tc)

        omega_array = np.stack([
            np.ones_like(self.omega),
            self.omega,
            self.omega*self.omega,
        ])

        if self.eos_type == "SRK":
            m_array = np.array([0.48, 1.574, -0.176])

        elif self.eos_type == "PR":
            m_array = np.array([0.37464, 1.54226, -0.26992])

        else:
            raise ValueError("EOS type must be vdW, SRK, or PR.")

        m = m_array @ omega_array
        return (1 + m * (1 - np.sqrt(self.t / self.Tc)))**2
  
        
    def calc_a_pure(self):
        return self.eos_coeff[0] * (const.R_universal*const.R_universal) * (self.Tc*self.Tc) *self.alpha / (self.Pc)

    def calc_b_pure(self):
        return self.eos_coeff[1] * const.R_universal * self.Tc / self.Pc




if __name__ == "__main__":
    for eos_type in ["vdW", "SRK", "PR"]:
        state = EOS(
            comp_name=["H2", "CO2", "CH4"],
            p=1.0,                         # atm
            t_C=300.0,                     # °C
            eos_type=eos_type,
            mole_fraction_percentage=[50, 30, 20],
        )

        x = state.mole_fraction
        n = len(x)

        # 반복문으로 혼합계수 계산 → NumPy 계산과 비교
        a_check = 0.0
        b_check = 0.0

        for i in range(n):
            b_check += x[i] * state.b_pure[i]

            for j in range(n):
                a_check += (
                    x[i] * x[j]
                    * np.sqrt(state.a_pure[i] * state.a_pure[j])
                    * (1.0 - state.kij[i, j])
                )

        np.testing.assert_allclose(state.a, a_check, rtol=1e-12, atol=0)
        np.testing.assert_allclose(state.b, b_check, rtol=1e-12, atol=0)

        print(f"\n===== {eos_type} =====")
        print(f"T = {state.t:.2f} K, P = {state.p:.2f} atm")
        print(f"{'Species':<10} {'x':>8} {'a_pure':>15} {'b_pure':>15}")

        for name, xi, ai, bi in zip(
            state.comp, x, state.a_pure, state.b_pure
        ):
            print(f"{name:<10} {xi:8.4f} {ai:15.7e} {bi:15.7e}")

        print("kij:")
        print(state.kij)
        print(f"a_mix = {state.a:.8e} Pa·m⁶/mol²")
        print(f"b_mix = {state.b:.8e} m³/mol")

