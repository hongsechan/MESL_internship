from thermo_properties import ThermoProperties
from eos_deaparture import EOSDeparture
import numpy as np

class RealGasThermo(ThermoProperties):
    def __init__(
        self, comp_name, p, t_C, mole_fraction_percentage,
        eos_type="id",
    ):

        super().__init__(
            comp_name=comp_name,
            p=p,
            t_C=t_C,
            mole_fraction_percentage=mole_fraction_percentage,
        )

        eos_names = np.array(["id", "vdW", "SRK", "PR"])
        matches = np.flatnonzero(
            np.char.lower(eos_names) == eos_type.strip().lower()
        )

        if matches.size == 0:
            raise ValueError("eos_type은 id, vdW, SRK, PR 중 하나여야 합니다.")

        self.idx_eos = int(matches[0])
        self.eos_type = eos_names[self.idx_eos]

        self.h_ideal = self.h
        self.s_ideal = self.s
        self.cp_ideal = self.Cp_mix
        self.g_ideal = self.g

        self.eos = None
        self.z = 1.0
        self.h_departure = 0.0
        self.s_departure = 0.0
        self.cp_departure = 0.0

        # id는 이상기체 결과를 그대로 사용
        if self.eos_type == "id":
            return

        # EOS의 몰 기준 잔차 계산
        self.eos = EOSDeparture(
            comp_name=comp_name,
            p=p,
            t_C=t_C,
            mole_fraction_percentage=mole_fraction_percentage,
            eos_type=self.eos_type,
        )
        self.z = self.eos.z

        self.h_departure = self.eos.enthalpy_departure / self.M_mix
        self.s_departure = self.eos.entropy_departure / self.M_mix
        self.cp_departure = self.eos.cp_departure / self.M_mix

        self.h = self.h_ideal + self.h_departure
        self.s = self.s_ideal + self.s_departure
        self.Cp_mix = self.cp_ideal + self.cp_departure

        self.g = self.g_ideal + (
            self.eos.enthalpy_departure
            - self.t * self.eos.entropy_departure
        )



if __name__ == "__main__":
    for eos_type in ("id", "vdW", "SRK", "PR"):
        state = RealGasThermo(
            comp_name=["H2", "CO2", "CH4"],
            p=10.0,  # atm
            t_C=300.0,  # °C
            mole_fraction_percentage=[50.0, 30.0, 20.0],
            eos_type=eos_type,
        )

        print(f"\n===== {eos_type} =====")
        print(f"Z  = {state.z:.8f}")
        print(
            f"h  = {state.h_ideal:.6f} + {state.h_departure:.6f}"
            f" = {state.h:.6f} kJ/kg"
        )
        print(
            f"s  = {state.s_ideal:.6f} + {state.s_departure:.6f}"
            f" = {state.s:.6f} kJ/(kg·K)"
        )
        print(
            f"Cp = {state.cp_ideal:.6f} + {state.cp_departure:.6f}"
            f" = {state.Cp_mix:.6f} kJ/(kg·K)"
        )