import numpy as np
from thermo_properties import ThermoProperties
import constants as const
import reaction as rxn

class WGSReformer():
    def __init__(self, mass_flow_CO = None, mass_flow_H2O = None, stream_CO = None, stream_H2O = None):
        if mass_flow_CO is None:
            self.mass_flow_CO = float(input("Enter the mass flow rate of CO (kg/s): "))
        else:
            self.mass_flow_CO = mass_flow_CO

        if mass_flow_H2O is None:
            self.mass_flow_H2O = float(input("Enter the mass flow rate of H2O (kg/s): "))    
        else:
            self.mass_flow_H2O = mass_flow_H2O

        if stream_CO is None:
            self.stream_CO = ThermoProperties()
        else:
            self.stream_CO = stream_CO

        if stream_H2O is None:
            self.stream_H2O = ThermoProperties()
        else:
            self.stream_H2O = stream_H2O


        self.mole_inletflow_CO = self.mass_flow_CO / self.stream_CO.M_mix
        self.mole_inletflow_H2O = self.mass_flow_H2O / self.stream_H2O.M_mix
        self.p_out = min(self.stream_CO.p, self.stream_H2O.p)

        self.inlet_mole_flows = np.array([
            0, 
            0, 
            self.mole_inletflow_CO, 
            self.mole_inletflow_H2O])

        self.stream_basis = ThermoProperties(comp_name=["H2", "CO2", "CO", "H2O"],
                                            p=self.p_out,
                                            t_C = 25,
                                            mole_fraction_percentage=[25, 25, 25, 25])    
            
        self.reaction_coefficients = rxn.WGS_reaction_coefficients[self.stream_basis.idx_comp]
        self.inlet_energy = (self.mass_flow_CO * self.stream_CO.h) + (self.mass_flow_H2O * self.stream_H2O.h)

        # 평형 출구 온도와 반응진행도
        self.T_out, self.x = self.solve_equilibrium()
        self.T_out_C = self.T_out - const.kelvin_offset
        # 출구 몰유량과 몰분율
        self.outlet_mole_flows = self.calculate_outlet_mole_flows(self.x)
        self.outlet_total_moles = np.sum(self.outlet_mole_flows)
        self.outlet_mole_fraction = self.outlet_mole_flows / self.outlet_total_moles

        #출구 질량유량
        self.outlet_mass_flows = self.outlet_mole_flows * self.stream_basis.M[self.stream_basis.idx_comp]
        self.outlet_total_mass_flow = np.sum(self.outlet_mass_flows)

    # 검사함수
    def validate_state(self, T, x):
        if T <= 0:
            raise ValueError("Temperature must be greater than 0 K.")

        x_max = min(self.mole_inletflow_CO, self.mole_inletflow_H2O)
        if not 0 < x < x_max:
            raise ValueError(f"Extent of reaction x must be between 0 and {x_max} kmol/s.")

    # 출구 물질 몰
    def calculate_outlet_mole_flows(self, x):   
        return self.inlet_mole_flows + x*self.reaction_coefficients   
    
    # 에너지 계산
    def outlet_energy(self, T, x):
        mole_outflow = self.calculate_outlet_mole_flows(x)
        outlet_energy = np.sum(mole_outflow * self.stream_basis.h_abs_at(T))
        return outlet_energy

    def energy_residual(self, T, x):
        return self.outlet_energy(T, x) - self.inlet_energy
    
    # 평형 계산    
    def equilibrium_constant(self, T):
        g_abs_at_T = self.stream_basis.g_abs_at(T)

        delta_G = np.sum(self.reaction_coefficients * g_abs_at_T)
        K_eq = np.exp(-delta_G / (const.R_universal * T))
        return K_eq
    
    def reaction_quotient(self, x):
        mole_outflow = self.calculate_outlet_mole_flows(x)
        activity_mask = self.reaction_coefficients != 0

        if np.any(mole_outflow[activity_mask] <= 0):
            raise ValueError("Mole outflow for a reactant or product is non-positive, cannot compute reaction quotient.")
        
        total_moles = np.sum(mole_outflow)
        if total_moles <= 0:
            raise ValueError("Total moles in the outlet flow is non-positive, cannot compute mole fractions.")
        
        mole_fraction = mole_outflow / total_moles
        p_term = self.p_out / const.P_ref

        activity = mole_fraction * p_term 
        Q = np.prod(activity[activity_mask] ** self.reaction_coefficients[activity_mask]) 
        return Q

    def equilibrium_residual(self, T, x):
        K_eq = self.equilibrium_constant(T)
        Q = self.reaction_quotient(x)
        return np.log(Q) - np.log(K_eq)

    # 연관 잔차 계산
    def coupled_residual(self, T, x):
        self.validate_state(T, x)

        energy_scale = max(abs(self.inlet_energy), 1.0) 

        return np.array([self.energy_residual(T, x) / energy_scale, self.equilibrium_residual(T, x) ])
    
    # 수치적 야코비안 계산
    def numerical_jacobian(self, T, x):
        dT = max(T*1e-4, 1e-3)

        x_max = min(self.mole_inletflow_CO, self.mole_inletflow_H2O)
        dx_base = max(x_max*1e-5, 1e-10)
        distance_to_bounds = min(x, x_max - x)
        dx = min(dx_base, distance_to_bounds * 0.5)

        df_dT = (self.coupled_residual(T + dT, x) - self.coupled_residual(T - dT, x)) / (2 * dT)
        df_dx = (self.coupled_residual(T, x + dx) - self.coupled_residual(T, x - dx)) / (2 * dx)

        jacobian_matrix = np.column_stack((df_dT, df_dx))
        return jacobian_matrix

    # 뉴턴-랩슨 방법을 이용한 연속 방정식 풀이    
    def solve_equilibrium(self, T_guess=None, x_guess=None, tol=1e-8, max_iter=100):
        if T_guess is None:
            T_guess = (self.stream_CO.t + self.stream_H2O.t) / 2

        if x_guess is None:
            x_guess = min(self.mole_inletflow_CO, self.mole_inletflow_H2O)*0.5
        
        variables = np.array([T_guess, x_guess])

        for _ in range(max_iter):
            T, x = variables
            residuals = self.coupled_residual(T, x)
            residual_norm = np.linalg.norm(residuals)

            if residual_norm < tol:
                return T, x

            else:
                jacobian = self.numerical_jacobian(T, x)

                try:
                    delta = np.linalg.solve(jacobian, -residuals)
                except np.linalg.LinAlgError:
                    raise ValueError("Jacobian is singular. Cannot proceed with Newton-Raphson method.")


                alpha = 1.0
                acceptable_step = False

                while alpha > 1e-10:
                    test_variables = variables + alpha * delta

                    test_T = test_variables[0]
                    test_x = test_variables[1]

                    try:
                        self.validate_state(test_T, test_x)
                    except ValueError:
                        alpha *= 0.5
                        continue

                    test_residuals = self.coupled_residual(test_T, test_x)
                    test_residual_norm = np.linalg.norm(test_residuals)

                    if test_residual_norm < residual_norm:
                        acceptable_step = True
                        break

                    else:
                        alpha *= 0.5

                if not acceptable_step:
                    raise ValueError("Failed to find an acceptable step size.")

                variables = test_variables


        raise ValueError("Newton-Raphson method did not converge within the maximum number of iterations.")
    

