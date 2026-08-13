import numpy as np
from thermo_properties import ThermoProperties
import constants as const
import reaction as rxn

class WGS_SMR_Reformer():
    def __init__(self, mass_flow_CO = None, mass_flow_H2O=None, mass_flow_CH4=None, stream_CO=None, stream_H2O=None, stream_CH4=None):
        if mass_flow_CO is None:
            self.mass_flow_CO = float(input("Enter the mass flow rate of CO (kg/s): "))
        else:
            self.mass_flow_CO = mass_flow_CO    

        if mass_flow_CH4 is None:
            self.mass_flow_CH4 = float(input("Enter the mass flow rate of CH4 (kg/s): "))
        else:
            self.mass_flow_CH4 = mass_flow_CH4

        if mass_flow_H2O is None:
            self.mass_flow_H2O = float(input("Enter the mass flow rate of H2O (kg/s): "))    
        else:
            self.mass_flow_H2O = mass_flow_H2O

        if stream_CO is None:
            self.stream_CO = ThermoProperties()
        else:
            self.stream_CO = stream_CO
            
        if stream_CH4 is None:
            self.stream_CH4 = ThermoProperties()
        else:
            self.stream_CH4 = stream_CH4

        if stream_H2O is None:
            self.stream_H2O = ThermoProperties()
        else:
            self.stream_H2O = stream_H2O

        self.mole_inletflow_CO = self.mass_flow_CO / self.stream_CO.M_mix
        self.mole_inletflow_H2O = self.mass_flow_H2O / self.stream_H2O.M_mix
        self.mole_inletflow_CH4 = self.mass_flow_CH4 / self.stream_CH4.M_mix

        self.p_out = min(self.stream_CH4.p, self.stream_H2O.p, self.stream_CO.p)

        self.inlet_mole_flows = np.array([
            0, 
            0,
            self.mole_inletflow_CO,
            self.mole_inletflow_H2O, 
            self.mole_inletflow_CH4,
            ])

        self.stream_basis = ThermoProperties(comp_name=["H2", "CO2","CO", "H2O", "CH4"],
                                            p=self.p_out,
                                            t_C=25,
                                            mole_fraction_percentage=[20, 20, 20, 20, 20])    

        self.WGS_reaction_coefficients = rxn.WGS_reaction_coefficients[self.stream_basis.idx_comp] 
        self.SMR_reaction_coefficients = rxn.SMR_reaction_coefficients[self.stream_basis.idx_comp] 

        self.inlet_energy = (self.mass_flow_CO * self.stream_CO.h) + (self.mass_flow_H2O * self.stream_H2O.h) + (self.mass_flow_CH4 * self.stream_CH4.h) 


        # 평형 출구 온도와 반응진행도
        self.T_out, self.x, self.y = self.solve_equilibrium()
        self.T_out_C = self.T_out - const.kelvin_offset

        # 출구 몰유량과 몰분율
        self.outlet_mole_flows = self.calculate_outlet_mole_flows(self.x, self.y)
        self.outlet_total_moles = np.sum(self.outlet_mole_flows)
        self.outlet_mole_fraction = self.outlet_mole_flows / self.outlet_total_moles

        #출구 질량유량
        self.outlet_mass_flows = self.outlet_mole_flows * self.stream_basis.M[self.stream_basis.idx_comp]
        self.outlet_total_mass_flow = np.sum(self.outlet_mass_flows)


    # 출구 물질 몰
    def calculate_outlet_mole_flows(self, x, y):
        outlet_mole_flows = self.inlet_mole_flows + x * self.WGS_reaction_coefficients + y * self.SMR_reaction_coefficients
        return outlet_mole_flows


    # 검사함수
    def validate_state(self, T, x, y):  # x: WGS reaction extent, y: SMR reaction extent
        if T <= 0:
            raise ValueError("Temperature must be greater than 0 K.")
        
        mole_outflow = self.calculate_outlet_mole_flows(x, y)

        if np.any(mole_outflow <= 0):
            raise ValueError(
            "All outlet mole flows must be greater than zero."
        )
    

    #에너지 계산
    def outlet_energy(self, T, x, y):
        mole_outflow = self.calculate_outlet_mole_flows(x, y)
        outlet_energy = np.sum(mole_outflow*self.stream_basis.h_abs_at(T))
        return outlet_energy

    def energy_residual(self, T, x, y):
        return self.outlet_energy(T, x, y) - self.inlet_energy     

    # 평형 계산

    def log_equilibrium_constant(self, rxn_coeff, T): 
        delta_G = np.sum(self.stream_basis.g_abs_at(T)*rxn_coeff)
        ln_K = -delta_G/(const.R_universal*T)
        return ln_K

    def log_reaction_quotient(self, rxn_coeff, x, y):
        mole_outflow = self.calculate_outlet_mole_flows(x, y)
        activity_mask = rxn_coeff !=0   
        total_moles = np.sum(mole_outflow)
        p_term = self.p_out/const.P_ref

        react_coeff = rxn_coeff [activity_mask]
        mole_react = mole_outflow [activity_mask]

        mole_fraction = mole_react / total_moles
        activity = mole_fraction * p_term

        ln_Q = np.sum(react_coeff*np.log(activity))
        return ln_Q

    def calculate_equilibrium_residual(self, rxn_coeff, T, x, y):
        ln_K = self.log_equilibrium_constant(rxn_coeff, T)
        ln_Q = self.log_reaction_quotient(rxn_coeff, x, y)
        return ln_Q -ln_K
    
    # 연관 잔차 정리
    def coupled_residual(self, T, x, y):

        WGS_equilibrium_residual = self.calculate_equilibrium_residual(self.WGS_reaction_coefficients, T, x, y)
        SMR_equilibrium_residual = self.calculate_equilibrium_residual(self.SMR_reaction_coefficients, T, x, y)
        energy_scale = max(abs(self.inlet_energy), 1.0)
        return np.array([self.energy_residual(T, x, y)/energy_scale, WGS_equilibrium_residual, SMR_equilibrium_residual])

    #수치적 야코비안 계산
    def numerical_jacobian(self, T, x, y):
        dT = max(T*1e-4, 1e-3)

        x_lower = max(0.0, -3.0*y)
        x_upper = min(self.mole_inletflow_H2O - y, self.mole_inletflow_CO + y)
        dx_distance = min(x-x_lower, x_upper-x)
        dx_base = max((x_upper-x_lower)*1e-5, 1e-10)
        dx = min(dx_base, dx_distance*0.5)
        
        y_lower = max( - x/3.0, x - self.mole_inletflow_CO)
        y_upper = min(self.mole_inletflow_H2O - x, self.mole_inletflow_CH4)
        dy_distance = min(y - y_lower, y_upper - y)
        dy_base = max((y_upper-y_lower)*1e-5, 1e-10)
        dy = min(dy_base, dy_distance*0.5)

        df_dT = (self.coupled_residual(T + dT, x, y) - self.coupled_residual(T - dT, x, y)) / (2 * dT)
        df_dx = (self.coupled_residual(T, x + dx, y) - self.coupled_residual(T, x - dx, y)) / (2 * dx)
        df_dy = (self.coupled_residual(T, x, y + dy) - self.coupled_residual(T, x,  y - dy)) / (2 * dy)
        jacobian_matrix = np.column_stack((df_dT, df_dx, df_dy))
        return jacobian_matrix

    # 뉴텁-랩슨 방법
    def solve_equilibrium(self, tol = 1e-8, max_iter = 100):
        T_guess = (self.stream_CO.t + self.stream_H2O.t +self.stream_CH4.t)/3
        y_guess = min(self.mole_inletflow_H2O, self.mole_inletflow_CH4)*0.5
        x_guess = min(self.mole_inletflow_CO + y_guess, self.mole_inletflow_H2O - y_guess)*0.5

        variables = np.array([T_guess, x_guess, y_guess])

        for _ in range(max_iter):
            T, x, y = variables
            self.validate_state(T, x, y)

            residuals = self.coupled_residual(T, x, y)
            residual_norm = np.linalg.norm(residuals)    
             
            if residual_norm < tol:
                return T, x, y

            else : 
                jacobian = self.numerical_jacobian(T, x, y)
                try:
                    delta = np.linalg.solve(jacobian, -residuals)
                except np.linalg.LinAlgError: 
                    raise ValueError("Jacobian is singular. Cannot proceed with Newton-Raphson method.")

            alpha = 1.0
            accept = False

            while alpha > 1e-10:
                test_variables = variables + alpha * delta    

                try:    
                    self.validate_state(test_variables[0], test_variables[1], test_variables[2])
                except ValueError:
                    alpha *=0.5
                    continue

                test_residuals = self.coupled_residual(test_variables[0], test_variables[1], test_variables[2])
                test_residual_norm= np.linalg.norm(test_residuals)

                if test_residual_norm < residual_norm:
                    accept = True
                    break
                else:
                    alpha *=0.5 

            if not accept :
                raise  ValueError("Failed to find an acceptable step size.")

            variables = test_variables   

        raise ValueError("Newton-Raphson method did not converge within the maximum number of iterations.")

                    
            
