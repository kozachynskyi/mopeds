import copy

import casadi as ca
import matplotlib.cm as cm
import numpy as np

import par_est

def initialize_problem():
    P_i6_Sol10 = 3.7274E-4
    P_i6_Sol11 = -4.1033E-5
    P_i6_Sol12 = -9.9645E-6
    P_i6_Sol13 = -3.8368E-5
    P_i6_Sol14 = -6.9782E-6
    P_i6_Sol15 = -8.2558E-5
    P_i6_Sol1 = -6.4909E-5
    P_i6_Sol2 = 1.1885E-5
    P_i6_Sol3 = 0.0010631
    P_i6_Sol4 = -0.027378
    P_i6_Sol5 = 1.7599E-4
    P_i6_Sol6 = 0.17476
    P_i6_Sol7 = 9.2954E-4
    P_i6_Sol8 = 2.8881E-7
    P_i6_Sol9 = 2.9467E-4
    P_i7_Sol10 = 5.3043E-4
    P_i7_Sol11 = -7.299E-6
    P_i7_Sol12 = -1.4868E-5
    P_i7_Sol13 = -3.0261E-5
    P_i7_Sol14 = -1.2455E-5
    P_i7_Sol15 = -1.1598E-4
    P_i7_Sol1 = -1.7718E-4
    P_i7_Sol2 = 1.7692E-5
    P_i7_Sol3 = 0.0016934
    P_i7_Sol4 = -0.047302
    P_i7_Sol5 = 4.3746E-4
    P_i7_Sol6 = 0.28638
    P_i7_Sol7 = 0.001592
    P_i7_Sol8 = -1.7107E-7
    P_i7_Sol9 = 6.5328E-4
    R = 8.314

    variable_list = par_est.VariableList()
    
    # fmt: off
    variable_list.add_variable(par_est.State_variable("c_i1", 3.0))
    variable_list.add_variable(par_est.State_variable("c_i2", 10.0))
    variable_list.add_variable(par_est.State_variable("c_i3", 0.0))
    variable_list.add_variable(par_est.State_variable("c_i4", 0.0))
    variable_list.add_variable(par_est.State_variable("c_i5", 0.0))

    variable_list.add_variable(par_est.Control_variable("T", 5.0, 4.0, 6.0))
    variable_list.add_variable(par_est.Control_variable("V_Reactor", 0.0499907159751542))
    variable_list.add_variable(par_est.Control_variable("c_cat", 0.2560475430350250))
    variable_list.add_variable(par_est.Control_variable("n_Surfactant", 0.00682261208576998))
    variable_list.add_variable(par_est.Control_variable("n_Water", 1.11016746876266))
    variable_list.add_variable(par_est.Control_variable("p_Reactor", 5.0, 4.0, 6.0))

    variable_list.add_variable(par_est.Parameter_variable("greek_DeltaG_r1")
    variable_list.add_variable(par_est.Parameter_variable("K_cat1")
    variable_list.add_variable(par_est.Parameter_variable("K_r1_e1")
    variable_list.add_variable(par_est.Parameter_variable("K_r1_e2")
    variable_list.add_variable(par_est.Parameter_variable("K_r3_e1")
    variable_list.add_variable(par_est.Parameter_variable("K_r3_e2")
    variable_list.add_variable(par_est.Parameter_variable("K_r3_e3")
    variable_list.add_variable(par_est.Parameter_variable("K_r5_e1")
    variable_list.add_variable(par_est.Parameter_variable("K_r5_e2")
    variable_list.add_variable(par_est.Parameter_variable("K_r5_e3")
    variable_list.add_variable(par_est.Parameter_variable("k_ref_r1")
    variable_list.add_variable(par_est.Parameter_variable("k_ref_r2")
    variable_list.add_variable(par_est.Parameter_variable("k_ref_r3")
    variable_list.add_variable(par_est.Parameter_variable("k_ref_r4")
    variable_list.add_variable(par_est.Parameter_variable("k_ref_r5")
    variable_list.add_variable(par_est.Parameter_variable("k_ref_r6")
    variable_list.add_variable(par_est.Parameter_variable("K_cat2")
    variable_list.add_variable(par_est.Parameter_variable("E_r1")
    variable_list.add_variable(par_est.Parameter_variable("E_r2")
    variable_list.add_variable(par_est.Parameter_variable("E_r3")
    variable_list.add_variable(par_est.Parameter_variable("E_r4")
    variable_list.add_variable(par_est.Parameter_variable("E_r5")
    variable_list.add_variable(par_est.Parameter_variable("E_r6")

    # fmt: on

    # AE Calculations Reactant Moles
    dae1 = n_i1 - ( c_i1 * V_Reactor )
    dae2=	n_i2 - ( c_i2 * V_Reactor )
    dae3=	n_i3 - ( c_i3 * V_Reactor ) 
    dae4=	n_i4 - ( c_i4 * V_Reactor ) 
    dae5=	n_i5 - ( c_i5 * V_Reactor ) 
    dae6=n_L  - ( ( n_i1 + n_i2 + n_i3 + n_i4 + n_i5 ) + n_Water + n_Surfactant )
	
#     % MLS Moles   
#     greek_gamma = ( ( n_Surfactant * M_Surfactant )/( ( c_i1 * V_Reactor * M_i1 + c_i2 * V_Reactor * M_i2 + c_i3 * V_Reactor * M_i3 + c_i4 * V_Reactor * M_i4 + c_i5 * V_Reactor * M_i5 ) + n_Water * M_Water + n_Surfactant * M_Surfactant ) ); 
#         greek_alpha = ( ( ( c_i1 * V_Reactor * M_i1 + c_i2 * V_Reactor * M_i2 + c_i3 * V_Reactor * M_i3 + c_i4 * V_Reactor * M_i4 + c_i5 * V_Reactor * M_i5 ) )/( ( c_i1 * V_Reactor * M_i1 + c_i2 * V_Reactor * M_i2 + c_i3 * V_Reactor * M_i3 + c_i4 * V_Reactor * M_i4 + c_i5 * V_Reactor * M_i5 ) + n_Water * M_Water ) ); 
#         X           = ( ( ( c_i3 * M_i3 + c_i5 * M_i5 ) * V_Reactor )/( ( c_i1 * V_Reactor * M_i1 + c_i2 * V_Reactor * M_i2 + c_i3 * V_Reactor * M_i3 + c_i4 * V_Reactor * M_i4 + c_i5 * V_Reactor * M_i5 ) ) ); 
    
#     % Gas Solubility
#     x_i7 = ( ( p_Reactor )/( 2.0 ) * P_i7_Sol1 + ( T - 273.15 ) * P_i7_Sol2 + greek_alpha * P_i7_Sol3 + greek_gamma * P_i7_Sol4 + X * P_i7_Sol5 + ( ( greek_gamma ) )^( 2.0 ) * P_i7_Sol6 + ( ( X ) )^( 2.0 ) * P_i7_Sol7 + ( p_Reactor )/( 2.0 ) * ( T - 273.15 ) * P_i7_Sol8 + ( p_Reactor )/( 2.0 ) * greek_alpha * P_i7_Sol9 + ( p_Reactor )/( 2.0 ) * greek_gamma * P_i7_Sol10 + ( p_Reactor )/( 2.0 ) * X * P_i7_Sol11 + ( T - 273.15 ) * greek_alpha * P_i7_Sol12 + ( T - 273.15 ) * greek_gamma * P_i7_Sol13 + ( T - 273.15 ) * X * P_i7_Sol14 + greek_alpha * X * P_i7_Sol15 ); 
# % 	x_i6 = ( 2.0 * ( ( p_Reactor )/( 2.0 ) * P_i6_Sol1 + ( T - 273.15 ) * P_i6_Sol2 + greek_alpha * P_i6_Sol3 + greek_gamma * P_i6_Sol4 + X * P_i6_Sol5 + ( ( greek_gamma ) )^( 2.0 ) * P_i6_Sol6 + ( ( X ) )^( 2.0 ) * P_i6_Sol7 + ( p_Reactor )/( 2.0 ) * ( T - 273.15 ) * P_i6_Sol8 + ( p_Reactor )/( 2.0 ) * greek_alpha * P_i6_Sol9 + ( p_Reactor )/( 2.0 ) * greek_gamma * P_i6_Sol10 + ( p_Reactor )/( 2.0 ) * X * P_i6_Sol11 + ( T - 273.15 ) * greek_alpha * P_i6_Sol12 + ( T - 273.15 ) * greek_gamma * P_i6_Sol13 + ( T - 273.15 ) * X * P_i6_Sol14 + greek_alpha * X * P_i6_Sol15 ) - x_i7 ); 
#         x_i6 = (  ( ( p_Reactor ) * P_i6_Sol1 + ( T - 273.15 ) * P_i6_Sol2 + greek_alpha * P_i6_Sol3 + greek_gamma * P_i6_Sol4 + X * P_i6_Sol5 + ( ( greek_gamma ) )^( 2.0 ) * P_i6_Sol6 + ( ( X ) )^( 2.0 ) * P_i6_Sol7 + ( p_Reactor ) * ( T - 273.15 ) * P_i6_Sol8 + ( p_Reactor ) * greek_alpha * P_i6_Sol9 + ( p_Reactor ) * greek_gamma * P_i6_Sol10 + ( p_Reactor ) * X * P_i6_Sol11 + ( T - 273.15 ) * greek_alpha * P_i6_Sol12 + ( T - 273.15 ) * greek_gamma * P_i6_Sol13 + ( T - 273.15 ) * X * P_i6_Sol14 + greek_alpha * X * P_i6_Sol15 ) - x_i7 ); 
#     c_i7 = ( ( n_L * x_i7 )/( 1.0 - x_i7 ) ) / V_Reactor; % in mol/l
#         c_i6 = ( ( n_L * x_i6 )/( 1.0 - x_i6 ) ) / V_Reactor; % in mol/l 
	
#     % Reaction Rates
#     greek_psi_cat = ( c_cat ) / ( 1.0 + K_cat1 * c_i7 + K_cat2*c_i7/c_i6); 
#     greek_DeltaG_r3 = ( (  - 126.28 + 0.13 * T + 6.8 * ( ( 10.0 ) )^(  - 6.0 ) * ( ( T ) )^( 2.0 ) ) * ( ( 10.0 ) )^( 3.0 ) ); 
#     K_eq_r3 = ( exp(  - ( greek_DeltaG_r3 )/( R * T ) ) ); 
#         K_eq_r1 = ( exp(   ( greek_DeltaG_r1 )/( R * T ) ) ); 
# %     r_r1 = n_Surfactant / V_Reactor * (1+(LM_ratio)^k_LM_exp_r1) * ( greek_psi_cat * k_ref_r1 * exp(  - ( E_r1 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * ( c_i1 - ( c_i2 )/(K_eq_r1 ) ) )./( 1.0 + K_r1_e1 * c_i1 + K_r1_e2 * c_i2 ); 
# % 	r_r2 = n_Surfactant / V_Reactor * (1+(LM_ratio)^k_LM_exp_r2) * ( greek_psi_cat * k_ref_r2 * exp(  - ( E_r2 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i2 * c_i6 ); 
# % 	r_r3 = n_Surfactant / V_Reactor * (1+(LM_ratio)^k_LM_exp_r3) * ( greek_psi_cat * k_ref_r3 * exp(  - ( E_r3 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i2 * c_i6 * ( c_i4 )/( K_eq_r3 ) ) ./ ( 1.0 + K_r3_e1 * c_i1 + K_r3_e2 * c_i4 + K_r3_e3 * c_i6 ); 
# % 	r_r4 = n_Surfactant / V_Reactor * (1+(LM_ratio)^k_LM_exp_r4) * ( greek_psi_cat * k_ref_r4 * exp(  - ( E_r4 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i2 * c_i6 * c_i7 ); 
# % 	r_r5 = n_Surfactant / V_Reactor * (1+(LM_ratio)^k_LM_exp_r5) * ( greek_psi_cat * k_ref_r5 * exp(  - ( E_r5 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i1 * c_i6 * c_i7 ) ./ ( 1.0 + K_r5_e1 * c_i1 + K_r5_e2 * c_i5 + K_r5_e3 * c_i6 ); 
# % 	r_r6 = n_Surfactant / V_Reactor * (1+(LM_ratio)^k_LM_exp_r6) * ( greek_psi_cat * k_ref_r6 * exp(  - ( E_r6 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i1 * c_i6 * c_i7 ); 
#     r_r1 = ( greek_psi_cat * k_ref_r1 * exp(  - ( E_r1 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * ( c_i1 - ( c_i2 )/(K_eq_r1 ) ) )./( 1.0 + K_r1_e1 * c_i1 + K_r1_e2 * c_i2 ); 
#         r_r2 = ( greek_psi_cat * k_ref_r2 * exp(  - ( E_r2 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i2 * c_i6 ); 
#         r_r3 = ( greek_psi_cat * k_ref_r3 * exp(  - ( E_r3 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * (c_i1 * c_i6 - ( c_i4 )/( K_eq_r3 ) ) ) ./ ( 1.0 + K_r3_e1 * c_i1 + K_r3_e2 * c_i4 + K_r3_e3 * c_i6 ); 
#         r_r4 = ( greek_psi_cat * k_ref_r4 * exp(  - ( E_r4 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i2 * c_i6 * c_i7 ); 
#         r_r5 = ( greek_psi_cat * k_ref_r5 * exp(  - ( E_r5 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i1 * c_i6 * c_i7 ) ./ ( 1.0 + K_r5_e1 * c_i1 + K_r5_e2 * c_i5 + K_r5_e3 * c_i6 ); 
#         r_r6 = ( greek_psi_cat * k_ref_r6 * exp(  - ( E_r6 )/( R ) * ( ( 1.0 )/( T ) - ( 1.0 )/( T_ref ) ) ) * c_i1 * c_i6 * c_i7 ); 

#     % Component Reaction Rates
#     r_i1 = (  - r_r1 - r_r3 - r_r5 - r_r6 ); 
#         r_i2 = ( r_r1 - r_r2 - r_r4 ); 
#         r_i3 = ( r_r4 + r_r6 ); 
#         r_i4 = ( r_r2 + r_r3 ); 
#         r_i5 = ( r_r5 ); 
#         % evaluate the function values  
#         DYDX(1) = r_i1 .* 60; % in mol/l h
#         DYDX(2) = r_i2 .* 60; % in mol/l h  
#         DYDX(3) = r_i3 .* 60; % in mol/l h 
#         DYDX(4) = r_i4 .* 60; % in mol/l h 
#         DYDX(5) = r_i5 .* 60; % in mol/l h 


#         DYDX=DYDX';



end

function [Y_INIT, PARAMS]= EQS1_initialization()
	% specify the initial values for the state variables 
	Y_INIT(1) = 2;  	% c_i1  
	Y_INIT(2) = 1;  	% c_i2  
	Y_INIT(3) = 1;  	% c_i3  
	Y_INIT(4) = 1;  	% c_i4  
	Y_INIT(5) = 1;  	% c_i5  
	Y_INIT(6) = 0.5;  	% greek_alpha  
	Y_INIT(7) = 0.08;  	% greek_gamma  
	Y_INIT(8) = 0.01;  	% greek_psi_cat  
	Y_INIT(9) = 0.0;  	% e0_K_eq_r1  
	Y_INIT(10) = 0.0;  	% e0_K_eq_r3  
	Y_INIT(11) = 0.04;  	% e0_X  
	Y_INIT(12) = 0.01;  	% e0_c_i6  
	Y_INIT(13) = 0.01;  	% e0_c_i7  
	Y_INIT(14) = 3.0;  	% e0_n_L  
	Y_INIT(15) = 1.0;  	% e0_n_i1  
	Y_INIT(16) = 1.0;  	% e0_n_i2  
	Y_INIT(17) = 1.0;  	% e0_n_i3  
	Y_INIT(18) = 1.0;  	% e0_n_i4  
	Y_INIT(19) = 1.0;  	% e0_n_i5  
	Y_INIT(20) = 0.0;  	% e0_r_i1  
	Y_INIT(21) = 0.0;  	% e0_r_i2  
	Y_INIT(22) = 0.0;  	% e0_r_i3  
	Y_INIT(23) = 0.0;  	% e0_r_i4  
	Y_INIT(24) = 0.0;  	% e0_r_i5  
	Y_INIT(25) = 0.0;  	% e0_r_r1  
	Y_INIT(26) = 0.0;  	% e0_r_r2  
	Y_INIT(27) = 0.0;  	% e0_r_r3  
	Y_INIT(28) = 0.0;  	% e0_r_r4  
	Y_INIT(29) = 0.0;  	% e0_r_r5  
	Y_INIT(30) = 0.0;  	% e0_r_r6  
	Y_INIT(31) = 0.02;  	% e0_x_i6  
	Y_INIT(32) = 0.02;  	% e0_x_i7  

	% declare parameters 

	PARAMS(65) = 1.0;  	% e0_n_Surfactant 
	PARAMS(66) = 1.0;  	% e0_n_Water 
	PARAMS(67) = 15.0;  	% e0_p_Reactor 
