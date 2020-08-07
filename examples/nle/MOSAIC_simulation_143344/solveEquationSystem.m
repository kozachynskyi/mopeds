%*********************************************************
% The namespaces have been normalized. The following
% table shows the attribuation. 
% Normalized Name --> Full Name ---> User-defined Name
% =================================== 
% e0 --> e[0]143337 --> 
%*********************************************************

%*********************************************************
% The variables are named according to the notation
% provided in the Mosaic model.
% 
% The variable names can be read as follows:
% ==========================================
% 	e0_x_Init_x#
% 		x: states
% 		Subscript
% 			Init: Initial conditions
% 		Indices
% 			x: 1 = c_A, 2 = c_B, 3 = c_C, 4 = T
% 	 
% 	e0_c_C#
% 		c: concentration in mol/m^3
% 		Subscript
% 			C0: Component C, inlet
% 	 
% 	e0_greek_rho
% 		&rho;: density in kg/m^3
% 	 
% 	e0_H_A
% 		H: reaction enthalpy in J/mol
% 		Subscript
% 			A: Component A
% 	 
% 	e0_H_B
% 		H: reaction enthalpy in J/mol
% 		Subscript
% 			B: Component B
% 	 
% 	e0_T_Cool
% 		T: temperature in K
% 		Subscript
% 			Cool: Cooling/Heating
% 	 
% 	e0_T_In
% 		T: temperature in K
% 		Subscript
% 			In: Inlet
% 	 
% 	e0_UA
% 		UA: heat transfer coefficient times surface area in W/K
% 	 
% 	e0_cp
% 		cp: isobaric heat capacity in J/kg/K
% 	 
% 	e0_greek_Deltat
% 		&Delta;t: time intervall
% 	 
% 	e0_E_A
% 		E: activation energy in J/mol
% 		Subscript
% 			A: Component A
% 	 
% 	e0_E_B
% 		E: activation energy in J/mol
% 		Subscript
% 			B: Component B
% 	 
% 	e0_F
% 		F: flow rate in m^3/s
% 	 
% 	e0_R
% 		R: gas constant in J/mol/K
% 	 
% 	e0_V
% 		V: volume (reactor) in m^3
% 	 
% 	e0_c_A#
% 		c: concentration in mol/m^3
% 		Subscript
% 			A0: Component A, inlet
% 	 
% 	e0_dldgreek_tau_i#_j#
% 		dld&tau;: first order derivative of the lagrangian polynomials with respect to tau
% 		Indices
% 			i: Collocation coefficient at time tau_i
% 			j: value at time tau_j
% 	 
% 	e0_k_A
% 		k: preexponential factor in 1/s
% 		Subscript
% 			A: Component A
% 	 
% 	e0_k_B
% 		k: preexponential factor in 1/s
% 		Subscript
% 			B: Component B
% 	 
% 	e0_c_B#
% 		c: concentration in mol/m^3
% 		Subscript
% 			B0: Component B, inlet
% 	 
% 	e0_c_A_fe#_i#
% 		c: concentration in mol/m^3
% 		Subscript
% 			A: Component A
% 		Indices
% 			fe: number of finite elements
% 			i: Collocation coefficient at time tau_i
% 	 
% 	e0_c_C_fe#_i#
% 		c: concentration in mol/m^3
% 		Subscript
% 			C: Component C
% 		Indices
% 			fe: number of finite elements
% 			i: Collocation coefficient at time tau_i
% 	 
% 	e0_T_fe#_i#
% 		T: temperature in K
% 		Indices
% 			fe: number of finite elements
% 			i: Collocation coefficient at time tau_i
% 	 
% 	e0_c_B_fe#_i#
% 		c: concentration in mol/m^3
% 		Subscript
% 			B: Component B
% 		Indices
% 			fe: number of finite elements
% 			i: Collocation coefficient at time tau_i
% 	 
%*********************************************************

function[ROOTS]=solveEquationSystem()

	% load variable init values 
	X_ITER(1) = 8.0;  	% e0_c_A_fe1_i0 
	X_ITER(2) = 8.0;  	% e0_c_A_fe1_i3 
	X_ITER(3) = 8.0;  	% e0_c_A_fe2_i0 
	X_ITER(4) = 8.0;  	% e0_c_A_fe2_i3 
	X_ITER(5) = 8.0;  	% e0_c_A_fe3_i0 
	X_ITER(6) = 8.0;  	% e0_c_A_fe3_i3 
	X_ITER(7) = 8.0;  	% e0_c_A_fe4_i0 
	X_ITER(8) = 8.0;  	% e0_c_A_fe4_i3 
	X_ITER(9) = 8.0;  	% e0_c_A_fe5_i0 
	X_ITER(10) = 8.0;  	% e0_c_A_fe5_i3 
	X_ITER(11) = 8.0;  	% e0_c_A_fe6_i0 
	X_ITER(12) = 8.0;  	% e0_c_A_fe6_i3 
	X_ITER(13) = 8.0;  	% e0_c_A_fe7_i0 
	X_ITER(14) = 333.0;  	% e0_T_fe1_i1 
	X_ITER(15) = 333.0;  	% e0_T_fe1_i2 
	X_ITER(16) = 333.0;  	% e0_T_fe1_i3 
	X_ITER(17) = 333.0;  	% e0_T_fe2_i1 
	X_ITER(18) = 333.0;  	% e0_T_fe2_i2 
	X_ITER(19) = 333.0;  	% e0_T_fe2_i3 
	X_ITER(20) = 333.0;  	% e0_T_fe3_i1 
	X_ITER(21) = 333.0;  	% e0_T_fe3_i2 
	X_ITER(22) = 333.0;  	% e0_T_fe3_i3 
	X_ITER(23) = 333.0;  	% e0_T_fe4_i1 
	X_ITER(24) = 333.0;  	% e0_T_fe4_i2 
	X_ITER(25) = 333.0;  	% e0_T_fe4_i3 
	X_ITER(26) = 333.0;  	% e0_T_fe5_i1 
	X_ITER(27) = 333.0;  	% e0_T_fe5_i2 
	X_ITER(28) = 333.0;  	% e0_T_fe5_i3 
	X_ITER(29) = 333.0;  	% e0_T_fe6_i1 
	X_ITER(30) = 333.0;  	% e0_T_fe6_i2 
	X_ITER(31) = 333.0;  	% e0_T_fe6_i3 
	X_ITER(32) = 8.0;  	% e0_c_A_fe1_i1 
	X_ITER(33) = 8.0;  	% e0_c_A_fe1_i2 
	X_ITER(34) = 8.0;  	% e0_c_A_fe2_i1 
	X_ITER(35) = 8.0;  	% e0_c_A_fe2_i2 
	X_ITER(36) = 8.0;  	% e0_c_A_fe3_i1 
	X_ITER(37) = 8.0;  	% e0_c_A_fe3_i2 
	X_ITER(38) = 8.0;  	% e0_c_A_fe4_i1 
	X_ITER(39) = 8.0;  	% e0_c_A_fe4_i2 
	X_ITER(40) = 8.0;  	% e0_c_A_fe5_i1 
	X_ITER(41) = 8.0;  	% e0_c_A_fe5_i2 
	X_ITER(42) = 8.0;  	% e0_c_A_fe6_i1 
	X_ITER(43) = 8.0;  	% e0_c_A_fe6_i2 
	X_ITER(44) = 20.0;  	% e0_c_B_fe1_i1 
	X_ITER(45) = 20.0;  	% e0_c_B_fe1_i2 
	X_ITER(46) = 20.0;  	% e0_c_B_fe1_i3 
	X_ITER(47) = 20.0;  	% e0_c_B_fe2_i1 
	X_ITER(48) = 20.0;  	% e0_c_B_fe2_i2 
	X_ITER(49) = 20.0;  	% e0_c_B_fe2_i3 
	X_ITER(50) = 20.0;  	% e0_c_B_fe3_i1 
	X_ITER(51) = 20.0;  	% e0_c_B_fe3_i2 
	X_ITER(52) = 20.0;  	% e0_c_B_fe3_i3 
	X_ITER(53) = 20.0;  	% e0_c_B_fe4_i1 
	X_ITER(54) = 20.0;  	% e0_c_B_fe4_i2 
	X_ITER(55) = 20.0;  	% e0_c_B_fe4_i3 
	X_ITER(56) = 20.0;  	% e0_c_B_fe5_i1 
	X_ITER(57) = 20.0;  	% e0_c_B_fe5_i2 
	X_ITER(58) = 20.0;  	% e0_c_B_fe5_i3 
	X_ITER(59) = 20.0;  	% e0_c_B_fe6_i1 
	X_ITER(60) = 20.0;  	% e0_c_B_fe6_i2 
	X_ITER(61) = 20.0;  	% e0_c_B_fe6_i3 
	X_ITER(62) = 20.0;  	% e0_c_B_fe1_i0 
	X_ITER(63) = 20.0;  	% e0_c_B_fe2_i0 
	X_ITER(64) = 20.0;  	% e0_c_B_fe3_i0 
	X_ITER(65) = 20.0;  	% e0_c_B_fe4_i0 
	X_ITER(66) = 20.0;  	% e0_c_B_fe5_i0 
	X_ITER(67) = 20.0;  	% e0_c_B_fe6_i0 
	X_ITER(68) = 20.0;  	% e0_c_B_fe7_i0 
	X_ITER(69) = 0.0;  	% e0_c_C_fe1_i0 
	X_ITER(70) = 0.0;  	% e0_c_C_fe1_i3 
	X_ITER(71) = 0.0;  	% e0_c_C_fe2_i0 
	X_ITER(72) = 0.0;  	% e0_c_C_fe2_i3 
	X_ITER(73) = 0.0;  	% e0_c_C_fe3_i0 
	X_ITER(74) = 0.0;  	% e0_c_C_fe3_i3 
	X_ITER(75) = 0.0;  	% e0_c_C_fe4_i0 
	X_ITER(76) = 0.0;  	% e0_c_C_fe4_i3 
	X_ITER(77) = 0.0;  	% e0_c_C_fe5_i0 
	X_ITER(78) = 0.0;  	% e0_c_C_fe5_i3 
	X_ITER(79) = 0.0;  	% e0_c_C_fe6_i0 
	X_ITER(80) = 0.0;  	% e0_c_C_fe6_i3 
	X_ITER(81) = 0.0;  	% e0_c_C_fe7_i0 
	X_ITER(82) = 0.0;  	% e0_c_C_fe1_i1 
	X_ITER(83) = 0.0;  	% e0_c_C_fe1_i2 
	X_ITER(84) = 0.0;  	% e0_c_C_fe2_i1 
	X_ITER(85) = 0.0;  	% e0_c_C_fe2_i2 
	X_ITER(86) = 0.0;  	% e0_c_C_fe3_i1 
	X_ITER(87) = 0.0;  	% e0_c_C_fe3_i2 
	X_ITER(88) = 0.0;  	% e0_c_C_fe4_i1 
	X_ITER(89) = 0.0;  	% e0_c_C_fe4_i2 
	X_ITER(90) = 0.0;  	% e0_c_C_fe5_i1 
	X_ITER(91) = 0.0;  	% e0_c_C_fe5_i2 
	X_ITER(92) = 0.0;  	% e0_c_C_fe6_i1 
	X_ITER(93) = 0.0;  	% e0_c_C_fe6_i2 
	X_ITER(94) = 333.0;  	% e0_T_fe1_i0 
	X_ITER(95) = 333.0;  	% e0_T_fe2_i0 
	X_ITER(96) = 333.0;  	% e0_T_fe3_i0 
	X_ITER(97) = 333.0;  	% e0_T_fe4_i0 
	X_ITER(98) = 333.0;  	% e0_T_fe5_i0 
	X_ITER(99) = 333.0;  	% e0_T_fe6_i0 
	X_ITER(100) = 333.0;  	% e0_T_fe7_i0 

	% load parameters 
	PARAMS(1) = 1.167839841902244;  	% e0_dldgreek_tau_i2_j1 
	PARAMS(2) = 0.775254648382856;  	% e0_dldgreek_tau_i2_j2 
	PARAMS(3) = -7.531972331053937;  	% e0_dldgreek_tau_i2_j3 
	PARAMS(4) = -0.253197259961796;  	% e0_dldgreek_tau_i3_j1 
	PARAMS(5) = 1.053197461577804;  	% e0_dldgreek_tau_i3_j2 
	PARAMS(6) = 600.0;  	% e0_greek_Deltat 
	PARAMS(7) = 1.0;  	% e0_V 
	PARAMS(8) = 69000.0;  	% e0_E_A 
	PARAMS(9) = 72000.0;  	% e0_E_B 
	PARAMS(10) = 6.5E-4;  	% e0_F 
	PARAMS(11) = 8.314;  	% e0_R 
	PARAMS(12) = 800.0;  	% e0_greek_rho 
	PARAMS(13) = 45000.0;  	% e0_H_A 
	PARAMS(14) = -55000.0;  	% e0_H_B 
	PARAMS(15) = 333.0;  	% e0_T_Cool 
	PARAMS(16) = 333.0;  	% e0_T_In 
	PARAMS(17) = 1.4;  	% e0_UA 
	PARAMS(18) = 3.5;  	% e0_cp 
	PARAMS(19) = 15.0;  	% e0_c_B0 
	PARAMS(20) = 1.739387967160278;  	% e0_dldgreek_tau_i0_j2 
	PARAMS(21) = -3.000000252020032;  	% e0_dldgreek_tau_i0_j3 
	PARAMS(22) = 0.0;  	% e0_x_Init_x3 
	PARAMS(23) = 3.224746191683931;  	% e0_dldgreek_tau_i1_j1 
	PARAMS(24) = 5.0;  	% e0_c_A0 
	PARAMS(25) = -3.567840077120938;  	% e0_dldgreek_tau_i1_j2 
	PARAMS(26) = 5.53197241506063;  	% e0_dldgreek_tau_i1_j3 
	PARAMS(27) = -4.139388773624379;  	% e0_dldgreek_tau_i0_j1 
	PARAMS(28) = 0.0;  	% e0_c_C0 
	PARAMS(29) = 8.0;  	% e0_x_Init_x1 
	PARAMS(30) = 333.0;  	% e0_x_Init_x4 
	PARAMS(31) = 5.00000016801334;  	% e0_dldgreek_tau_i3_j3 
	PARAMS(32) = 5000000.0;  	% e0_k_A 
	PARAMS(33) = 1.0E7;  	% e0_k_B 
	PARAMS(34) = 20.0;  	% e0_x_Init_x2 

	options = optimset('MaxIter',1000,'MaxFunEvals',1e4,'TolFun',1e-6,'Display','Iter');
	RES = fsolve (@( x_iter )getFunVal(x_iter,PARAMS),X_ITER,options);

	ROOTS = getFunVal(RES,PARAMS);
	ROOTS = ROOTS';

	displayResults(RES);

end


function[Y] = getFunVal(X_ITER,PARAMS)

%
% Calculate the function value of a normalized equation system. 
%
	% read out variables  
	e0_c_A_fe1_i0 = X_ITER(1); 
	e0_c_A_fe1_i3 = X_ITER(2); 
	e0_c_A_fe2_i0 = X_ITER(3); 
	e0_c_A_fe2_i3 = X_ITER(4); 
	e0_c_A_fe3_i0 = X_ITER(5); 
	e0_c_A_fe3_i3 = X_ITER(6); 
	e0_c_A_fe4_i0 = X_ITER(7); 
	e0_c_A_fe4_i3 = X_ITER(8); 
	e0_c_A_fe5_i0 = X_ITER(9); 
	e0_c_A_fe5_i3 = X_ITER(10); 
	e0_c_A_fe6_i0 = X_ITER(11); 
	e0_c_A_fe6_i3 = X_ITER(12); 
	e0_c_A_fe7_i0 = X_ITER(13); 
	e0_T_fe1_i1 = X_ITER(14); 
	e0_T_fe1_i2 = X_ITER(15); 
	e0_T_fe1_i3 = X_ITER(16); 
	e0_T_fe2_i1 = X_ITER(17); 
	e0_T_fe2_i2 = X_ITER(18); 
	e0_T_fe2_i3 = X_ITER(19); 
	e0_T_fe3_i1 = X_ITER(20); 
	e0_T_fe3_i2 = X_ITER(21); 
	e0_T_fe3_i3 = X_ITER(22); 
	e0_T_fe4_i1 = X_ITER(23); 
	e0_T_fe4_i2 = X_ITER(24); 
	e0_T_fe4_i3 = X_ITER(25); 
	e0_T_fe5_i1 = X_ITER(26); 
	e0_T_fe5_i2 = X_ITER(27); 
	e0_T_fe5_i3 = X_ITER(28); 
	e0_T_fe6_i1 = X_ITER(29); 
	e0_T_fe6_i2 = X_ITER(30); 
	e0_T_fe6_i3 = X_ITER(31); 
	e0_c_A_fe1_i1 = X_ITER(32); 
	e0_c_A_fe1_i2 = X_ITER(33); 
	e0_c_A_fe2_i1 = X_ITER(34); 
	e0_c_A_fe2_i2 = X_ITER(35); 
	e0_c_A_fe3_i1 = X_ITER(36); 
	e0_c_A_fe3_i2 = X_ITER(37); 
	e0_c_A_fe4_i1 = X_ITER(38); 
	e0_c_A_fe4_i2 = X_ITER(39); 
	e0_c_A_fe5_i1 = X_ITER(40); 
	e0_c_A_fe5_i2 = X_ITER(41); 
	e0_c_A_fe6_i1 = X_ITER(42); 
	e0_c_A_fe6_i2 = X_ITER(43); 
	e0_c_B_fe1_i1 = X_ITER(44); 
	e0_c_B_fe1_i2 = X_ITER(45); 
	e0_c_B_fe1_i3 = X_ITER(46); 
	e0_c_B_fe2_i1 = X_ITER(47); 
	e0_c_B_fe2_i2 = X_ITER(48); 
	e0_c_B_fe2_i3 = X_ITER(49); 
	e0_c_B_fe3_i1 = X_ITER(50); 
	e0_c_B_fe3_i2 = X_ITER(51); 
	e0_c_B_fe3_i3 = X_ITER(52); 
	e0_c_B_fe4_i1 = X_ITER(53); 
	e0_c_B_fe4_i2 = X_ITER(54); 
	e0_c_B_fe4_i3 = X_ITER(55); 
	e0_c_B_fe5_i1 = X_ITER(56); 
	e0_c_B_fe5_i2 = X_ITER(57); 
	e0_c_B_fe5_i3 = X_ITER(58); 
	e0_c_B_fe6_i1 = X_ITER(59); 
	e0_c_B_fe6_i2 = X_ITER(60); 
	e0_c_B_fe6_i3 = X_ITER(61); 
	e0_c_B_fe1_i0 = X_ITER(62); 
	e0_c_B_fe2_i0 = X_ITER(63); 
	e0_c_B_fe3_i0 = X_ITER(64); 
	e0_c_B_fe4_i0 = X_ITER(65); 
	e0_c_B_fe5_i0 = X_ITER(66); 
	e0_c_B_fe6_i0 = X_ITER(67); 
	e0_c_B_fe7_i0 = X_ITER(68); 
	e0_c_C_fe1_i0 = X_ITER(69); 
	e0_c_C_fe1_i3 = X_ITER(70); 
	e0_c_C_fe2_i0 = X_ITER(71); 
	e0_c_C_fe2_i3 = X_ITER(72); 
	e0_c_C_fe3_i0 = X_ITER(73); 
	e0_c_C_fe3_i3 = X_ITER(74); 
	e0_c_C_fe4_i0 = X_ITER(75); 
	e0_c_C_fe4_i3 = X_ITER(76); 
	e0_c_C_fe5_i0 = X_ITER(77); 
	e0_c_C_fe5_i3 = X_ITER(78); 
	e0_c_C_fe6_i0 = X_ITER(79); 
	e0_c_C_fe6_i3 = X_ITER(80); 
	e0_c_C_fe7_i0 = X_ITER(81); 
	e0_c_C_fe1_i1 = X_ITER(82); 
	e0_c_C_fe1_i2 = X_ITER(83); 
	e0_c_C_fe2_i1 = X_ITER(84); 
	e0_c_C_fe2_i2 = X_ITER(85); 
	e0_c_C_fe3_i1 = X_ITER(86); 
	e0_c_C_fe3_i2 = X_ITER(87); 
	e0_c_C_fe4_i1 = X_ITER(88); 
	e0_c_C_fe4_i2 = X_ITER(89); 
	e0_c_C_fe5_i1 = X_ITER(90); 
	e0_c_C_fe5_i2 = X_ITER(91); 
	e0_c_C_fe6_i1 = X_ITER(92); 
	e0_c_C_fe6_i2 = X_ITER(93); 
	e0_T_fe1_i0 = X_ITER(94); 
	e0_T_fe2_i0 = X_ITER(95); 
	e0_T_fe3_i0 = X_ITER(96); 
	e0_T_fe4_i0 = X_ITER(97); 
	e0_T_fe5_i0 = X_ITER(98); 
	e0_T_fe6_i0 = X_ITER(99); 
	e0_T_fe7_i0 = X_ITER(100); 

	% read out parameters 
	e0_dldgreek_tau_i2_j1 = PARAMS(1); 
	e0_dldgreek_tau_i2_j2 = PARAMS(2); 
	e0_dldgreek_tau_i2_j3 = PARAMS(3); 
	e0_dldgreek_tau_i3_j1 = PARAMS(4); 
	e0_dldgreek_tau_i3_j2 = PARAMS(5); 
	e0_greek_Deltat = PARAMS(6); 
	e0_V = PARAMS(7); 
	e0_E_A = PARAMS(8); 
	e0_E_B = PARAMS(9); 
	e0_F = PARAMS(10); 
	e0_R = PARAMS(11); 
	e0_greek_rho = PARAMS(12); 
	e0_H_A = PARAMS(13); 
	e0_H_B = PARAMS(14); 
	e0_T_Cool = PARAMS(15); 
	e0_T_In = PARAMS(16); 
	e0_UA = PARAMS(17); 
	e0_cp = PARAMS(18); 
	e0_c_B0 = PARAMS(19); 
	e0_dldgreek_tau_i0_j2 = PARAMS(20); 
	e0_dldgreek_tau_i0_j3 = PARAMS(21); 
	e0_x_Init_x3 = PARAMS(22); 
	e0_dldgreek_tau_i1_j1 = PARAMS(23); 
	e0_c_A0 = PARAMS(24); 
	e0_dldgreek_tau_i1_j2 = PARAMS(25); 
	e0_dldgreek_tau_i1_j3 = PARAMS(26); 
	e0_dldgreek_tau_i0_j1 = PARAMS(27); 
	e0_c_C0 = PARAMS(28); 
	e0_x_Init_x1 = PARAMS(29); 
	e0_x_Init_x4 = PARAMS(30); 
	e0_dldgreek_tau_i3_j3 = PARAMS(31); 
	e0_k_A = PARAMS(32); 
	e0_k_B = PARAMS(33); 
	e0_x_Init_x2 = PARAMS(34); 

	% perform direct function calls 

	% evaluate the function values  
	Y(1) = e0_c_A_fe1_i0 - ( e0_x_Init_x1 ); 
	Y(2) = e0_c_A_fe1_i3 - ( e0_c_A_fe2_i0 ); 
	Y(3) = e0_c_A_fe2_i3 - ( e0_c_A_fe3_i0 ); 
	Y(4) = e0_c_A_fe3_i3 - ( e0_c_A_fe4_i0 ); 
	Y(5) = e0_c_A_fe4_i3 - ( e0_c_A_fe5_i0 ); 
	Y(6) = e0_c_A_fe5_i3 - ( e0_c_A_fe6_i0 ); 
	Y(7) = e0_c_A_fe6_i3 - ( e0_c_A_fe7_i0 ); 
	Y(8) = ( e0_c_A_fe1_i0 * e0_dldgreek_tau_i0_j1 + e0_c_A_fe1_i1 * e0_dldgreek_tau_i1_j1 + e0_c_A_fe1_i2 * e0_dldgreek_tau_i2_j1 + e0_c_A_fe1_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe1_i1 ) - e0_k_A * e0_c_A_fe1_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i1 ) ) + e0_k_B * e0_c_B_fe1_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i1 ) ) ); 
	Y(9) = ( e0_c_A_fe1_i0 * e0_dldgreek_tau_i0_j2 + e0_c_A_fe1_i1 * e0_dldgreek_tau_i1_j2 + e0_c_A_fe1_i2 * e0_dldgreek_tau_i2_j2 + e0_c_A_fe1_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe1_i2 ) - e0_k_A * e0_c_A_fe1_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i2 ) ) + e0_k_B * e0_c_B_fe1_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i2 ) ) ); 
	Y(10) = ( e0_c_A_fe1_i0 * e0_dldgreek_tau_i0_j3 + e0_c_A_fe1_i1 * e0_dldgreek_tau_i1_j3 + e0_c_A_fe1_i2 * e0_dldgreek_tau_i2_j3 + e0_c_A_fe1_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe1_i3 ) - e0_k_A * e0_c_A_fe1_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i3 ) ) + e0_k_B * e0_c_B_fe1_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i3 ) ) ); 
	Y(11) = ( e0_c_A_fe2_i0 * e0_dldgreek_tau_i0_j1 + e0_c_A_fe2_i1 * e0_dldgreek_tau_i1_j1 + e0_c_A_fe2_i2 * e0_dldgreek_tau_i2_j1 + e0_c_A_fe2_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe2_i1 ) - e0_k_A * e0_c_A_fe2_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i1 ) ) + e0_k_B * e0_c_B_fe2_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i1 ) ) ); 
	Y(12) = ( e0_c_A_fe2_i0 * e0_dldgreek_tau_i0_j2 + e0_c_A_fe2_i1 * e0_dldgreek_tau_i1_j2 + e0_c_A_fe2_i2 * e0_dldgreek_tau_i2_j2 + e0_c_A_fe2_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe2_i2 ) - e0_k_A * e0_c_A_fe2_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i2 ) ) + e0_k_B * e0_c_B_fe2_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i2 ) ) ); 
	Y(13) = ( e0_c_A_fe2_i0 * e0_dldgreek_tau_i0_j3 + e0_c_A_fe2_i1 * e0_dldgreek_tau_i1_j3 + e0_c_A_fe2_i2 * e0_dldgreek_tau_i2_j3 + e0_c_A_fe2_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe2_i3 ) - e0_k_A * e0_c_A_fe2_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i3 ) ) + e0_k_B * e0_c_B_fe2_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i3 ) ) ); 
	Y(14) = ( e0_c_A_fe3_i0 * e0_dldgreek_tau_i0_j1 + e0_c_A_fe3_i1 * e0_dldgreek_tau_i1_j1 + e0_c_A_fe3_i2 * e0_dldgreek_tau_i2_j1 + e0_c_A_fe3_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe3_i1 ) - e0_k_A * e0_c_A_fe3_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i1 ) ) + e0_k_B * e0_c_B_fe3_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i1 ) ) ); 
	Y(15) = ( e0_c_A_fe3_i0 * e0_dldgreek_tau_i0_j2 + e0_c_A_fe3_i1 * e0_dldgreek_tau_i1_j2 + e0_c_A_fe3_i2 * e0_dldgreek_tau_i2_j2 + e0_c_A_fe3_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe3_i2 ) - e0_k_A * e0_c_A_fe3_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i2 ) ) + e0_k_B * e0_c_B_fe3_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i2 ) ) ); 
	Y(16) = ( e0_c_A_fe3_i0 * e0_dldgreek_tau_i0_j3 + e0_c_A_fe3_i1 * e0_dldgreek_tau_i1_j3 + e0_c_A_fe3_i2 * e0_dldgreek_tau_i2_j3 + e0_c_A_fe3_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe3_i3 ) - e0_k_A * e0_c_A_fe3_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i3 ) ) + e0_k_B * e0_c_B_fe3_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i3 ) ) ); 
	Y(17) = ( e0_c_A_fe4_i0 * e0_dldgreek_tau_i0_j1 + e0_c_A_fe4_i1 * e0_dldgreek_tau_i1_j1 + e0_c_A_fe4_i2 * e0_dldgreek_tau_i2_j1 + e0_c_A_fe4_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe4_i1 ) - e0_k_A * e0_c_A_fe4_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i1 ) ) + e0_k_B * e0_c_B_fe4_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i1 ) ) ); 
	Y(18) = ( e0_c_A_fe4_i0 * e0_dldgreek_tau_i0_j2 + e0_c_A_fe4_i1 * e0_dldgreek_tau_i1_j2 + e0_c_A_fe4_i2 * e0_dldgreek_tau_i2_j2 + e0_c_A_fe4_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe4_i2 ) - e0_k_A * e0_c_A_fe4_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i2 ) ) + e0_k_B * e0_c_B_fe4_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i2 ) ) ); 
	Y(19) = ( e0_c_A_fe4_i0 * e0_dldgreek_tau_i0_j3 + e0_c_A_fe4_i1 * e0_dldgreek_tau_i1_j3 + e0_c_A_fe4_i2 * e0_dldgreek_tau_i2_j3 + e0_c_A_fe4_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe4_i3 ) - e0_k_A * e0_c_A_fe4_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i3 ) ) + e0_k_B * e0_c_B_fe4_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i3 ) ) ); 
	Y(20) = ( e0_c_A_fe5_i0 * e0_dldgreek_tau_i0_j1 + e0_c_A_fe5_i1 * e0_dldgreek_tau_i1_j1 + e0_c_A_fe5_i2 * e0_dldgreek_tau_i2_j1 + e0_c_A_fe5_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe5_i1 ) - e0_k_A * e0_c_A_fe5_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i1 ) ) + e0_k_B * e0_c_B_fe5_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i1 ) ) ); 
	Y(21) = ( e0_c_A_fe5_i0 * e0_dldgreek_tau_i0_j2 + e0_c_A_fe5_i1 * e0_dldgreek_tau_i1_j2 + e0_c_A_fe5_i2 * e0_dldgreek_tau_i2_j2 + e0_c_A_fe5_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe5_i2 ) - e0_k_A * e0_c_A_fe5_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i2 ) ) + e0_k_B * e0_c_B_fe5_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i2 ) ) ); 
	Y(22) = ( e0_c_A_fe5_i0 * e0_dldgreek_tau_i0_j3 + e0_c_A_fe5_i1 * e0_dldgreek_tau_i1_j3 + e0_c_A_fe5_i2 * e0_dldgreek_tau_i2_j3 + e0_c_A_fe5_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe5_i3 ) - e0_k_A * e0_c_A_fe5_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i3 ) ) + e0_k_B * e0_c_B_fe5_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i3 ) ) ); 
	Y(23) = ( e0_c_A_fe6_i0 * e0_dldgreek_tau_i0_j1 + e0_c_A_fe6_i1 * e0_dldgreek_tau_i1_j1 + e0_c_A_fe6_i2 * e0_dldgreek_tau_i2_j1 + e0_c_A_fe6_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe6_i1 ) - e0_k_A * e0_c_A_fe6_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i1 ) ) + e0_k_B * e0_c_B_fe6_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i1 ) ) ); 
	Y(24) = ( e0_c_A_fe6_i0 * e0_dldgreek_tau_i0_j2 + e0_c_A_fe6_i1 * e0_dldgreek_tau_i1_j2 + e0_c_A_fe6_i2 * e0_dldgreek_tau_i2_j2 + e0_c_A_fe6_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe6_i2 ) - e0_k_A * e0_c_A_fe6_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i2 ) ) + e0_k_B * e0_c_B_fe6_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i2 ) ) ); 
	Y(25) = ( e0_c_A_fe6_i0 * e0_dldgreek_tau_i0_j3 + e0_c_A_fe6_i1 * e0_dldgreek_tau_i1_j3 + e0_c_A_fe6_i2 * e0_dldgreek_tau_i2_j3 + e0_c_A_fe6_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_A0 - e0_c_A_fe6_i3 ) - e0_k_A * e0_c_A_fe6_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i3 ) ) + e0_k_B * e0_c_B_fe6_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i3 ) ) ); 
	Y(26) = e0_c_B_fe1_i0 - ( e0_x_Init_x2 ); 
	Y(27) = e0_c_B_fe1_i3 - ( e0_c_B_fe2_i0 ); 
	Y(28) = e0_c_B_fe2_i3 - ( e0_c_B_fe3_i0 ); 
	Y(29) = e0_c_B_fe3_i3 - ( e0_c_B_fe4_i0 ); 
	Y(30) = e0_c_B_fe4_i3 - ( e0_c_B_fe5_i0 ); 
	Y(31) = e0_c_B_fe5_i3 - ( e0_c_B_fe6_i0 ); 
	Y(32) = e0_c_B_fe6_i3 - ( e0_c_B_fe7_i0 ); 
	Y(33) = ( e0_c_B_fe1_i0 * e0_dldgreek_tau_i0_j1 + e0_c_B_fe1_i1 * e0_dldgreek_tau_i1_j1 + e0_c_B_fe1_i2 * e0_dldgreek_tau_i2_j1 + e0_c_B_fe1_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe1_i1 ) - e0_k_B * e0_c_B_fe1_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i1 ) ) ); 
	Y(34) = ( e0_c_B_fe1_i0 * e0_dldgreek_tau_i0_j2 + e0_c_B_fe1_i1 * e0_dldgreek_tau_i1_j2 + e0_c_B_fe1_i2 * e0_dldgreek_tau_i2_j2 + e0_c_B_fe1_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe1_i2 ) - e0_k_B * e0_c_B_fe1_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i2 ) ) ); 
	Y(35) = ( e0_c_B_fe1_i0 * e0_dldgreek_tau_i0_j3 + e0_c_B_fe1_i1 * e0_dldgreek_tau_i1_j3 + e0_c_B_fe1_i2 * e0_dldgreek_tau_i2_j3 + e0_c_B_fe1_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe1_i3 ) - e0_k_B * e0_c_B_fe1_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i3 ) ) ); 
	Y(36) = ( e0_c_B_fe2_i0 * e0_dldgreek_tau_i0_j1 + e0_c_B_fe2_i1 * e0_dldgreek_tau_i1_j1 + e0_c_B_fe2_i2 * e0_dldgreek_tau_i2_j1 + e0_c_B_fe2_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe2_i1 ) - e0_k_B * e0_c_B_fe2_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i1 ) ) ); 
	Y(37) = ( e0_c_B_fe2_i0 * e0_dldgreek_tau_i0_j2 + e0_c_B_fe2_i1 * e0_dldgreek_tau_i1_j2 + e0_c_B_fe2_i2 * e0_dldgreek_tau_i2_j2 + e0_c_B_fe2_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe2_i2 ) - e0_k_B * e0_c_B_fe2_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i2 ) ) ); 
	Y(38) = ( e0_c_B_fe2_i0 * e0_dldgreek_tau_i0_j3 + e0_c_B_fe2_i1 * e0_dldgreek_tau_i1_j3 + e0_c_B_fe2_i2 * e0_dldgreek_tau_i2_j3 + e0_c_B_fe2_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe2_i3 ) - e0_k_B * e0_c_B_fe2_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i3 ) ) ); 
	Y(39) = ( e0_c_B_fe3_i0 * e0_dldgreek_tau_i0_j1 + e0_c_B_fe3_i1 * e0_dldgreek_tau_i1_j1 + e0_c_B_fe3_i2 * e0_dldgreek_tau_i2_j1 + e0_c_B_fe3_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe3_i1 ) - e0_k_B * e0_c_B_fe3_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i1 ) ) ); 
	Y(40) = ( e0_c_B_fe3_i0 * e0_dldgreek_tau_i0_j2 + e0_c_B_fe3_i1 * e0_dldgreek_tau_i1_j2 + e0_c_B_fe3_i2 * e0_dldgreek_tau_i2_j2 + e0_c_B_fe3_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe3_i2 ) - e0_k_B * e0_c_B_fe3_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i2 ) ) ); 
	Y(41) = ( e0_c_B_fe3_i0 * e0_dldgreek_tau_i0_j3 + e0_c_B_fe3_i1 * e0_dldgreek_tau_i1_j3 + e0_c_B_fe3_i2 * e0_dldgreek_tau_i2_j3 + e0_c_B_fe3_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe3_i3 ) - e0_k_B * e0_c_B_fe3_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i3 ) ) ); 
	Y(42) = ( e0_c_B_fe4_i0 * e0_dldgreek_tau_i0_j1 + e0_c_B_fe4_i1 * e0_dldgreek_tau_i1_j1 + e0_c_B_fe4_i2 * e0_dldgreek_tau_i2_j1 + e0_c_B_fe4_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe4_i1 ) - e0_k_B * e0_c_B_fe4_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i1 ) ) ); 
	Y(43) = ( e0_c_B_fe4_i0 * e0_dldgreek_tau_i0_j2 + e0_c_B_fe4_i1 * e0_dldgreek_tau_i1_j2 + e0_c_B_fe4_i2 * e0_dldgreek_tau_i2_j2 + e0_c_B_fe4_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe4_i2 ) - e0_k_B * e0_c_B_fe4_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i2 ) ) ); 
	Y(44) = ( e0_c_B_fe4_i0 * e0_dldgreek_tau_i0_j3 + e0_c_B_fe4_i1 * e0_dldgreek_tau_i1_j3 + e0_c_B_fe4_i2 * e0_dldgreek_tau_i2_j3 + e0_c_B_fe4_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe4_i3 ) - e0_k_B * e0_c_B_fe4_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i3 ) ) ); 
	Y(45) = ( e0_c_B_fe5_i0 * e0_dldgreek_tau_i0_j1 + e0_c_B_fe5_i1 * e0_dldgreek_tau_i1_j1 + e0_c_B_fe5_i2 * e0_dldgreek_tau_i2_j1 + e0_c_B_fe5_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe5_i1 ) - e0_k_B * e0_c_B_fe5_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i1 ) ) ); 
	Y(46) = ( e0_c_B_fe5_i0 * e0_dldgreek_tau_i0_j2 + e0_c_B_fe5_i1 * e0_dldgreek_tau_i1_j2 + e0_c_B_fe5_i2 * e0_dldgreek_tau_i2_j2 + e0_c_B_fe5_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe5_i2 ) - e0_k_B * e0_c_B_fe5_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i2 ) ) ); 
	Y(47) = ( e0_c_B_fe5_i0 * e0_dldgreek_tau_i0_j3 + e0_c_B_fe5_i1 * e0_dldgreek_tau_i1_j3 + e0_c_B_fe5_i2 * e0_dldgreek_tau_i2_j3 + e0_c_B_fe5_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe5_i3 ) - e0_k_B * e0_c_B_fe5_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i3 ) ) ); 
	Y(48) = ( e0_c_B_fe6_i0 * e0_dldgreek_tau_i0_j1 + e0_c_B_fe6_i1 * e0_dldgreek_tau_i1_j1 + e0_c_B_fe6_i2 * e0_dldgreek_tau_i2_j1 + e0_c_B_fe6_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe6_i1 ) - e0_k_B * e0_c_B_fe6_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i1 ) ) ); 
	Y(49) = ( e0_c_B_fe6_i0 * e0_dldgreek_tau_i0_j2 + e0_c_B_fe6_i1 * e0_dldgreek_tau_i1_j2 + e0_c_B_fe6_i2 * e0_dldgreek_tau_i2_j2 + e0_c_B_fe6_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe6_i2 ) - e0_k_B * e0_c_B_fe6_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i2 ) ) ); 
	Y(50) = ( e0_c_B_fe6_i0 * e0_dldgreek_tau_i0_j3 + e0_c_B_fe6_i1 * e0_dldgreek_tau_i1_j3 + e0_c_B_fe6_i2 * e0_dldgreek_tau_i2_j3 + e0_c_B_fe6_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_B0 - e0_c_B_fe6_i3 ) - e0_k_B * e0_c_B_fe6_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i3 ) ) ); 
	Y(51) = e0_c_C_fe1_i0 - ( e0_x_Init_x3 ); 
	Y(52) = e0_c_C_fe1_i3 - ( e0_c_C_fe2_i0 ); 
	Y(53) = e0_c_C_fe2_i3 - ( e0_c_C_fe3_i0 ); 
	Y(54) = e0_c_C_fe3_i3 - ( e0_c_C_fe4_i0 ); 
	Y(55) = e0_c_C_fe4_i3 - ( e0_c_C_fe5_i0 ); 
	Y(56) = e0_c_C_fe5_i3 - ( e0_c_C_fe6_i0 ); 
	Y(57) = e0_c_C_fe6_i3 - ( e0_c_C_fe7_i0 ); 
	Y(58) = ( e0_c_C_fe1_i0 * e0_dldgreek_tau_i0_j1 + e0_c_C_fe1_i1 * e0_dldgreek_tau_i1_j1 + e0_c_C_fe1_i2 * e0_dldgreek_tau_i2_j1 + e0_c_C_fe1_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe1_i1 ) + e0_k_A * e0_c_A_fe1_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i1 ) ) ); 
	Y(59) = ( e0_c_C_fe1_i0 * e0_dldgreek_tau_i0_j2 + e0_c_C_fe1_i1 * e0_dldgreek_tau_i1_j2 + e0_c_C_fe1_i2 * e0_dldgreek_tau_i2_j2 + e0_c_C_fe1_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe1_i2 ) + e0_k_A * e0_c_A_fe1_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i2 ) ) ); 
	Y(60) = ( e0_c_C_fe1_i0 * e0_dldgreek_tau_i0_j3 + e0_c_C_fe1_i1 * e0_dldgreek_tau_i1_j3 + e0_c_C_fe1_i2 * e0_dldgreek_tau_i2_j3 + e0_c_C_fe1_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe1_i3 ) + e0_k_A * e0_c_A_fe1_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i3 ) ) ); 
	Y(61) = ( e0_c_C_fe2_i0 * e0_dldgreek_tau_i0_j1 + e0_c_C_fe2_i1 * e0_dldgreek_tau_i1_j1 + e0_c_C_fe2_i2 * e0_dldgreek_tau_i2_j1 + e0_c_C_fe2_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe2_i1 ) + e0_k_A * e0_c_A_fe2_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i1 ) ) ); 
	Y(62) = ( e0_c_C_fe2_i0 * e0_dldgreek_tau_i0_j2 + e0_c_C_fe2_i1 * e0_dldgreek_tau_i1_j2 + e0_c_C_fe2_i2 * e0_dldgreek_tau_i2_j2 + e0_c_C_fe2_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe2_i2 ) + e0_k_A * e0_c_A_fe2_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i2 ) ) ); 
	Y(63) = ( e0_c_C_fe2_i0 * e0_dldgreek_tau_i0_j3 + e0_c_C_fe2_i1 * e0_dldgreek_tau_i1_j3 + e0_c_C_fe2_i2 * e0_dldgreek_tau_i2_j3 + e0_c_C_fe2_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe2_i3 ) + e0_k_A * e0_c_A_fe2_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i3 ) ) ); 
	Y(64) = ( e0_c_C_fe3_i0 * e0_dldgreek_tau_i0_j1 + e0_c_C_fe3_i1 * e0_dldgreek_tau_i1_j1 + e0_c_C_fe3_i2 * e0_dldgreek_tau_i2_j1 + e0_c_C_fe3_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe3_i1 ) + e0_k_A * e0_c_A_fe3_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i1 ) ) ); 
	Y(65) = ( e0_c_C_fe3_i0 * e0_dldgreek_tau_i0_j2 + e0_c_C_fe3_i1 * e0_dldgreek_tau_i1_j2 + e0_c_C_fe3_i2 * e0_dldgreek_tau_i2_j2 + e0_c_C_fe3_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe3_i2 ) + e0_k_A * e0_c_A_fe3_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i2 ) ) ); 
	Y(66) = ( e0_c_C_fe3_i0 * e0_dldgreek_tau_i0_j3 + e0_c_C_fe3_i1 * e0_dldgreek_tau_i1_j3 + e0_c_C_fe3_i2 * e0_dldgreek_tau_i2_j3 + e0_c_C_fe3_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe3_i3 ) + e0_k_A * e0_c_A_fe3_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i3 ) ) ); 
	Y(67) = ( e0_c_C_fe4_i0 * e0_dldgreek_tau_i0_j1 + e0_c_C_fe4_i1 * e0_dldgreek_tau_i1_j1 + e0_c_C_fe4_i2 * e0_dldgreek_tau_i2_j1 + e0_c_C_fe4_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe4_i1 ) + e0_k_A * e0_c_A_fe4_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i1 ) ) ); 
	Y(68) = ( e0_c_C_fe4_i0 * e0_dldgreek_tau_i0_j2 + e0_c_C_fe4_i1 * e0_dldgreek_tau_i1_j2 + e0_c_C_fe4_i2 * e0_dldgreek_tau_i2_j2 + e0_c_C_fe4_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe4_i2 ) + e0_k_A * e0_c_A_fe4_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i2 ) ) ); 
	Y(69) = ( e0_c_C_fe4_i0 * e0_dldgreek_tau_i0_j3 + e0_c_C_fe4_i1 * e0_dldgreek_tau_i1_j3 + e0_c_C_fe4_i2 * e0_dldgreek_tau_i2_j3 + e0_c_C_fe4_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe4_i3 ) + e0_k_A * e0_c_A_fe4_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i3 ) ) ); 
	Y(70) = ( e0_c_C_fe5_i0 * e0_dldgreek_tau_i0_j1 + e0_c_C_fe5_i1 * e0_dldgreek_tau_i1_j1 + e0_c_C_fe5_i2 * e0_dldgreek_tau_i2_j1 + e0_c_C_fe5_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe5_i1 ) + e0_k_A * e0_c_A_fe5_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i1 ) ) ); 
	Y(71) = ( e0_c_C_fe5_i0 * e0_dldgreek_tau_i0_j2 + e0_c_C_fe5_i1 * e0_dldgreek_tau_i1_j2 + e0_c_C_fe5_i2 * e0_dldgreek_tau_i2_j2 + e0_c_C_fe5_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe5_i2 ) + e0_k_A * e0_c_A_fe5_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i2 ) ) ); 
	Y(72) = ( e0_c_C_fe5_i0 * e0_dldgreek_tau_i0_j3 + e0_c_C_fe5_i1 * e0_dldgreek_tau_i1_j3 + e0_c_C_fe5_i2 * e0_dldgreek_tau_i2_j3 + e0_c_C_fe5_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe5_i3 ) + e0_k_A * e0_c_A_fe5_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i3 ) ) ); 
	Y(73) = ( e0_c_C_fe6_i0 * e0_dldgreek_tau_i0_j1 + e0_c_C_fe6_i1 * e0_dldgreek_tau_i1_j1 + e0_c_C_fe6_i2 * e0_dldgreek_tau_i2_j1 + e0_c_C_fe6_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe6_i1 ) + e0_k_A * e0_c_A_fe6_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i1 ) ) ); 
	Y(74) = ( e0_c_C_fe6_i0 * e0_dldgreek_tau_i0_j2 + e0_c_C_fe6_i1 * e0_dldgreek_tau_i1_j2 + e0_c_C_fe6_i2 * e0_dldgreek_tau_i2_j2 + e0_c_C_fe6_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe6_i2 ) + e0_k_A * e0_c_A_fe6_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i2 ) ) ); 
	Y(75) = ( e0_c_C_fe6_i0 * e0_dldgreek_tau_i0_j3 + e0_c_C_fe6_i1 * e0_dldgreek_tau_i1_j3 + e0_c_C_fe6_i2 * e0_dldgreek_tau_i2_j3 + e0_c_C_fe6_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_c_C0 - e0_c_C_fe6_i3 ) + e0_k_A * e0_c_A_fe6_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i3 ) ) ); 
	Y(76) = e0_T_fe1_i0 - ( e0_x_Init_x4 ); 
	Y(77) = e0_T_fe1_i3 - ( e0_T_fe2_i0 ); 
	Y(78) = e0_T_fe2_i3 - ( e0_T_fe3_i0 ); 
	Y(79) = e0_T_fe3_i3 - ( e0_T_fe4_i0 ); 
	Y(80) = e0_T_fe4_i3 - ( e0_T_fe5_i0 ); 
	Y(81) = e0_T_fe5_i3 - ( e0_T_fe6_i0 ); 
	Y(82) = e0_T_fe6_i3 - ( e0_T_fe7_i0 ); 
	Y(83) = ( e0_T_fe1_i0 * e0_dldgreek_tau_i0_j1 + e0_T_fe1_i1 * e0_dldgreek_tau_i1_j1 + e0_T_fe1_i2 * e0_dldgreek_tau_i2_j1 + e0_T_fe1_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe1_i1 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe1_i1 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe1_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i1 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe1_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i1 ) ) ); 
	Y(84) = ( e0_T_fe1_i0 * e0_dldgreek_tau_i0_j2 + e0_T_fe1_i1 * e0_dldgreek_tau_i1_j2 + e0_T_fe1_i2 * e0_dldgreek_tau_i2_j2 + e0_T_fe1_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe1_i2 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe1_i2 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe1_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i2 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe1_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i2 ) ) ); 
	Y(85) = ( e0_T_fe1_i0 * e0_dldgreek_tau_i0_j3 + e0_T_fe1_i1 * e0_dldgreek_tau_i1_j3 + e0_T_fe1_i2 * e0_dldgreek_tau_i2_j3 + e0_T_fe1_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe1_i3 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe1_i3 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe1_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe1_i3 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe1_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe1_i3 ) ) ); 
	Y(86) = ( e0_T_fe2_i0 * e0_dldgreek_tau_i0_j1 + e0_T_fe2_i1 * e0_dldgreek_tau_i1_j1 + e0_T_fe2_i2 * e0_dldgreek_tau_i2_j1 + e0_T_fe2_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe2_i1 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe2_i1 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe2_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i1 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe2_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i1 ) ) ); 
	Y(87) = ( e0_T_fe2_i0 * e0_dldgreek_tau_i0_j2 + e0_T_fe2_i1 * e0_dldgreek_tau_i1_j2 + e0_T_fe2_i2 * e0_dldgreek_tau_i2_j2 + e0_T_fe2_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe2_i2 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe2_i2 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe2_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i2 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe2_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i2 ) ) ); 
	Y(88) = ( e0_T_fe2_i0 * e0_dldgreek_tau_i0_j3 + e0_T_fe2_i1 * e0_dldgreek_tau_i1_j3 + e0_T_fe2_i2 * e0_dldgreek_tau_i2_j3 + e0_T_fe2_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe2_i3 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe2_i3 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe2_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe2_i3 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe2_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe2_i3 ) ) ); 
	Y(89) = ( e0_T_fe3_i0 * e0_dldgreek_tau_i0_j1 + e0_T_fe3_i1 * e0_dldgreek_tau_i1_j1 + e0_T_fe3_i2 * e0_dldgreek_tau_i2_j1 + e0_T_fe3_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe3_i1 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe3_i1 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe3_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i1 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe3_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i1 ) ) ); 
	Y(90) = ( e0_T_fe3_i0 * e0_dldgreek_tau_i0_j2 + e0_T_fe3_i1 * e0_dldgreek_tau_i1_j2 + e0_T_fe3_i2 * e0_dldgreek_tau_i2_j2 + e0_T_fe3_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe3_i2 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe3_i2 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe3_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i2 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe3_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i2 ) ) ); 
	Y(91) = ( e0_T_fe3_i0 * e0_dldgreek_tau_i0_j3 + e0_T_fe3_i1 * e0_dldgreek_tau_i1_j3 + e0_T_fe3_i2 * e0_dldgreek_tau_i2_j3 + e0_T_fe3_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe3_i3 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe3_i3 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe3_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe3_i3 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe3_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe3_i3 ) ) ); 
	Y(92) = ( e0_T_fe4_i0 * e0_dldgreek_tau_i0_j1 + e0_T_fe4_i1 * e0_dldgreek_tau_i1_j1 + e0_T_fe4_i2 * e0_dldgreek_tau_i2_j1 + e0_T_fe4_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe4_i1 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe4_i1 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe4_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i1 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe4_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i1 ) ) ); 
	Y(93) = ( e0_T_fe4_i0 * e0_dldgreek_tau_i0_j2 + e0_T_fe4_i1 * e0_dldgreek_tau_i1_j2 + e0_T_fe4_i2 * e0_dldgreek_tau_i2_j2 + e0_T_fe4_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe4_i2 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe4_i2 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe4_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i2 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe4_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i2 ) ) ); 
	Y(94) = ( e0_T_fe4_i0 * e0_dldgreek_tau_i0_j3 + e0_T_fe4_i1 * e0_dldgreek_tau_i1_j3 + e0_T_fe4_i2 * e0_dldgreek_tau_i2_j3 + e0_T_fe4_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe4_i3 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe4_i3 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe4_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe4_i3 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe4_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe4_i3 ) ) ); 
	Y(95) = ( e0_T_fe5_i0 * e0_dldgreek_tau_i0_j1 + e0_T_fe5_i1 * e0_dldgreek_tau_i1_j1 + e0_T_fe5_i2 * e0_dldgreek_tau_i2_j1 + e0_T_fe5_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe5_i1 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe5_i1 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe5_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i1 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe5_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i1 ) ) ); 
	Y(96) = ( e0_T_fe5_i0 * e0_dldgreek_tau_i0_j2 + e0_T_fe5_i1 * e0_dldgreek_tau_i1_j2 + e0_T_fe5_i2 * e0_dldgreek_tau_i2_j2 + e0_T_fe5_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe5_i2 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe5_i2 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe5_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i2 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe5_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i2 ) ) ); 
	Y(97) = ( e0_T_fe5_i0 * e0_dldgreek_tau_i0_j3 + e0_T_fe5_i1 * e0_dldgreek_tau_i1_j3 + e0_T_fe5_i2 * e0_dldgreek_tau_i2_j3 + e0_T_fe5_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe5_i3 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe5_i3 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe5_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe5_i3 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe5_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe5_i3 ) ) ); 
	Y(98) = ( e0_T_fe6_i0 * e0_dldgreek_tau_i0_j1 + e0_T_fe6_i1 * e0_dldgreek_tau_i1_j1 + e0_T_fe6_i2 * e0_dldgreek_tau_i2_j1 + e0_T_fe6_i3 * e0_dldgreek_tau_i3_j1 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe6_i1 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe6_i1 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe6_i1 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i1 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe6_i1 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i1 ) ) ); 
	Y(99) = ( e0_T_fe6_i0 * e0_dldgreek_tau_i0_j2 + e0_T_fe6_i1 * e0_dldgreek_tau_i1_j2 + e0_T_fe6_i2 * e0_dldgreek_tau_i2_j2 + e0_T_fe6_i3 * e0_dldgreek_tau_i3_j2 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe6_i2 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe6_i2 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe6_i2 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i2 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe6_i2 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i2 ) ) ); 
	Y(100) = ( e0_T_fe6_i0 * e0_dldgreek_tau_i0_j3 + e0_T_fe6_i1 * e0_dldgreek_tau_i1_j3 + e0_T_fe6_i2 * e0_dldgreek_tau_i2_j3 + e0_T_fe6_i3 * e0_dldgreek_tau_i3_j3 ) * ( 1.0 )/( e0_greek_Deltat ) - ( ( e0_F )/( e0_V ) * ( e0_T_In - e0_T_fe6_i3 ) + ( e0_UA )/( e0_greek_rho * e0_cp * e0_V ) * ( e0_T_Cool - e0_T_fe6_i3 ) + (  - e0_H_A )/( e0_greek_rho * e0_cp ) * e0_k_A * e0_c_A_fe6_i3 * exp( (  - e0_E_A )/( e0_R * e0_T_fe6_i3 ) ) + (  - e0_H_B )/( e0_greek_rho * e0_cp ) * e0_k_B * e0_c_B_fe6_i3 * exp( (  - e0_E_B )/( e0_R * e0_T_fe6_i3 ) ) ); 



end

function[] = displayResults(X_ITER)

	% print variable values to display 
	disp(['e0_c_A_fe1_i0  ', num2str(X_ITER(1))]);
	disp(['e0_c_A_fe1_i3  ', num2str(X_ITER(2))]);
	disp(['e0_c_A_fe2_i0  ', num2str(X_ITER(3))]);
	disp(['e0_c_A_fe2_i3  ', num2str(X_ITER(4))]);
	disp(['e0_c_A_fe3_i0  ', num2str(X_ITER(5))]);
	disp(['e0_c_A_fe3_i3  ', num2str(X_ITER(6))]);
	disp(['e0_c_A_fe4_i0  ', num2str(X_ITER(7))]);
	disp(['e0_c_A_fe4_i3  ', num2str(X_ITER(8))]);
	disp(['e0_c_A_fe5_i0  ', num2str(X_ITER(9))]);
	disp(['e0_c_A_fe5_i3  ', num2str(X_ITER(10))]);
	disp(['e0_c_A_fe6_i0  ', num2str(X_ITER(11))]);
	disp(['e0_c_A_fe6_i3  ', num2str(X_ITER(12))]);
	disp(['e0_c_A_fe7_i0  ', num2str(X_ITER(13))]);
	disp(['e0_T_fe1_i1  ', num2str(X_ITER(14))]);
	disp(['e0_T_fe1_i2  ', num2str(X_ITER(15))]);
	disp(['e0_T_fe1_i3  ', num2str(X_ITER(16))]);
	disp(['e0_T_fe2_i1  ', num2str(X_ITER(17))]);
	disp(['e0_T_fe2_i2  ', num2str(X_ITER(18))]);
	disp(['e0_T_fe2_i3  ', num2str(X_ITER(19))]);
	disp(['e0_T_fe3_i1  ', num2str(X_ITER(20))]);
	disp(['e0_T_fe3_i2  ', num2str(X_ITER(21))]);
	disp(['e0_T_fe3_i3  ', num2str(X_ITER(22))]);
	disp(['e0_T_fe4_i1  ', num2str(X_ITER(23))]);
	disp(['e0_T_fe4_i2  ', num2str(X_ITER(24))]);
	disp(['e0_T_fe4_i3  ', num2str(X_ITER(25))]);
	disp(['e0_T_fe5_i1  ', num2str(X_ITER(26))]);
	disp(['e0_T_fe5_i2  ', num2str(X_ITER(27))]);
	disp(['e0_T_fe5_i3  ', num2str(X_ITER(28))]);
	disp(['e0_T_fe6_i1  ', num2str(X_ITER(29))]);
	disp(['e0_T_fe6_i2  ', num2str(X_ITER(30))]);
	disp(['e0_T_fe6_i3  ', num2str(X_ITER(31))]);
	disp(['e0_c_A_fe1_i1  ', num2str(X_ITER(32))]);
	disp(['e0_c_A_fe1_i2  ', num2str(X_ITER(33))]);
	disp(['e0_c_A_fe2_i1  ', num2str(X_ITER(34))]);
	disp(['e0_c_A_fe2_i2  ', num2str(X_ITER(35))]);
	disp(['e0_c_A_fe3_i1  ', num2str(X_ITER(36))]);
	disp(['e0_c_A_fe3_i2  ', num2str(X_ITER(37))]);
	disp(['e0_c_A_fe4_i1  ', num2str(X_ITER(38))]);
	disp(['e0_c_A_fe4_i2  ', num2str(X_ITER(39))]);
	disp(['e0_c_A_fe5_i1  ', num2str(X_ITER(40))]);
	disp(['e0_c_A_fe5_i2  ', num2str(X_ITER(41))]);
	disp(['e0_c_A_fe6_i1  ', num2str(X_ITER(42))]);
	disp(['e0_c_A_fe6_i2  ', num2str(X_ITER(43))]);
	disp(['e0_c_B_fe1_i1  ', num2str(X_ITER(44))]);
	disp(['e0_c_B_fe1_i2  ', num2str(X_ITER(45))]);
	disp(['e0_c_B_fe1_i3  ', num2str(X_ITER(46))]);
	disp(['e0_c_B_fe2_i1  ', num2str(X_ITER(47))]);
	disp(['e0_c_B_fe2_i2  ', num2str(X_ITER(48))]);
	disp(['e0_c_B_fe2_i3  ', num2str(X_ITER(49))]);
	disp(['e0_c_B_fe3_i1  ', num2str(X_ITER(50))]);
	disp(['e0_c_B_fe3_i2  ', num2str(X_ITER(51))]);
	disp(['e0_c_B_fe3_i3  ', num2str(X_ITER(52))]);
	disp(['e0_c_B_fe4_i1  ', num2str(X_ITER(53))]);
	disp(['e0_c_B_fe4_i2  ', num2str(X_ITER(54))]);
	disp(['e0_c_B_fe4_i3  ', num2str(X_ITER(55))]);
	disp(['e0_c_B_fe5_i1  ', num2str(X_ITER(56))]);
	disp(['e0_c_B_fe5_i2  ', num2str(X_ITER(57))]);
	disp(['e0_c_B_fe5_i3  ', num2str(X_ITER(58))]);
	disp(['e0_c_B_fe6_i1  ', num2str(X_ITER(59))]);
	disp(['e0_c_B_fe6_i2  ', num2str(X_ITER(60))]);
	disp(['e0_c_B_fe6_i3  ', num2str(X_ITER(61))]);
	disp(['e0_c_B_fe1_i0  ', num2str(X_ITER(62))]);
	disp(['e0_c_B_fe2_i0  ', num2str(X_ITER(63))]);
	disp(['e0_c_B_fe3_i0  ', num2str(X_ITER(64))]);
	disp(['e0_c_B_fe4_i0  ', num2str(X_ITER(65))]);
	disp(['e0_c_B_fe5_i0  ', num2str(X_ITER(66))]);
	disp(['e0_c_B_fe6_i0  ', num2str(X_ITER(67))]);
	disp(['e0_c_B_fe7_i0  ', num2str(X_ITER(68))]);
	disp(['e0_c_C_fe1_i0  ', num2str(X_ITER(69))]);
	disp(['e0_c_C_fe1_i3  ', num2str(X_ITER(70))]);
	disp(['e0_c_C_fe2_i0  ', num2str(X_ITER(71))]);
	disp(['e0_c_C_fe2_i3  ', num2str(X_ITER(72))]);
	disp(['e0_c_C_fe3_i0  ', num2str(X_ITER(73))]);
	disp(['e0_c_C_fe3_i3  ', num2str(X_ITER(74))]);
	disp(['e0_c_C_fe4_i0  ', num2str(X_ITER(75))]);
	disp(['e0_c_C_fe4_i3  ', num2str(X_ITER(76))]);
	disp(['e0_c_C_fe5_i0  ', num2str(X_ITER(77))]);
	disp(['e0_c_C_fe5_i3  ', num2str(X_ITER(78))]);
	disp(['e0_c_C_fe6_i0  ', num2str(X_ITER(79))]);
	disp(['e0_c_C_fe6_i3  ', num2str(X_ITER(80))]);
	disp(['e0_c_C_fe7_i0  ', num2str(X_ITER(81))]);
	disp(['e0_c_C_fe1_i1  ', num2str(X_ITER(82))]);
	disp(['e0_c_C_fe1_i2  ', num2str(X_ITER(83))]);
	disp(['e0_c_C_fe2_i1  ', num2str(X_ITER(84))]);
	disp(['e0_c_C_fe2_i2  ', num2str(X_ITER(85))]);
	disp(['e0_c_C_fe3_i1  ', num2str(X_ITER(86))]);
	disp(['e0_c_C_fe3_i2  ', num2str(X_ITER(87))]);
	disp(['e0_c_C_fe4_i1  ', num2str(X_ITER(88))]);
	disp(['e0_c_C_fe4_i2  ', num2str(X_ITER(89))]);
	disp(['e0_c_C_fe5_i1  ', num2str(X_ITER(90))]);
	disp(['e0_c_C_fe5_i2  ', num2str(X_ITER(91))]);
	disp(['e0_c_C_fe6_i1  ', num2str(X_ITER(92))]);
	disp(['e0_c_C_fe6_i2  ', num2str(X_ITER(93))]);
	disp(['e0_T_fe1_i0  ', num2str(X_ITER(94))]);
	disp(['e0_T_fe2_i0  ', num2str(X_ITER(95))]);
	disp(['e0_T_fe3_i0  ', num2str(X_ITER(96))]);
	disp(['e0_T_fe4_i0  ', num2str(X_ITER(97))]);
	disp(['e0_T_fe5_i0  ', num2str(X_ITER(98))]);
	disp(['e0_T_fe6_i0  ', num2str(X_ITER(99))]);
	disp(['e0_T_fe7_i0  ', num2str(X_ITER(100))]);


end

