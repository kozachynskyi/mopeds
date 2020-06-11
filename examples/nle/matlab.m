%*********************************************************
% The namespaces have been normalized. The following
% table shows the attribuation. 
% Normalized Name --> Full Name ---> User-defined Name
% =================================== 
% e0 --> e[0]126726 --> 
%*********************************************************

%*********************************************************
% The variables are named according to the notation
% provided in the Mosaic model.
% 
% The variable names can be read as follows:
% ==========================================
% 	e0_T_j#
% 		T: temperature
% 		Indices
% 			j: Number of streams
% 	 
% 	e0_A_r#_par#
% 		A: parameter
% 		Indices
% 			r: Number of Rections
% 			par: number of coefficients
% 	 
% 	e0_x_i#_j#
% 		x: mole fraction
% 		Indices
% 			i: Number of components
% 			j: Number of streams
% 	 
% 	e0_p_j#
% 		p: pressure
% 		Indices
% 			j: Number of streams
% 	 
% 	e0_V
% 		V: Volume
% 	 
% 	e0_greek_rho_r#
% 		&rho;: kg catalyst /m³
% 		Indices
% 			r: Number of Rections
% 	 
% 	e0_F_j#
% 		F: Flow
% 		Indices
% 			j: Number of streams
% 	 
% 	e0_K_r#
% 		K: equilibrium constant
% 		Indices
% 			r: Number of Rections
% 	 
% 	e0_X_i#
% 		X: conversion
% 		Indices
% 			i: Number of components
% 	 
% 	e0_Y_i#
% 		Y: yield
% 		Indices
% 			i: Number of components
% 	 
% 	e0_h_j#
% 		h: specific enthalpy
% 		Indices
% 			j: Number of streams
% 	 
% 	e0_h_o_j#_i#
% 		h: specific enthalpy
% 		Superscripts
% 			o: pure component
% 		Indices
% 			j: Number of streams
% 			i: Number of components
% 	 
% 	e0_p_i#_j#
% 		p: pressure
% 		Indices
% 			i: Number of components
% 			j: Number of streams
% 	 
% 	e0_HU
% 		HU: Hold up
% 	 
% 	e0_U
% 		U: internal energy
% 	 
% 	e0_u_j#
% 		u: specific internal energy
% 		Indices
% 			j: Number of streams
% 	 
% 	e0_HU_i#
% 		HU: Hold up
% 		Indices
% 			i: Number of components
% 	 
% 	e0_v_j#
% 		v: specific volume
% 		Indices
% 			j: Number of streams
% 	 
% 	e0_K_i#
% 		K: equilibrium constant
% 		Indices
% 			i: Number of components
% 	 
% 	e0_k_r#
% 		k: arrhenius factor
% 		Indices
% 			r: Number of Rections
% 	 
% 	e0_r_r#
% 		r: reaction rate mol/s
% 		Indices
% 			r: Number of Rections
% 	 
% 	e0_k_par#
% 		k: arrhenius factor
% 		Indices
% 			par: number of coefficients
% 	 
% 	e0_Q
% 		Q: Heat Stream
% 	 
% 	e0_greek_nu_r#_i#
% 		&nu;: stochiometric coefficient
% 		Indices
% 			r: Number of Rections
% 			i: Number of components
% 	 
% 	e0_A_par#_i#
% 		A: parameter
% 		Indices
% 			par: number of coefficients
% 			i: Number of components
% 	 
% 	e0_T_f
% 		T: temperature
% 		Superscripts
% 			f: formation
% 	 
% 	e0_h_f_i#
% 		h: specific enthalpy
% 		Superscripts
% 			f: formation
% 		Indices
% 			i: Number of components
% 	 
% 	e0_E_par#_r#
% 		E: activation energy
% 		Indices
% 			par: number of coefficients
% 			r: Number of Rections
% 	 
% 	e0_A_r#_i#
% 		A: parameter
% 		Indices
% 			r: Number of Rections
% 			i: Number of components
% 	 
% 	e0_E_r#_i#
% 		E: activation energy
% 		Indices
% 			r: Number of Rections
% 			i: Number of components
% 	 
% 	e0_A_r#
% 		A: parameter
% 		Indices
% 			r: Number of Rections
% 	 
% 	e0_E_r#
% 		E: activation energy
% 		Indices
% 			r: Number of Rections
% 	 
% 	e0_A_eq_r#_par#
% 		A: parameter
% 		Superscripts
% 			eq: equilibrium
% 		Indices
% 			r: Number of Rections
% 			par: number of coefficients
% 	 
% 	e0_R
% 		R: universal gas constant
% 	 
%*********************************************************

function[ROOTS]=solveEquationSystem()

	% load variable init values 
	X_ITER(1) = 0.52397;  	% e0_x_i2_j2  
	X_ITER(2) = 0.21455;  	% e0_x_i1_j2  
	X_ITER(3) = 0.0074497;  	% e0_HU_i3  
	X_ITER(4) = 0.063022;  	% e0_HU_i2  
	X_ITER(5) = 23.0699;  	% e0_K_r3  
	X_ITER(6) = 0.025805;  	% e0_HU_i1  
	X_ITER(7) = -0.083762;  	% e0_u_j2  
	X_ITER(8) = -10.0748;  	% e0_U  
	X_ITER(9) = 0.12028;  	% e0_HU  
	X_ITER(10) = 0.7498;  	% e0_p_i7_j2  
	X_ITER(11) = 0.31874;  	% e0_p_i6_j2  
	X_ITER(12) = 8.3803;  	% e0_p_i5_j2  
	X_ITER(13) = 0.52869;  	% e0_p_i4_j2  
	X_ITER(14) = 3.0968;  	% e0_p_i3_j2  
	X_ITER(15) = 26.1983;  	% e0_p_i2_j2  
	X_ITER(16) = 10.7273;  	% e0_p_i1_j2  
	X_ITER(17) = -168023.1327725076;  	% e0_h_o_j2_i7  
	X_ITER(18) = -184427.01699715087;  	% e0_h_o_j1_i7  
	X_ITER(19) = -190537.8047525016;  	% e0_h_o_j2_i6  
	X_ITER(20) = -201159.24396779988;  	% e0_h_o_j1_i6  
	X_ITER(21) = 4.6492;  	% e0_Q  
	X_ITER(22) = -385214.2159451514;  	% e0_h_o_j2_i5  
	X_ITER(23) = -393705.0747228285;  	% e0_h_o_j1_i5  
	X_ITER(24) = -104597.15130755084;  	% e0_h_o_j2_i4  
	X_ITER(25) = -110675.6190793542;  	% e0_h_o_j1_i4  
	X_ITER(26) = -234887.44250385024;  	% e0_h_o_j2_i3  
	X_ITER(27) = -0.9322;  	% e0_F_j2  
	X_ITER(28) = -241998.63319940824;  	% e0_h_o_j1_i3  
	X_ITER(29) = 5877.734015746798;  	% e0_h_o_j2_i2  
	X_ITER(30) = 0.0098568;  	% e0_r_r1  
	X_ITER(31) = -144.29466839718256;  	% e0_h_o_j1_i2  
	X_ITER(32) = 0.6905973165554351;  	% e0_k_par5  
	X_ITER(33) = 5910.160589771112;  	% e0_h_o_j2_i1  
	X_ITER(34) = 0.033901;  	% e0_r_r2  
	X_ITER(35) = -145.74383548430316;  	% e0_h_o_j1_i1  
	X_ITER(36) = 613.6897535075899;  	% e0_k_par4  
	X_ITER(37) = 31.24090746419173;  	% e0_k_par3  
	X_ITER(38) = 3453.38;  	% e0_k_par2  
	X_ITER(39) = 16154.96060463936;  	% e0_k_par1  
	X_ITER(40) = -0.079605;  	% e0_h_j2  
	X_ITER(41) = 0.0076384;  	% e0_K_r1  
	X_ITER(42) = -0.078857;  	% e0_h_j1  
	X_ITER(43) = 0.014996;  	% e0_x_i7_j2  
	X_ITER(44) = 4.2463E-17;  	% e0_x_i7_j1  
	X_ITER(45) = 0.0063748;  	% e0_x_i6_j2  
	X_ITER(46) = 0.013979;  	% e0_r_r3  
	X_ITER(47) = 3.4674E-5;  	% e0_K_r2  
	X_ITER(48) = 0.15069620388555666;  	% e0_k_r3  
	X_ITER(49) = 0.16761;  	% e0_x_i5_j2  
	X_ITER(50) = 2.132122117309319E13;  	% e0_K_i6  
	X_ITER(51) = 2.2476341621984563E9;  	% e0_K_i3  
	X_ITER(52) = 0.010574;  	% e0_x_i4_j2  
	X_ITER(53) = 8.314E-4;  	% e0_v_j2  
	X_ITER(54) = 0.0018037;  	% e0_HU_i7  
	X_ITER(55) = 7.6676E-4;  	% e0_HU_i6  
	X_ITER(56) = 0.02016;  	% e0_HU_i5  
	X_ITER(57) = 0.0012718;  	% e0_HU_i4  
	X_ITER(58) = 0.18593;  	% e0_X_i2  
	X_ITER(59) = 0.069896;  	% e0_Y_i7  
	X_ITER(60) = 0.061937;  	% e0_x_i3_j2  
	X_ITER(61) = 0.21879;  	% e0_X_i5  

	% load parameters 
	PARAMS(1) = 3.8E-7;  	% e0_A_r3_par4 
	PARAMS(2) = -65610.0;  	% e0_A_r3_par5 
	PARAMS(3) = 0.0;  	% e0_greek_nu_r1_i1 
	PARAMS(4) = -26.64;  	% e0_A_r3_par6 
	PARAMS(5) = 0.0;  	% e0_greek_nu_r2_i1 
	PARAMS(6) = 0.0;  	% e0_greek_nu_r3_i1 
	PARAMS(7) = 0.2;  	% e0_x_i1_j1 
	PARAMS(8) = -1.0;  	% e0_greek_nu_r1_i2 
	PARAMS(9) = -3.0;  	% e0_greek_nu_r2_i2 
	PARAMS(10) = 0.6;  	% e0_x_i2_j1 
	PARAMS(11) = 0.0;  	% e0_greek_nu_r3_i2 
	PARAMS(12) = 0.0;  	% e0_x_i3_j1 
	PARAMS(13) = 3.7E-5;  	% e0_A_par3_i1 
	PARAMS(14) = -0.0169;  	% e0_A_par4_i1 
	PARAMS(15) = 0.0;  	% e0_greek_nu_r2_i7 
	PARAMS(16) = 1.0;  	% e0_greek_nu_r3_i7 
	PARAMS(17) = 5.45E-12;  	% e0_A_par1_i1 
	PARAMS(18) = -2.44E-8;  	% e0_A_par2_i1 
	PARAMS(19) = 1.0;  	% e0_greek_nu_r1_i4 
	PARAMS(20) = 0.0;  	% e0_greek_nu_r2_i4 
	PARAMS(21) = 0.0;  	% e0_greek_nu_r3_i4 
	PARAMS(22) = -1.0;  	% e0_greek_nu_r1_i5 
	PARAMS(23) = -1.0;  	% e0_greek_nu_r2_i5 
	PARAMS(24) = 0.0;  	% e0_greek_nu_r3_i5 
	PARAMS(25) = 0.0;  	% e0_greek_nu_r1_i6 
	PARAMS(26) = 1.0;  	% e0_greek_nu_r2_i6 
	PARAMS(27) = 1.0;  	% e0_greek_nu_r1_i3 
	PARAMS(28) = 1.0;  	% e0_greek_nu_r2_i3 
	PARAMS(29) = 1.0;  	% e0_greek_nu_r3_i3 
	PARAMS(30) = -2.50993E-7;  	% e0_A_par2_i3 
	PARAMS(31) = 2.2E-4;  	% e0_A_par3_i3 
	PARAMS(32) = 3.707;  	% e0_A_r3_par2 
	PARAMS(33) = -0.002783;  	% e0_A_r3_par3 
	PARAMS(34) = 4019.0;  	% e0_A_r3_par1 
	PARAMS(35) = -241830.0;  	% e0_h_f_i3 
	PARAMS(36) = -1.54836E-11;  	% e0_A_par1_i4 
	PARAMS(37) = 8.314;  	% e0_R 
	PARAMS(38) = -0.074;  	% e0_A_par4_i3 
	PARAMS(39) = 42.04061276;  	% e0_A_par5_i3 
	PARAMS(40) = -0.00322;  	% e0_A_par4_i4 
	PARAMS(41) = 29.59614774;  	% e0_A_par5_i4 
	PARAMS(42) = 2.10847E-8;  	% e0_A_par2_i4 
	PARAMS(43) = 6.07E-7;  	% e0_A_par3_i4 
	PARAMS(44) = 0.0;  	% e0_x_i4_j1 
	PARAMS(45) = 0.2;  	% e0_x_i5_j1 
	PARAMS(46) = 1.0E-16;  	% e0_x_i6_j1 
	PARAMS(47) = 31.5;  	% e0_A_par5_i1 
	PARAMS(48) = 500.0;  	% e0_T_j2 
	PARAMS(49) = 1.07368E-10;  	% e0_A_par1_i3 
	PARAMS(50) = -7.52E-12;  	% e0_A_par1_i2 
	PARAMS(51) = 2.7E-8;  	% e0_A_par2_i2 
	PARAMS(52) = 298.15;  	% e0_T_f 
	PARAMS(53) = 0.0;  	% e0_h_f_i1 
	PARAMS(54) = 26.2;  	% e0_A_par5_i2 
	PARAMS(55) = 0.0;  	% e0_h_f_i2 
	PARAMS(56) = -3.17E-5;  	% e0_A_par3_i2 
	PARAMS(57) = 0.0162;  	% e0_A_par4_i2 
	PARAMS(58) = -2073.0;  	% e0_A_eq_r1_par1 
	PARAMS(59) = 2.029;  	% e0_A_eq_r1_par2 
	PARAMS(60) = 0.1;  	% e0_V 
	PARAMS(61) = -200940.0;  	% e0_h_f_i6 
	PARAMS(62) = 3066.0;  	% e0_A_eq_r2_par1 
	PARAMS(63) = 3.87561E-10;  	% e0_A_par1_i7 
	PARAMS(64) = -10.592;  	% e0_A_eq_r2_par2 
	PARAMS(65) = -0.121;  	% e0_A_par4_i6 
	PARAMS(66) = 45.78156724;  	% e0_A_par5_i6 
	PARAMS(67) = -0.0526;  	% e0_A_par4_i7 
	PARAMS(68) = 45.14491836;  	% e0_A_par5_i7 
	PARAMS(69) = -8.56412E-7;  	% e0_A_par2_i7 
	PARAMS(70) = 6.29E-4;  	% e0_A_par3_i7 
	PARAMS(71) = 50.0;  	% e0_p_j2 
	PARAMS(72) = -184100.0;  	% e0_h_f_i7 
	PARAMS(73) = 1.07;  	% e0_A_par1_r1 
	PARAMS(74) = 1.93148E-11;  	% e0_A_par1_i5 
	PARAMS(75) = -1.32293E-8;  	% e0_A_par2_i5 
	PARAMS(76) = -110530.0;  	% e0_h_f_i4 
	PARAMS(77) = 19.6353642;  	% e0_A_par5_i5 
	PARAMS(78) = -393520.0;  	% e0_h_f_i5 
	PARAMS(79) = -4.18E-5;  	% e0_A_par3_i5 
	PARAMS(80) = 0.0718;  	% e0_A_par4_i5 
	PARAMS(81) = 5.57E-4;  	% e0_A_par3_i6 
	PARAMS(82) = 2.65055E-10;  	% e0_A_par1_i6 
	PARAMS(83) = -6.5285E-7;  	% e0_A_par2_i6 
	PARAMS(84) = 92000.0;  	% e0_E_r3_i3 
	PARAMS(85) = 223.2;  	% e0_A_r3_i6 
	PARAMS(86) = 1.0;  	% e0_F_j1 
	PARAMS(87) = -55060.0;  	% e0_E_r3 
	PARAMS(88) = 105100.0;  	% e0_E_r3_i6 
	PARAMS(89) = 85190.0;  	% e0_A_r3 
	PARAMS(90) = 293.15;  	% e0_T_j1 
	PARAMS(91) = -2.0;  	% e0_greek_nu_r3_i6 
	PARAMS(92) = 0.0;  	% e0_greek_nu_r1_i7 
	PARAMS(93) = 1775.0;  	% e0_greek_rho_r1 
	PARAMS(94) = 100.0;  	% e0_greek_rho_r3 
	PARAMS(95) = 40000.0;  	% e0_E_par1_r1 
	PARAMS(96) = 0.499;  	% e0_A_par3_r1 
	PARAMS(97) = 17197.0;  	% e0_E_par3_r1 
	PARAMS(98) = 3453.38;  	% e0_A_par2_r1 
	PARAMS(99) = 0.0;  	% e0_E_par2_r1 
	PARAMS(100) = 1.22E10;  	% e0_A_par5_r1 
	PARAMS(101) = -98084.0;  	% e0_E_par5_r1 
	PARAMS(102) = 6.62E-11;  	% e0_A_par4_r1 
	PARAMS(103) = 124119.0;  	% e0_E_par4_r1 
	PARAMS(104) = 0.5498;  	% e0_A_r3_i3 

	options = optimset('MaxIter',1000,'TolFun',1e-6,'Display','Iter');
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
	e0_x_i2_j2 = X_ITER(1); 
	e0_x_i1_j2 = X_ITER(2); 
	e0_HU_i3 = X_ITER(3); 
	e0_HU_i2 = X_ITER(4); 
	e0_K_r3 = X_ITER(5); 
	e0_HU_i1 = X_ITER(6); 
	e0_u_j2 = X_ITER(7); 
	e0_U = X_ITER(8); 
	e0_HU = X_ITER(9); 
	e0_p_i7_j2 = X_ITER(10); 
	e0_p_i6_j2 = X_ITER(11); 
	e0_p_i5_j2 = X_ITER(12); 
	e0_p_i4_j2 = X_ITER(13); 
	e0_p_i3_j2 = X_ITER(14); 
	e0_p_i2_j2 = X_ITER(15); 
	e0_p_i1_j2 = X_ITER(16); 
	e0_h_o_j2_i7 = X_ITER(17); 
	e0_h_o_j1_i7 = X_ITER(18); 
	e0_h_o_j2_i6 = X_ITER(19); 
	e0_h_o_j1_i6 = X_ITER(20); 
	e0_Q = X_ITER(21); 
	e0_h_o_j2_i5 = X_ITER(22); 
	e0_h_o_j1_i5 = X_ITER(23); 
	e0_h_o_j2_i4 = X_ITER(24); 
	e0_h_o_j1_i4 = X_ITER(25); 
	e0_h_o_j2_i3 = X_ITER(26); 
	e0_F_j2 = X_ITER(27); 
	e0_h_o_j1_i3 = X_ITER(28); 
	e0_h_o_j2_i2 = X_ITER(29); 
	e0_r_r1 = X_ITER(30); 
	e0_h_o_j1_i2 = X_ITER(31); 
	e0_k_par5 = X_ITER(32); 
	e0_h_o_j2_i1 = X_ITER(33); 
	e0_r_r2 = X_ITER(34); 
	e0_h_o_j1_i1 = X_ITER(35); 
	e0_k_par4 = X_ITER(36); 
	e0_k_par3 = X_ITER(37); 
	e0_k_par2 = X_ITER(38); 
	e0_k_par1 = X_ITER(39); 
	e0_h_j2 = X_ITER(40); 
	e0_K_r1 = X_ITER(41); 
	e0_h_j1 = X_ITER(42); 
	e0_x_i7_j2 = X_ITER(43); 
	e0_x_i7_j1 = X_ITER(44); 
	e0_x_i6_j2 = X_ITER(45); 
	e0_r_r3 = X_ITER(46); 
	e0_K_r2 = X_ITER(47); 
	e0_k_r3 = X_ITER(48); 
	e0_x_i5_j2 = X_ITER(49); 
	e0_K_i6 = X_ITER(50); 
	e0_K_i3 = X_ITER(51); 
	e0_x_i4_j2 = X_ITER(52); 
	e0_v_j2 = X_ITER(53); 
	e0_HU_i7 = X_ITER(54); 
	e0_HU_i6 = X_ITER(55); 
	e0_HU_i5 = X_ITER(56); 
	e0_HU_i4 = X_ITER(57); 
	e0_X_i2 = X_ITER(58); 
	e0_Y_i7 = X_ITER(59); 
	e0_x_i3_j2 = X_ITER(60); 
	e0_X_i5 = X_ITER(61); 

	% read out parameters 
	e0_A_r3_par4 = PARAMS(1); 
	e0_A_r3_par5 = PARAMS(2); 
	e0_greek_nu_r1_i1 = PARAMS(3); 
	e0_A_r3_par6 = PARAMS(4); 
	e0_greek_nu_r2_i1 = PARAMS(5); 
	e0_greek_nu_r3_i1 = PARAMS(6); 
	e0_x_i1_j1 = PARAMS(7); 
	e0_greek_nu_r1_i2 = PARAMS(8); 
	e0_greek_nu_r2_i2 = PARAMS(9); 
	e0_x_i2_j1 = PARAMS(10); 
	e0_greek_nu_r3_i2 = PARAMS(11); 
	e0_x_i3_j1 = PARAMS(12); 
	e0_A_par3_i1 = PARAMS(13); 
	e0_A_par4_i1 = PARAMS(14); 
	e0_greek_nu_r2_i7 = PARAMS(15); 
	e0_greek_nu_r3_i7 = PARAMS(16); 
	e0_A_par1_i1 = PARAMS(17); 
	e0_A_par2_i1 = PARAMS(18); 
	e0_greek_nu_r1_i4 = PARAMS(19); 
	e0_greek_nu_r2_i4 = PARAMS(20); 
	e0_greek_nu_r3_i4 = PARAMS(21); 
	e0_greek_nu_r1_i5 = PARAMS(22); 
	e0_greek_nu_r2_i5 = PARAMS(23); 
	e0_greek_nu_r3_i5 = PARAMS(24); 
	e0_greek_nu_r1_i6 = PARAMS(25); 
	e0_greek_nu_r2_i6 = PARAMS(26); 
	e0_greek_nu_r1_i3 = PARAMS(27); 
	e0_greek_nu_r2_i3 = PARAMS(28); 
	e0_greek_nu_r3_i3 = PARAMS(29); 
	e0_A_par2_i3 = PARAMS(30); 
	e0_A_par3_i3 = PARAMS(31); 
	e0_A_r3_par2 = PARAMS(32); 
	e0_A_r3_par3 = PARAMS(33); 
	e0_A_r3_par1 = PARAMS(34); 
	e0_h_f_i3 = PARAMS(35); 
	e0_A_par1_i4 = PARAMS(36); 
	e0_R = PARAMS(37); 
	e0_A_par4_i3 = PARAMS(38); 
	e0_A_par5_i3 = PARAMS(39); 
	e0_A_par4_i4 = PARAMS(40); 
	e0_A_par5_i4 = PARAMS(41); 
	e0_A_par2_i4 = PARAMS(42); 
	e0_A_par3_i4 = PARAMS(43); 
	e0_x_i4_j1 = PARAMS(44); 
	e0_x_i5_j1 = PARAMS(45); 
	e0_x_i6_j1 = PARAMS(46); 
	e0_A_par5_i1 = PARAMS(47); 
	e0_T_j2 = PARAMS(48); 
	e0_A_par1_i3 = PARAMS(49); 
	e0_A_par1_i2 = PARAMS(50); 
	e0_A_par2_i2 = PARAMS(51); 
	e0_T_f = PARAMS(52); 
	e0_h_f_i1 = PARAMS(53); 
	e0_A_par5_i2 = PARAMS(54); 
	e0_h_f_i2 = PARAMS(55); 
	e0_A_par3_i2 = PARAMS(56); 
	e0_A_par4_i2 = PARAMS(57); 
	e0_A_eq_r1_par1 = PARAMS(58); 
	e0_A_eq_r1_par2 = PARAMS(59); 
	e0_V = PARAMS(60); 
	e0_h_f_i6 = PARAMS(61); 
	e0_A_eq_r2_par1 = PARAMS(62); 
	e0_A_par1_i7 = PARAMS(63); 
	e0_A_eq_r2_par2 = PARAMS(64); 
	e0_A_par4_i6 = PARAMS(65); 
	e0_A_par5_i6 = PARAMS(66); 
	e0_A_par4_i7 = PARAMS(67); 
	e0_A_par5_i7 = PARAMS(68); 
	e0_A_par2_i7 = PARAMS(69); 
	e0_A_par3_i7 = PARAMS(70); 
	e0_p_j2 = PARAMS(71); 
	e0_h_f_i7 = PARAMS(72); 
	e0_A_par1_r1 = PARAMS(73); 
	e0_A_par1_i5 = PARAMS(74); 
	e0_A_par2_i5 = PARAMS(75); 
	e0_h_f_i4 = PARAMS(76); 
	e0_A_par5_i5 = PARAMS(77); 
	e0_h_f_i5 = PARAMS(78); 
	e0_A_par3_i5 = PARAMS(79); 
	e0_A_par4_i5 = PARAMS(80); 
	e0_A_par3_i6 = PARAMS(81); 
	e0_A_par1_i6 = PARAMS(82); 
	e0_A_par2_i6 = PARAMS(83); 
	e0_E_r3_i3 = PARAMS(84); 
	e0_A_r3_i6 = PARAMS(85); 
	e0_F_j1 = PARAMS(86); 
	e0_E_r3 = PARAMS(87); 
	e0_E_r3_i6 = PARAMS(88); 
	e0_A_r3 = PARAMS(89); 
	e0_T_j1 = PARAMS(90); 
	e0_greek_nu_r3_i6 = PARAMS(91); 
	e0_greek_nu_r1_i7 = PARAMS(92); 
	e0_greek_rho_r1 = PARAMS(93); 
	e0_greek_rho_r3 = PARAMS(94); 
	e0_E_par1_r1 = PARAMS(95); 
	e0_A_par3_r1 = PARAMS(96); 
	e0_E_par3_r1 = PARAMS(97); 
	e0_A_par2_r1 = PARAMS(98); 
	e0_E_par2_r1 = PARAMS(99); 
	e0_A_par5_r1 = PARAMS(100); 
	e0_E_par5_r1 = PARAMS(101); 
	e0_A_par4_r1 = PARAMS(102); 
	e0_E_par4_r1 = PARAMS(103); 
	e0_A_r3_i3 = PARAMS(104); 

	% evaluate the function values  
	Y(1) = e0_K_r2 - (power((10.0),(e0_A_eq_r2_par1)/(e0_T_j2) + e0_A_eq_r2_par2)); 
	Y(2) = e0_K_r1 - (power((10.0),(e0_A_eq_r1_par1)/(e0_T_j2) + e0_A_eq_r1_par2)); 
	Y(3) = e0_K_r3 - (exp((e0_A_r3_par1)/(e0_T_j2) + e0_A_r3_par2 * log(e0_T_j2) + e0_A_r3_par3 * e0_T_j2 + e0_A_r3_par4 * power((e0_T_j2),2.0) + (e0_A_r3_par5)/(power((e0_T_j2),3.0)) + e0_A_r3_par6)); 
	Y(4) = 1.0 - ((e0_x_i1_j1 + e0_x_i2_j1 + e0_x_i3_j1 + e0_x_i4_j1 + e0_x_i5_j1 + e0_x_i6_j1 + e0_x_i7_j1)); 
	Y(5) = 1.0 - ((e0_x_i1_j2 + e0_x_i2_j2 + e0_x_i3_j2 + e0_x_i4_j2 + e0_x_i5_j2 + e0_x_i6_j2 + e0_x_i7_j2)); 
	Y(6) = e0_h_j1 - (((e0_x_i1_j1 * e0_h_o_j1_i1 + e0_x_i2_j1 * e0_h_o_j1_i2 + e0_x_i3_j1 * e0_h_o_j1_i3 + e0_x_i4_j1 * e0_h_o_j1_i4 + e0_x_i5_j1 * e0_h_o_j1_i5 + e0_x_i6_j1 * e0_h_o_j1_i6 + e0_x_i7_j1 * e0_h_o_j1_i7))/(1000000.0)); 
	Y(7) = e0_h_j2 - (((e0_x_i1_j2 * e0_h_o_j2_i1 + e0_x_i2_j2 * e0_h_o_j2_i2 + e0_x_i3_j2 * e0_h_o_j2_i3 + e0_x_i4_j2 * e0_h_o_j2_i4 + e0_x_i5_j2 * e0_h_o_j2_i5 + e0_x_i6_j2 * e0_h_o_j2_i6 + e0_x_i7_j2 * e0_h_o_j2_i7))/(1000000.0)); 
	Y(8) = e0_p_i1_j2 - (e0_p_j2 * e0_x_i1_j2); 
	Y(9) = e0_p_i2_j2 - (e0_p_j2 * e0_x_i2_j2); 
	Y(10) = e0_p_i3_j2 - (e0_p_j2 * e0_x_i3_j2); 
	Y(11) = e0_p_i4_j2 - (e0_p_j2 * e0_x_i4_j2); 
	Y(12) = e0_p_i5_j2 - (e0_p_j2 * e0_x_i5_j2); 
	Y(13) = e0_p_i6_j2 - (e0_p_j2 * e0_x_i6_j2); 
	Y(14) = e0_p_i7_j2 - (e0_p_j2 * e0_x_i7_j2); 
	Y(15) = e0_U - (e0_u_j2 * e0_HU * 1000.0); 
	Y(16) = e0_p_j2 * 100.0 * e0_V - (e0_HU * e0_R * e0_T_j2); 
	Y(17) = e0_HU_i1 - (e0_HU * e0_x_i1_j2); 
	Y(18) = e0_HU_i2 - (e0_HU * e0_x_i2_j2); 
	Y(19) = e0_HU_i3 - (e0_HU * e0_x_i3_j2); 
	Y(20) = e0_HU_i4 - (e0_HU * e0_x_i4_j2); 
	Y(21) = e0_HU_i5 - (e0_HU * e0_x_i5_j2); 
	Y(22) = e0_HU_i6 - (e0_HU * e0_x_i6_j2); 
	Y(23) = e0_HU_i7 - (e0_HU * e0_x_i7_j2); 
	Y(24) = e0_v_j2 - ((e0_V)/(e0_HU * 1000.0)); 
	Y(25) = e0_r_r3 - (e0_greek_rho_r3 * e0_k_r3 * power((e0_K_i6),2.0) * (power((e0_p_i6_j2),2.0) - (e0_p_i7_j2 * e0_p_i3_j2)/(e0_K_r3))/(power((1.0 + e0_K_i6 * e0_p_i6_j2 + e0_K_i3 * e0_p_i3_j2),2.0)) * e0_V); 
	Y(26) = e0_r_r2 - (e0_greek_rho_r1 * e0_k_par1 * (e0_p_i5_j2 * e0_p_i2_j2 * (1.0 - (e0_p_i3_j2 * e0_p_i6_j2)/(power((e0_p_i2_j2),3.0) * e0_p_i5_j2 * e0_K_r2)))/(power((1.0 + e0_k_par2 * (e0_p_i3_j2)/(e0_p_i2_j2) + e0_k_par3 * power((e0_p_i2_j2),0.5) + e0_k_par4 * (e0_p_i3_j2)),3.0)) * e0_V); 
	Y(27) = e0_r_r1 - (e0_greek_rho_r1 * e0_k_par5 * (e0_p_i5_j2 * (1.0 - (e0_p_i4_j2 * e0_p_i3_j2)/(e0_p_i2_j2 * e0_p_i5_j2 * e0_K_r1)))/(1.0 + e0_k_par2 * (e0_p_i3_j2)/(e0_p_i2_j2) + e0_k_par3 * power((e0_p_i2_j2),0.5) + e0_k_par4 * (e0_p_i3_j2)) * e0_V); 
	Y(28) = e0_u_j2 * 10.0 - (e0_h_j2 * 10.0 - e0_p_j2 * e0_v_j2); 
	Y(29) = 0.0 - ((e0_F_j1 * 1000.0 * e0_h_j1 + e0_F_j2 * 1000.0 * e0_h_j2) + e0_Q); 
	Y(30) = 0.0 - ((e0_F_j1 * e0_x_i1_j1 + e0_F_j2 * e0_x_i1_j2) + (e0_greek_nu_r1_i1 * e0_r_r1 + e0_greek_nu_r2_i1 * e0_r_r2 + e0_greek_nu_r3_i1 * e0_r_r3)); 
	Y(31) = 0.0 - ((e0_F_j1 * e0_x_i2_j1 + e0_F_j2 * e0_x_i2_j2) + (e0_greek_nu_r1_i2 * e0_r_r1 + e0_greek_nu_r2_i2 * e0_r_r2 + e0_greek_nu_r3_i2 * e0_r_r3)); 
	Y(32) = 0.0 - ((e0_F_j1 * e0_x_i3_j1 + e0_F_j2 * e0_x_i3_j2) + (e0_greek_nu_r1_i3 * e0_r_r1 + e0_greek_nu_r2_i3 * e0_r_r2 + e0_greek_nu_r3_i3 * e0_r_r3)); 
	Y(33) = 0.0 - ((e0_F_j1 * e0_x_i4_j1 + e0_F_j2 * e0_x_i4_j2) + (e0_greek_nu_r1_i4 * e0_r_r1 + e0_greek_nu_r2_i4 * e0_r_r2 + e0_greek_nu_r3_i4 * e0_r_r3)); 
	Y(34) = 0.0 - ((e0_F_j1 * e0_x_i5_j1 + e0_F_j2 * e0_x_i5_j2) + (e0_greek_nu_r1_i5 * e0_r_r1 + e0_greek_nu_r2_i5 * e0_r_r2 + e0_greek_nu_r3_i5 * e0_r_r3)); 
	Y(35) = 0.0 - ((e0_F_j1 * e0_x_i6_j1 + e0_F_j2 * e0_x_i6_j2) + (e0_greek_nu_r1_i6 * e0_r_r1 + e0_greek_nu_r2_i6 * e0_r_r2 + e0_greek_nu_r3_i6 * e0_r_r3)); 
	Y(36) = 0.0 - ((e0_F_j1 * e0_x_i7_j1 + e0_F_j2 * e0_x_i7_j2) + (e0_greek_nu_r1_i7 * e0_r_r1 + e0_greek_nu_r2_i7 * e0_r_r2 + e0_greek_nu_r3_i7 * e0_r_r3)); 
	Y(37) = e0_X_i5 - ((e0_F_j1 * e0_x_i5_j1 + e0_F_j2 * e0_x_i5_j2)/(e0_F_j1 * e0_x_i5_j1)); 
	Y(38) = e0_X_i2 - ((e0_F_j1 * e0_x_i2_j1 + e0_F_j2 * e0_x_i2_j2)/(e0_F_j1 * e0_x_i2_j1)); 
	Y(39) = e0_Y_i7 - ( - 1.0 * (e0_F_j1 * e0_x_i7_j1 + e0_F_j2 * e0_x_i7_j2)/(e0_F_j1 * e0_x_i5_j1)); 
	Y(40) = e0_h_f_i1 + e0_A_par1_i1 * (power((e0_T_j1),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i1 * (power((e0_T_j1),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i1 * (power((e0_T_j1),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i1 * (power((e0_T_j1),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i1 * (power((e0_T_j1),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j1_i1); 
	Y(41) = e0_h_f_i2 + e0_A_par1_i2 * (power((e0_T_j1),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i2 * (power((e0_T_j1),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i2 * (power((e0_T_j1),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i2 * (power((e0_T_j1),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i2 * (power((e0_T_j1),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j1_i2); 
	Y(42) = e0_h_f_i3 + e0_A_par1_i3 * (power((e0_T_j1),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i3 * (power((e0_T_j1),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i3 * (power((e0_T_j1),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i3 * (power((e0_T_j1),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i3 * (power((e0_T_j1),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j1_i3); 
	Y(43) = e0_h_f_i4 + e0_A_par1_i4 * (power((e0_T_j1),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i4 * (power((e0_T_j1),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i4 * (power((e0_T_j1),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i4 * (power((e0_T_j1),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i4 * (power((e0_T_j1),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j1_i4); 
	Y(44) = e0_h_f_i5 + e0_A_par1_i5 * (power((e0_T_j1),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i5 * (power((e0_T_j1),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i5 * (power((e0_T_j1),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i5 * (power((e0_T_j1),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i5 * (power((e0_T_j1),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j1_i5); 
	Y(45) = e0_h_f_i6 + e0_A_par1_i6 * (power((e0_T_j1),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i6 * (power((e0_T_j1),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i6 * (power((e0_T_j1),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i6 * (power((e0_T_j1),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i6 * (power((e0_T_j1),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j1_i6); 
	Y(46) = e0_h_f_i7 + e0_A_par1_i7 * (power((e0_T_j1),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i7 * (power((e0_T_j1),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i7 * (power((e0_T_j1),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i7 * (power((e0_T_j1),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i7 * (power((e0_T_j1),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j1_i7); 
	Y(47) = e0_h_f_i1 + e0_A_par1_i1 * (power((e0_T_j2),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i1 * (power((e0_T_j2),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i1 * (power((e0_T_j2),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i1 * (power((e0_T_j2),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i1 * (power((e0_T_j2),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j2_i1); 
	Y(48) = e0_h_f_i2 + e0_A_par1_i2 * (power((e0_T_j2),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i2 * (power((e0_T_j2),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i2 * (power((e0_T_j2),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i2 * (power((e0_T_j2),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i2 * (power((e0_T_j2),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j2_i2); 
	Y(49) = e0_h_f_i3 + e0_A_par1_i3 * (power((e0_T_j2),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i3 * (power((e0_T_j2),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i3 * (power((e0_T_j2),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i3 * (power((e0_T_j2),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i3 * (power((e0_T_j2),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j2_i3); 
	Y(50) = e0_h_f_i4 + e0_A_par1_i4 * (power((e0_T_j2),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i4 * (power((e0_T_j2),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i4 * (power((e0_T_j2),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i4 * (power((e0_T_j2),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i4 * (power((e0_T_j2),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j2_i4); 
	Y(51) = e0_h_f_i5 + e0_A_par1_i5 * (power((e0_T_j2),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i5 * (power((e0_T_j2),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i5 * (power((e0_T_j2),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i5 * (power((e0_T_j2),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i5 * (power((e0_T_j2),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j2_i5); 
	Y(52) = e0_h_f_i6 + e0_A_par1_i6 * (power((e0_T_j2),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i6 * (power((e0_T_j2),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i6 * (power((e0_T_j2),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i6 * (power((e0_T_j2),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i6 * (power((e0_T_j2),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j2_i6); 
	Y(53) = e0_h_f_i7 + e0_A_par1_i7 * (power((e0_T_j2),5.0) - power((e0_T_f),5.0))/(5.0) + e0_A_par2_i7 * (power((e0_T_j2),4.0) - power((e0_T_f),4.0))/(4.0) + e0_A_par3_i7 * (power((e0_T_j2),3.0) - power((e0_T_f),3.0))/(3.0) + e0_A_par4_i7 * (power((e0_T_j2),2.0) - power((e0_T_f),2.0))/(2.0) + e0_A_par5_i7 * (power((e0_T_j2),1.0) - power((e0_T_f),1.0))/(1.0) - (e0_h_o_j2_i7); 
	Y(54) = e0_A_par1_r1 * exp((e0_E_par1_r1)/(e0_R * e0_T_j2)) - (e0_k_par1); 
	Y(55) = e0_A_par2_r1 * exp((e0_E_par2_r1)/(e0_R * e0_T_j2)) - (e0_k_par2); 
	Y(56) = e0_A_par3_r1 * exp((e0_E_par3_r1)/(e0_R * e0_T_j2)) - (e0_k_par3); 
	Y(57) = e0_A_par4_r1 * exp((e0_E_par4_r1)/(e0_R * e0_T_j2)) - (e0_k_par4); 
	Y(58) = e0_A_par5_r1 * exp((e0_E_par5_r1)/(e0_R * e0_T_j2)) - (e0_k_par5); 
	Y(59) = e0_A_r3_i3 * exp((e0_E_r3_i3)/(e0_R * e0_T_j2)) - (e0_K_i3); 
	Y(60) = e0_A_r3_i6 * exp((e0_E_r3_i6)/(e0_R * e0_T_j2)) - (e0_K_i6); 
	Y(61) = e0_A_r3 * exp((e0_E_r3)/(e0_R * e0_T_j2)) - (e0_k_r3); 



end

function[] = displayResults(X_ITER)

	% print variable values to display 
	disp(['e0_x_i2_j2 = ', num2str(X_ITER(1))]);
	disp(['e0_x_i1_j2 = ', num2str(X_ITER(2))]);
	disp(['e0_HU_i3 = ', num2str(X_ITER(3))]);
	disp(['e0_HU_i2 = ', num2str(X_ITER(4))]);
	disp(['e0_K_r3 = ', num2str(X_ITER(5))]);
	disp(['e0_HU_i1 = ', num2str(X_ITER(6))]);
	disp(['e0_u_j2 = ', num2str(X_ITER(7))]);
	disp(['e0_U = ', num2str(X_ITER(8))]);
	disp(['e0_HU = ', num2str(X_ITER(9))]);
	disp(['e0_p_i7_j2 = ', num2str(X_ITER(10))]);
	disp(['e0_p_i6_j2 = ', num2str(X_ITER(11))]);
	disp(['e0_p_i5_j2 = ', num2str(X_ITER(12))]);
	disp(['e0_p_i4_j2 = ', num2str(X_ITER(13))]);
	disp(['e0_p_i3_j2 = ', num2str(X_ITER(14))]);
	disp(['e0_p_i2_j2 = ', num2str(X_ITER(15))]);
	disp(['e0_p_i1_j2 = ', num2str(X_ITER(16))]);
	disp(['e0_h_o_j2_i7 = ', num2str(X_ITER(17))]);
	disp(['e0_h_o_j1_i7 = ', num2str(X_ITER(18))]);
	disp(['e0_h_o_j2_i6 = ', num2str(X_ITER(19))]);
	disp(['e0_h_o_j1_i6 = ', num2str(X_ITER(20))]);
	disp(['e0_Q = ', num2str(X_ITER(21))]);
	disp(['e0_h_o_j2_i5 = ', num2str(X_ITER(22))]);
	disp(['e0_h_o_j1_i5 = ', num2str(X_ITER(23))]);
	disp(['e0_h_o_j2_i4 = ', num2str(X_ITER(24))]);
	disp(['e0_h_o_j1_i4 = ', num2str(X_ITER(25))]);
	disp(['e0_h_o_j2_i3 = ', num2str(X_ITER(26))]);
	disp(['e0_F_j2 = ', num2str(X_ITER(27))]);
	disp(['e0_h_o_j1_i3 = ', num2str(X_ITER(28))]);
	disp(['e0_h_o_j2_i2 = ', num2str(X_ITER(29))]);
	disp(['e0_r_r1 = ', num2str(X_ITER(30))]);
	disp(['e0_h_o_j1_i2 = ', num2str(X_ITER(31))]);
	disp(['e0_k_par5 = ', num2str(X_ITER(32))]);
	disp(['e0_h_o_j2_i1 = ', num2str(X_ITER(33))]);
	disp(['e0_r_r2 = ', num2str(X_ITER(34))]);
	disp(['e0_h_o_j1_i1 = ', num2str(X_ITER(35))]);
	disp(['e0_k_par4 = ', num2str(X_ITER(36))]);
	disp(['e0_k_par3 = ', num2str(X_ITER(37))]);
	disp(['e0_k_par2 = ', num2str(X_ITER(38))]);
	disp(['e0_k_par1 = ', num2str(X_ITER(39))]);
	disp(['e0_h_j2 = ', num2str(X_ITER(40))]);
	disp(['e0_K_r1 = ', num2str(X_ITER(41))]);
	disp(['e0_h_j1 = ', num2str(X_ITER(42))]);
	disp(['e0_x_i7_j2 = ', num2str(X_ITER(43))]);
	disp(['e0_x_i7_j1 = ', num2str(X_ITER(44))]);
	disp(['e0_x_i6_j2 = ', num2str(X_ITER(45))]);
	disp(['e0_r_r3 = ', num2str(X_ITER(46))]);
	disp(['e0_K_r2 = ', num2str(X_ITER(47))]);
	disp(['e0_k_r3 = ', num2str(X_ITER(48))]);
	disp(['e0_x_i5_j2 = ', num2str(X_ITER(49))]);
	disp(['e0_K_i6 = ', num2str(X_ITER(50))]);
	disp(['e0_K_i3 = ', num2str(X_ITER(51))]);
	disp(['e0_x_i4_j2 = ', num2str(X_ITER(52))]);
	disp(['e0_v_j2 = ', num2str(X_ITER(53))]);
	disp(['e0_HU_i7 = ', num2str(X_ITER(54))]);
	disp(['e0_HU_i6 = ', num2str(X_ITER(55))]);
	disp(['e0_HU_i5 = ', num2str(X_ITER(56))]);
	disp(['e0_HU_i4 = ', num2str(X_ITER(57))]);
	disp(['e0_X_i2 = ', num2str(X_ITER(58))]);
	disp(['e0_Y_i7 = ', num2str(X_ITER(59))]);
	disp(['e0_x_i3_j2 = ', num2str(X_ITER(60))]);
	disp(['e0_X_i5 = ', num2str(X_ITER(61))]);


end


