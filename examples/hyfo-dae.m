%*********************************************************
% The namespaces have been normalized. The following
% table shows the attribuation. 
% Normalized Name --> Full Name ---> User-defined Name
% =================================== 
% e0 --> e[0]87906 --> 
%*********************************************************

%*********************************************************
% The variables are named according to the notation
% provided in the Mosaic model.
% 
% The variable names can be read as follows:
% ==========================================
% 	e0_V_Reactor
% 		V: Volume [l]
% 		Superscripts
% 			Reactor: Reactor related Variables
% 	 
% 	e0_E_r#
% 		E: Activation Energy [J/mol]
% 		Indices
% 			r: Number of reactions: 1 = Iso; 2 = HydA; 3 = HydB; 4 = HyfoA; 5 = HyfoB; 6 = HyfoC;
% 	 
% 	e0_k_ref_r#
% 		k: Preexponential factor
% 		Superscripts
% 			ref: Enhancement factor with reference temperature
% 		Indices
% 			r: Number of reactions: 1 = Iso; 2 = HydA; 3 = HydB; 4 = HyfoA; 5 = HyfoB; 6 = HyfoC;
% 	 
% 	e0_n_Surfactant
% 		n: Amount of moles [mol]
% 		Superscripts
% 			Surfactant: Surfactant
% 	 
% 	e0_n_Water
% 		n: Amount of moles [mol]
% 		Superscripts
% 			Water: Water as a pure Component
% 	 
% 	e0_M_i#
% 		M: Molar mass [g/mol]
% 		Indices
% 			i: Number of components: 1 = Doce; 2 = iso_Doce; 3 = iso-TDC; 4 = Doca; 5 = TDC; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant;
% 	 
% 	e0_M_Water
% 		M: Molar mass [g/mol]
% 		Superscripts
% 			Water: Water as a pure Component
% 	 
% 	e0_M_Surfactant
% 		M: Molar mass [g/mol]
% 		Superscripts
% 			Surfactant: Surfactant
% 	 
% 	e0_T
% 		T: Temperature [K]
% 	 
% 	e0_p_Reactor
% 		p: Pressure [bar]
% 		Superscripts
% 			Reactor: Reactor related Variables
% 	 
% 	e0_K_cat_e#
% 		K: Inhibition factor
% 		Subscript
% 			cat: catalyst
% 		Indices
% 			e: e = 1 = alpha; e = 2 = beta; e = 3 = delta
% 	 
% 	e0_c_cat
% 		c: Concentration [mol/l] or DIPPR 105 Parameter
% 		Subscript
% 			cat: catalyst
% 	 
% 	e0_R
% 		R: Ideal Gas Constant [J/(mol K)]
% 	 
% 	e0_greek_DeltaG_r#
% 		&Delta;G: Gibbs free enthalpy
% 		Indices
% 			r: Number of reactions: 1 = Iso; 2 = HydA; 3 = HydB; 4 = HyfoA; 5 = HyfoB; 6 = HyfoC;
% 	 
% 	e0_K_r#_e#
% 		K: Inhibition factor
% 		Indices
% 			r: Number of reactions: 1 = Iso; 2 = HydA; 3 = HydB; 4 = HyfoA; 5 = HyfoB; 6 = HyfoC;
% 			e: e = 1 = alpha; e = 2 = beta; e = 3 = delta
% 	 
% 	e0_K_LM
% 		K: Inhibition factor
% 		Superscripts
% 			LM: Ligand to Metal Ratio
% 	 
% 	e0_P_Surfactant
% 		P: Model Parameter
% 		Superscripts
% 			Surfactant: Surfactant
% 	 
% 	e0_P_trig_r#
% 		P: Model Parameter
% 		Superscripts
% 			trig: Trigger Tuning
% 		Indices
% 			r: Number of reactions: 1 = Iso; 2 = HydA; 3 = HydB; 4 = HyfoA; 5 = HyfoB; 6 = HyfoC;
% 	 
% 	e0_T_ref
% 		T: Temperature [K]
% 		Superscripts
% 			ref: Enhancement factor with reference temperature
% 	 
% 	e0_k_LM_r#
% 		k: Preexponential factor
% 		Superscripts
% 			LM: Ligand to Metal Ratio
% 		Indices
% 			r: Number of reactions: 1 = Iso; 2 = HydA; 3 = HydB; 4 = HyfoA; 5 = HyfoB; 6 = HyfoC;
% 	 
% 	e0_n_Cat
% 		n: Amount of moles [mol]
% 		Superscripts
% 			Cat: Catalyst
% 	 
% 	e0_n_Lig
% 		n: Amount of moles [mol]
% 		Superscripts
% 			Lig: Ligand
% 	 
% 	e0_P_trig_Hyfo
% 		P: Model Parameter
% 		Superscripts
% 			trig: Trigger Tuning
% 		Subscript
% 			Hyfo: Hydroformylation
% 	 
% 	e0_k_LM_Hyfo
% 		k: Preexponential factor
% 		Superscripts
% 			LM: Ligand to Metal Ratio
% 		Subscript
% 			Hyfo: Hydroformylation
% 	 
% 	e0_c_i#
% 		c: Concentration [mol/l] or DIPPR 105 Parameter
% 		Indices
% 			i: Number of components: 1 = Doce; 2 = iso_Doce; 3 = iso-TDC; 4 = Doca; 5 = TDC; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant;
% 	 
% 	e0_n_i#
% 		n: Amount of moles [mol]
% 		Indices
% 			i: Number of components: 1 = Doce; 2 = iso_Doce; 3 = iso-TDC; 4 = Doca; 5 = TDC; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant;
% 	 
% 	e0_r_r#
% 		r: Reaction rate [g/(l h)]
% 		Indices
% 			r: Number of reactions: 1 = Iso; 2 = HydA; 3 = HydB; 4 = HyfoA; 5 = HyfoB; 6 = HyfoC;
% 	 
% 	e0_r_i#
% 		r: Reaction rate [g/(l h)]
% 		Indices
% 			i: Number of components: 1 = Doce; 2 = iso_Doce; 3 = iso-TDC; 4 = Doca; 5 = TDC; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant;
% 	 
% 	e0_n_L
% 		n: Amount of moles [mol]
% 		Superscripts
% 			L: Liquid phase
% 	 
% 	e0_greek_alpha
% 		&alpha;: Oil to water ratio
% 	 
% 	e0_greek_gamma
% 		&gamma;: Surfactant ratio
% 	 
% 	e0_X
% 		X: Conversion
% 	 
% 	e0_x_i#
% 		x: Mol fraction [mol/mol]
% 		Indices
% 			i: Number of components: 1 = Doce; 2 = iso_Doce; 3 = iso-TDC; 4 = Doca; 5 = TDC; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant;
% 	 
% 	e0_greek_psi_cat
% 		&psi;: Reformulation factor
% 		Subscript
% 			cat: catalyst
% 	 
% 	e0_K_eq_r#
% 		K: Inhibition factor
% 		Superscripts
% 			eq: Equilibrium
% 		Indices
% 			r: Number of reactions: 1 = Iso; 2 = HydA; 3 = HydB; 4 = HyfoA; 5 = HyfoB; 6 = HyfoC;
% 	 
% 	e0_P_i#_Sol#
% 		P: Model Parameter
% 		Indices
% 			i: Number of components: 1 = Doce; 2 = iso_Doce; 3 = iso-TDC; 4 = Doca; 5 = TDC; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant;
% 			Sol: Number of Solubility Model Parameter
% 	 
%*********************************************************

function[X,Y]=solveEquationSystem()

	% declare interval of independent variable
	X_START=0.0;
	X_END=4.0;
	X_INTERVAL=[X_START X_END];	% e0_t

	% specify the initial values for the state variables 
	Y_INIT(1) = 2.2564697;  	% e0_c_i1  
	Y_INIT(2) = 0.08573363;  	% e0_c_i2  
	Y_INIT(3) = 1.0E-7;  	% e0_c_i3  
	Y_INIT(4) = 0.00824684;  	% e0_c_i4  
	Y_INIT(5) = 1.0E-7;  	% e0_c_i5  
	Y_INIT(6) = 1.967716;  	% e0_n_i1  
	Y_INIT(7) = 0.07476256;  	% e0_n_i2  
	Y_INIT(8) = 8.72E-8;  	% e0_n_i3  
	Y_INIT(9) = 0.007191517;  	% e0_n_i4  
	Y_INIT(10) = 8.72E-8;  	% e0_n_i5  
	Y_INIT(11) = 21.322771;  	% e0_n_L  
	Y_INIT(12) = 0.49993634;  	% e0_greek_alpha  
	Y_INIT(13) = 0.08001318;  	% e0_greek_gamma  
	Y_INIT(14) = 1.0E-7;  	% e0_X  
	Y_INIT(15) = 9.62E-4;  	% e0_x_i6  
	Y_INIT(16) = 9.57E-4;  	% e0_x_i7  
	Y_INIT(17) = 0.023546033;  	% e0_c_i6  
	Y_INIT(18) = 0.023419946;  	% e0_c_i7  
	Y_INIT(19) = 2.06E-4;  	% e0_greek_psi_cat  
	Y_INIT(20) = -77498.87;  	% e0_greek_DeltaG_r3  
	Y_INIT(21) = 9.91E10;  	% e0_K_eq_r3  
	Y_INIT(22) = 260178.73;  	% e0_K_eq_r1  
	Y_INIT(23) = 1.13E-4;  	% e0_r_r1  
	Y_INIT(24) = 3.04E-10;  	% e0_r_r2  
	Y_INIT(25) = 1.08E-5;  	% e0_r_r3  
	Y_INIT(26) = 3.05E-5;  	% e0_r_r4  
	Y_INIT(27) = 0.003999636;  	% e0_r_r5  
	Y_INIT(28) = 4.11E-10;  	% e0_r_r6  
	Y_INIT(29) = -0.004123436;  	% e0_r_i1  
	Y_INIT(30) = 8.25E-5;  	% e0_r_i2  
	Y_INIT(31) = 3.05E-5;  	% e0_r_i3  
	Y_INIT(32) = 1.08E-5;  	% e0_r_i4  
	Y_INIT(33) = 0.003999636;  	% e0_r_i5  

	% declare parameters 
	PARAMS(1) = 0.872033;  	% e0_V_Reactor 
	PARAMS(2) = 0.117;  	% e0_n_Surfactant 
	PARAMS(3) = 19.1561;  	% e0_n_Water 
	PARAMS(4) = 168.32;  	% e0_M_i1 
	PARAMS(5) = 168.32;  	% e0_M_i2 
	PARAMS(6) = 198.34;  	% e0_M_i3 
	PARAMS(7) = 170.34;  	% e0_M_i4 
	PARAMS(8) = 198.34;  	% e0_M_i5 
	PARAMS(9) = 18.0153;  	% e0_M_Water 
	PARAMS(10) = 513.0;  	% e0_M_Surfactant 
	PARAMS(11) = 368.15;  	% e0_T 
	PARAMS(12) = 15.0;  	% e0_p_Reactor 
	PARAMS(13) = -6.4909E-5;  	% e0_P_i6_Sol1 
	PARAMS(14) = 1.1885E-5;  	% e0_P_i6_Sol2 
	PARAMS(15) = 0.0010631;  	% e0_P_i6_Sol3 
	PARAMS(16) = -0.027378;  	% e0_P_i6_Sol4 
	PARAMS(17) = 1.7599E-4;  	% e0_P_i6_Sol5 
	PARAMS(18) = 0.17476;  	% e0_P_i6_Sol6 
	PARAMS(19) = 9.2954E-4;  	% e0_P_i6_Sol7 
	PARAMS(20) = 2.8881E-7;  	% e0_P_i6_Sol8 
	PARAMS(21) = 2.9467E-4;  	% e0_P_i6_Sol9 
	PARAMS(22) = 3.7274E-4;  	% e0_P_i6_Sol10 
	PARAMS(23) = -4.1033E-5;  	% e0_P_i6_Sol11 
	PARAMS(24) = -9.9645E-6;  	% e0_P_i6_Sol12 
	PARAMS(25) = -3.8368E-5;  	% e0_P_i6_Sol13 
	PARAMS(26) = -6.9782E-6;  	% e0_P_i6_Sol14 
	PARAMS(27) = -8.2558E-5;  	% e0_P_i6_Sol15 
	PARAMS(28) = -1.7718E-4;  	% e0_P_i7_Sol1 
	PARAMS(29) = 1.7692E-5;  	% e0_P_i7_Sol2 
	PARAMS(30) = 0.0016934;  	% e0_P_i7_Sol3 
	PARAMS(31) = -0.047302;  	% e0_P_i7_Sol4 
	PARAMS(32) = 4.3746E-4;  	% e0_P_i7_Sol5 
	PARAMS(33) = 0.28638;  	% e0_P_i7_Sol6 
	PARAMS(34) = 0.001592;  	% e0_P_i7_Sol7 
	PARAMS(35) = -1.7107E-7;  	% e0_P_i7_Sol8 
	PARAMS(36) = 6.5328E-4;  	% e0_P_i7_Sol9 
	PARAMS(37) = 5.3043E-4;  	% e0_P_i7_Sol10 
	PARAMS(38) = -7.299E-6;  	% e0_P_i7_Sol11 
	PARAMS(39) = -1.4868E-5;  	% e0_P_i7_Sol12 
	PARAMS(40) = -3.0261E-5;  	% e0_P_i7_Sol13 
	PARAMS(41) = -1.2455E-5;  	% e0_P_i7_Sol14 
	PARAMS(42) = -1.1598E-4;  	% e0_P_i7_Sol15 
	PARAMS(43) = 45087.07;  	% e0_K_cat_e1 
	PARAMS(44) = 189.31375;  	% e0_K_cat_e2 
	PARAMS(45) = 0.25682598;  	% e0_c_cat 
	PARAMS(46) = 8.314;  	% e0_R 
	PARAMS(47) = 38165.484;  	% e0_greek_DeltaG_r1 
	PARAMS(48) = 40749.277;  	% e0_E_r1 
	PARAMS(49) = 0.72770315;  	% e0_K_r1_e1 
	PARAMS(50) = 4.05E-5;  	% e0_K_r1_e2 
	PARAMS(51) = 2.7251527;  	% e0_K_LM 
	PARAMS(52) = 1.0315819;  	% e0_P_Surfactant 
	PARAMS(53) = 14.282191;  	% e0_P_trig_r1 
	PARAMS(54) = 363.15;  	% e0_T_ref 
	PARAMS(55) = 66.92345;  	% e0_k_LM_r1 
	PARAMS(56) = 4.242135;  	% e0_k_ref_r1 
	PARAMS(57) = 8.58E-4;  	% e0_n_Cat 
	PARAMS(58) = 0.0043;  	% e0_n_Lig 
	PARAMS(59) = 6285.8706;  	% e0_E_r2 
	PARAMS(60) = 0.005641299;  	% e0_k_ref_r2 
	PARAMS(61) = 104496.37;  	% e0_E_r3 
	PARAMS(62) = 0.47820565;  	% e0_K_r3_e1 
	PARAMS(63) = 13262.677;  	% e0_K_r3_e2 
	PARAMS(64) = 1028.9795;  	% e0_K_r3_e3 
	PARAMS(65) = 17428.53;  	% e0_k_ref_r3 
	PARAMS(66) = 107045.41;  	% e0_E_r4 
	PARAMS(67) = 11.312137;  	% e0_P_trig_Hyfo 
	PARAMS(68) = 15349.087;  	% e0_k_ref_r4 
	PARAMS(69) = 1.0487578;  	% e0_k_LM_Hyfo 
	PARAMS(70) = 57858.113;  	% e0_E_r5 
	PARAMS(71) = 0.023340752;  	% e0_K_r5_e1 
	PARAMS(72) = 895.06036;  	% e0_K_r5_e2 
	PARAMS(73) = 44226.242;  	% e0_K_r5_e3 
	PARAMS(74) = 9.94E7;  	% e0_k_ref_r5 
	PARAMS(75) = 32422.021;  	% e0_E_r6 
	PARAMS(76) = 0.010987442;  	% e0_k_ref_r6 

	M = daeSystemMM();

	OPTIONS = odeset('Mass',M);

	[X,Y] = ode15s(@(X,Y)daeSystemLHS(X,Y,PARAMS),X_INTERVAL,Y_INIT',OPTIONS);

	displayResults(X,Y);

end


function MASS = daeSystemMM()

	MASS = zeros(33, 33); %total number of equations
	MASS(1:5, 1:5)=eye(5,5); % number of odes

end


% evaluate the differential function.
function[DYDX] = daeSystemLHS(X,Y,PARAMS)

	% read out variables  
	e0_c_i1 = Y(1); 
	e0_c_i2 = Y(2); 
	e0_c_i3 = Y(3); 
	e0_c_i4 = Y(4); 
	e0_c_i5 = Y(5); 
	e0_n_i1 = Y(6); 
	e0_n_i2 = Y(7); 
	e0_n_i3 = Y(8); 
	e0_n_i4 = Y(9); 
	e0_n_i5 = Y(10); 
	e0_n_L = Y(11); 
	e0_greek_alpha = Y(12); 
	e0_greek_gamma = Y(13); 
	e0_X = Y(14); 
	e0_x_i6 = Y(15); 
	e0_x_i7 = Y(16); 
	e0_c_i6 = Y(17); 
	e0_c_i7 = Y(18); 
	e0_greek_psi_cat = Y(19); 
	e0_greek_DeltaG_r3 = Y(20); 
	e0_K_eq_r3 = Y(21); 
	e0_K_eq_r1 = Y(22); 
	e0_r_r1 = Y(23); 
	e0_r_r2 = Y(24); 
	e0_r_r3 = Y(25); 
	e0_r_r4 = Y(26); 
	e0_r_r5 = Y(27); 
	e0_r_r6 = Y(28); 
	e0_r_i1 = Y(29); 
	e0_r_i2 = Y(30); 
	e0_r_i3 = Y(31); 
	e0_r_i4 = Y(32); 
	e0_r_i5 = Y(33); 

	% read out differential variable
	e0_t = X;

	% read out parameters 
	e0_V_Reactor = PARAMS(1); 
	e0_n_Surfactant = PARAMS(2); 
	e0_n_Water = PARAMS(3); 
	e0_M_i1 = PARAMS(4); 
	e0_M_i2 = PARAMS(5); 
	e0_M_i3 = PARAMS(6); 
	e0_M_i4 = PARAMS(7); 
	e0_M_i5 = PARAMS(8); 
	e0_M_Water = PARAMS(9); 
	e0_M_Surfactant = PARAMS(10); 
	e0_T = PARAMS(11); 
	e0_p_Reactor = PARAMS(12); 
	e0_P_i6_Sol1 = PARAMS(13); 
	e0_P_i6_Sol2 = PARAMS(14); 
	e0_P_i6_Sol3 = PARAMS(15); 
	e0_P_i6_Sol4 = PARAMS(16); 
	e0_P_i6_Sol5 = PARAMS(17); 
	e0_P_i6_Sol6 = PARAMS(18); 
	e0_P_i6_Sol7 = PARAMS(19); 
	e0_P_i6_Sol8 = PARAMS(20); 
	e0_P_i6_Sol9 = PARAMS(21); 
	e0_P_i6_Sol10 = PARAMS(22); 
	e0_P_i6_Sol11 = PARAMS(23); 
	e0_P_i6_Sol12 = PARAMS(24); 
	e0_P_i6_Sol13 = PARAMS(25); 
	e0_P_i6_Sol14 = PARAMS(26); 
	e0_P_i6_Sol15 = PARAMS(27); 
	e0_P_i7_Sol1 = PARAMS(28); 
	e0_P_i7_Sol2 = PARAMS(29); 
	e0_P_i7_Sol3 = PARAMS(30); 
	e0_P_i7_Sol4 = PARAMS(31); 
	e0_P_i7_Sol5 = PARAMS(32); 
	e0_P_i7_Sol6 = PARAMS(33); 
	e0_P_i7_Sol7 = PARAMS(34); 
	e0_P_i7_Sol8 = PARAMS(35); 
	e0_P_i7_Sol9 = PARAMS(36); 
	e0_P_i7_Sol10 = PARAMS(37); 
	e0_P_i7_Sol11 = PARAMS(38); 
	e0_P_i7_Sol12 = PARAMS(39); 
	e0_P_i7_Sol13 = PARAMS(40); 
	e0_P_i7_Sol14 = PARAMS(41); 
	e0_P_i7_Sol15 = PARAMS(42); 
	e0_K_cat_e1 = PARAMS(43); 
	e0_K_cat_e2 = PARAMS(44); 
	e0_c_cat = PARAMS(45); 
	e0_R = PARAMS(46); 
	e0_greek_DeltaG_r1 = PARAMS(47); 
	e0_E_r1 = PARAMS(48); 
	e0_K_r1_e1 = PARAMS(49); 
	e0_K_r1_e2 = PARAMS(50); 
	e0_K_LM = PARAMS(51); 
	e0_P_Surfactant = PARAMS(52); 
	e0_P_trig_r1 = PARAMS(53); 
	e0_T_ref = PARAMS(54); 
	e0_k_LM_r1 = PARAMS(55); 
	e0_k_ref_r1 = PARAMS(56); 
	e0_n_Cat = PARAMS(57); 
	e0_n_Lig = PARAMS(58); 
	e0_E_r2 = PARAMS(59); 
	e0_k_ref_r2 = PARAMS(60); 
	e0_E_r3 = PARAMS(61); 
	e0_K_r3_e1 = PARAMS(62); 
	e0_K_r3_e2 = PARAMS(63); 
	e0_K_r3_e3 = PARAMS(64); 
	e0_k_ref_r3 = PARAMS(65); 
	e0_E_r4 = PARAMS(66); 
	e0_P_trig_Hyfo = PARAMS(67); 
	e0_k_ref_r4 = PARAMS(68); 
	e0_k_LM_Hyfo = PARAMS(69); 
	e0_E_r5 = PARAMS(70); 
	e0_K_r5_e1 = PARAMS(71); 
	e0_K_r5_e2 = PARAMS(72); 
	e0_K_r5_e3 = PARAMS(73); 
	e0_k_ref_r5 = PARAMS(74); 
	e0_E_r6 = PARAMS(75); 
	e0_k_ref_r6 = PARAMS(76); 

	% evaluate the function values  
	DYDX(1) = e0_r_i1 * 60.0; 
	DYDX(2) = e0_r_i2 * 60.0; 
	DYDX(3) = e0_r_i3 * 60.0; 
	DYDX(4) = e0_r_i4 * 60.0; 
	DYDX(5) = e0_r_i5 * 60.0; 
	DYDX(6) = e0_n_i1 - ( e0_c_i1 * e0_V_Reactor ); 
	DYDX(7) = e0_n_i2 - ( e0_c_i2 * e0_V_Reactor ); 
	DYDX(8) = e0_n_i3 - ( e0_c_i3 * e0_V_Reactor ); 
	DYDX(9) = e0_n_i4 - ( e0_c_i4 * e0_V_Reactor ); 
	DYDX(10) = e0_n_i5 - ( e0_c_i5 * e0_V_Reactor ); 
	DYDX(11) = e0_n_L - ( ( e0_n_i1 + e0_n_i2 + e0_n_i3 + e0_n_i4 + e0_n_i5 ) + e0_n_Water + e0_n_Surfactant ); 
	DYDX(12) = e0_greek_alpha - ( ( ( e0_c_i1 * e0_V_Reactor * e0_M_i1 + e0_c_i2 * e0_V_Reactor * e0_M_i2 + e0_c_i3 * e0_V_Reactor * e0_M_i3 + e0_c_i4 * e0_V_Reactor * e0_M_i4 + e0_c_i5 * e0_V_Reactor * e0_M_i5 ) )/( ( e0_c_i1 * e0_V_Reactor * e0_M_i1 + e0_c_i2 * e0_V_Reactor * e0_M_i2 + e0_c_i3 * e0_V_Reactor * e0_M_i3 + e0_c_i4 * e0_V_Reactor * e0_M_i4 + e0_c_i5 * e0_V_Reactor * e0_M_i5 ) + e0_n_Water * e0_M_Water ) ); 
	DYDX(13) = e0_greek_gamma - ( ( e0_n_Surfactant * e0_M_Surfactant )/( ( e0_c_i1 * e0_V_Reactor * e0_M_i1 + e0_c_i2 * e0_V_Reactor * e0_M_i2 + e0_c_i3 * e0_V_Reactor * e0_M_i3 + e0_c_i4 * e0_V_Reactor * e0_M_i4 + e0_c_i5 * e0_V_Reactor * e0_M_i5 ) + e0_n_Water * e0_M_Water + e0_n_Surfactant * e0_M_Surfactant ) ); 
	DYDX(14) = e0_X - ( ( ( e0_c_i3 * e0_M_i3 + e0_c_i5 * e0_M_i5 ) * e0_V_Reactor )/( ( e0_c_i1 * e0_V_Reactor * e0_M_i1 + e0_c_i2 * e0_V_Reactor * e0_M_i2 + e0_c_i3 * e0_V_Reactor * e0_M_i3 + e0_c_i4 * e0_V_Reactor * e0_M_i4 + e0_c_i5 * e0_V_Reactor * e0_M_i5 ) ) ); 
	DYDX(15) = e0_x_i6 - ( ( e0_p_Reactor * e0_P_i6_Sol1 + ( e0_T - 273.15 ) * e0_P_i6_Sol2 + e0_greek_alpha * e0_P_i6_Sol3 + e0_greek_gamma * e0_P_i6_Sol4 + e0_X * e0_P_i6_Sol5 + ( ( e0_greek_gamma ) )^( 2.0 ) * e0_P_i6_Sol6 + ( ( e0_X ) )^( 2.0 ) * e0_P_i6_Sol7 + e0_p_Reactor * ( e0_T - 273.15 ) * e0_P_i6_Sol8 + e0_p_Reactor * e0_greek_alpha * e0_P_i6_Sol9 + e0_p_Reactor * e0_greek_gamma * e0_P_i6_Sol10 + e0_p_Reactor * e0_X * e0_P_i6_Sol11 + ( e0_T - 273.15 ) * e0_greek_alpha * e0_P_i6_Sol12 + ( e0_T - 273.15 ) * e0_greek_gamma * e0_P_i6_Sol13 + ( e0_T - 273.15 ) * e0_X * e0_P_i6_Sol14 + e0_greek_alpha * e0_X * e0_P_i6_Sol15 ) - e0_x_i7 ); 
	DYDX(16) = e0_x_i7 - ( ( e0_p_Reactor )/( 2.0 ) * e0_P_i7_Sol1 + ( e0_T - 273.15 ) * e0_P_i7_Sol2 + e0_greek_alpha * e0_P_i7_Sol3 + e0_greek_gamma * e0_P_i7_Sol4 + e0_X * e0_P_i7_Sol5 + ( ( e0_greek_gamma ) )^( 2.0 ) * e0_P_i7_Sol6 + ( ( e0_X ) )^( 2.0 ) * e0_P_i7_Sol7 + ( e0_p_Reactor )/( 2.0 ) * ( e0_T - 273.15 ) * e0_P_i7_Sol8 + ( e0_p_Reactor )/( 2.0 ) * e0_greek_alpha * e0_P_i7_Sol9 + ( e0_p_Reactor )/( 2.0 ) * e0_greek_gamma * e0_P_i7_Sol10 + ( e0_p_Reactor )/( 2.0 ) * e0_X * e0_P_i7_Sol11 + ( e0_T - 273.15 ) * e0_greek_alpha * e0_P_i7_Sol12 + ( e0_T - 273.15 ) * e0_greek_gamma * e0_P_i7_Sol13 + ( e0_T - 273.15 ) * e0_X * e0_P_i7_Sol14 + e0_greek_alpha * e0_X * e0_P_i7_Sol15 ); 
	DYDX(17) = e0_c_i6 * e0_V_Reactor - ( ( e0_n_L * e0_x_i6 )/( 1.0 - e0_x_i6 ) ); 
	DYDX(18) = e0_c_i7 * e0_V_Reactor - ( ( e0_n_L * e0_x_i7 )/( 1.0 - e0_x_i7 ) ); 
	DYDX(19) = e0_greek_psi_cat * ( 1.0 + e0_K_cat_e1 * e0_c_i7 + e0_K_cat_e2 * ( e0_c_i7 )/( e0_c_i6 ) ) - ( e0_c_cat ); 
	DYDX(20) = e0_greek_DeltaG_r3 - ( (  - 126.28 + 0.13 * e0_T + 6.8 * ( ( 10.0 ) )^(  - 6.0 ) * ( ( e0_T ) )^( 2.0 ) ) * ( ( 10.0 ) )^( 3.0 ) ); 
	DYDX(21) = e0_K_eq_r3 - ( exp(  - ( e0_greek_DeltaG_r3 )/( e0_R * e0_T ) ) ); 
	DYDX(22) = e0_K_eq_r1 - ( exp( ( e0_greek_DeltaG_r1 )/( e0_R * e0_T ) ) ); 
	DYDX(23) = e0_r_r1 * ( 1.0 + e0_K_r1_e1 * e0_c_i1 + e0_K_r1_e2 * e0_c_i2 ) - ( ( ( ( e0_n_Surfactant )/( e0_V_Reactor ) ) )^( e0_P_Surfactant ) * ( 1.0 + ( e0_k_LM_r1 )/( 1.0 + exp(  - ( e0_K_LM - ( e0_n_Lig )/( e0_n_Cat ) ) * e0_P_trig_r1 ) ) ) * e0_greek_psi_cat * e0_k_ref_r1 * exp(  - ( e0_E_r1 )/( e0_R ) * ( ( 1.0 )/( e0_T ) - ( 1.0 )/( e0_T_ref ) ) ) * ( e0_c_i1 - ( e0_c_i2 )/( e0_K_eq_r1 ) ) ); 
	DYDX(24) = e0_r_r2 - ( ( ( ( e0_n_Surfactant )/( e0_V_Reactor ) ) )^( e0_P_Surfactant ) * e0_greek_psi_cat * e0_k_ref_r2 * exp(  - ( e0_E_r2 )/( e0_R ) * ( ( 1.0 )/( e0_T ) - ( 1.0 )/( e0_T_ref ) ) ) * e0_c_i2 * e0_c_i6 ); 
	DYDX(25) = e0_r_r3 * ( 1.0 + e0_K_r3_e1 * e0_c_i1 + e0_K_r3_e2 * e0_c_i4 + e0_K_r3_e3 * e0_c_i6 ) - ( ( ( ( e0_n_Surfactant )/( e0_V_Reactor ) ) )^( e0_P_Surfactant ) * e0_greek_psi_cat * e0_k_ref_r3 * exp(  - ( e0_E_r3 )/( e0_R ) * ( ( 1.0 )/( e0_T ) - ( 1.0 )/( e0_T_ref ) ) ) * ( e0_c_i2 * e0_c_i6 - ( e0_c_i4 )/( e0_K_eq_r3 ) ) ); 
	DYDX(26) = e0_r_r4 - ( ( ( ( e0_n_Surfactant )/( e0_V_Reactor ) ) )^( e0_P_Surfactant ) * ( 1.0 + ( e0_k_LM_Hyfo )/( 1.0 + exp(  - ( e0_K_LM - ( e0_n_Lig )/( e0_n_Cat ) ) * e0_P_trig_Hyfo ) ) ) * e0_greek_psi_cat * e0_k_ref_r4 * exp(  - ( e0_E_r4 )/( e0_R ) * ( ( 1.0 )/( e0_T ) - ( 1.0 )/( e0_T_ref ) ) ) * e0_c_i2 * e0_c_i6 * e0_c_i7 ); 
	DYDX(27) = e0_r_r5 * ( 1.0 + e0_K_r5_e1 * e0_c_i1 + e0_K_r5_e2 * e0_c_i5 + e0_K_r5_e3 * e0_c_i6 ) - ( ( ( ( e0_n_Surfactant )/( e0_V_Reactor ) ) )^( e0_P_Surfactant ) * ( 1.0 + ( e0_k_LM_Hyfo )/( 1.0 + exp(  - ( e0_K_LM - ( e0_n_Lig )/( e0_n_Cat ) ) * e0_P_trig_Hyfo ) ) ) * e0_greek_psi_cat * e0_k_ref_r5 * exp(  - ( e0_E_r5 )/( e0_R ) * ( ( 1.0 )/( e0_T ) - ( 1.0 )/( e0_T_ref ) ) ) * e0_c_i1 * e0_c_i6 * e0_c_i7 ); 
	DYDX(28) = e0_r_r6 - ( ( ( ( e0_n_Surfactant )/( e0_V_Reactor ) ) )^( e0_P_Surfactant ) * ( 1.0 + ( e0_k_LM_Hyfo )/( 1.0 + exp(  - ( e0_K_LM - ( e0_n_Lig )/( e0_n_Cat ) ) * e0_P_trig_Hyfo ) ) ) * e0_greek_psi_cat * e0_k_ref_r6 * exp(  - ( e0_E_r6 )/( e0_R ) * ( ( 1.0 )/( e0_T ) - ( 1.0 )/( e0_T_ref ) ) ) * e0_c_i1 * e0_c_i6 * e0_c_i7 ); 
	DYDX(29) = e0_r_i1 - (  - e0_r_r1 - e0_r_r3 - e0_r_r5 - e0_r_r6 ); 
	DYDX(30) = e0_r_i2 - ( e0_r_r1 - e0_r_r2 - e0_r_r4 ); 
	DYDX(31) = e0_r_i3 - ( e0_r_r4 + e0_r_r6 ); 
	DYDX(32) = e0_r_i4 - ( e0_r_r2 + e0_r_r3 ); 
	DYDX(33) = e0_r_i5 - ( e0_r_r5 ); 

	DYDX=DYDX';

end

function[] = displayResults(X,Y)

	% decide for a plot type: 
	%   0 	-> Plot the variables into individual figures 
	%   1 	-> Plot into sub figures 
	%   2 	-> Plot all selected into one figure 
	%   other 	-> Do not plot 
	plotType = 0; 

	% set a line width: 
	linewidth = 1.5; 

	% define which dependent variables should be plotted
	%   1 	-> Plot.
	%   other 	-> Do not plot.
	plotControl=[ 
		1	% e0_c_i1  
		1	% e0_c_i2  
		1	% e0_c_i3  
		1	% e0_c_i4  
		1	% e0_c_i5  
		1	% e0_n_i1  
		1	% e0_n_i2  
		1	% e0_n_i3  
		1	% e0_n_i4  
		1	% e0_n_i5  
		1	% e0_n_L  
		1	% e0_greek_alpha  
		1	% e0_greek_gamma  
		1	% e0_X  
		1	% e0_x_i6  
		1	% e0_x_i7  
		1	% e0_c_i6  
		1	% e0_c_i7  
		1	% e0_greek_psi_cat  
		1	% e0_greek_DeltaG_r3  
		1	% e0_K_eq_r3  
		1	% e0_K_eq_r1  
		1	% e0_r_r1  
		1	% e0_r_r2  
		1	% e0_r_r3  
		1	% e0_r_r4  
		1	% e0_r_r5  
		1	% e0_r_r6  
		1	% e0_r_i1  
		1	% e0_r_i2  
		1	% e0_r_i3  
		1	% e0_r_i4  
		1	% e0_r_i5  
		];

	% decide wether to normalize the y axis
	%   1 	-> Normalized
	%   other 	-> Individual maximum scale
	axisControl = 2;

	%====================================================

	% labels of the dependent variables
	yAxisLabels=[
		'e0.c_{i=1}       '	% e0_c_i1
		'e0.c_{i=2}       '	% e0_c_i2
		'e0.c_{i=3}       '	% e0_c_i3
		'e0.c_{i=4}       '	% e0_c_i4
		'e0.c_{i=5}       '	% e0_c_i5
		'e0.n_{i=1}       '	% e0_n_i1
		'e0.n_{i=2}       '	% e0_n_i2
		'e0.n_{i=3}       '	% e0_n_i3
		'e0.n_{i=4}       '	% e0_n_i4
		'e0.n_{i=5}       '	% e0_n_i5
		'e0.n^{L}         '	% e0_n_L
		'e0.&alpha;       '	% e0_greek_alpha
		'e0.&gamma;       '	% e0_greek_gamma
		'e0.X             '	% e0_X
		'e0.x_{i=6}       '	% e0_x_i6
		'e0.x_{i=7}       '	% e0_x_i7
		'e0.c_{i=6}       '	% e0_c_i6
		'e0.c_{i=7}       '	% e0_c_i7
		'e0.&psi;_{cat}   '	% e0_greek_psi_cat
		'e0.&Delta;G_{r=3}'	% e0_greek_DeltaG_r3
		'e0.K_{r=3}^{eq}  '	% e0_K_eq_r3
		'e0.K_{r=1}^{eq}  '	% e0_K_eq_r1
		'e0.r_{r=1}       '	% e0_r_r1
		'e0.r_{r=2}       '	% e0_r_r2
		'e0.r_{r=3}       '	% e0_r_r3
		'e0.r_{r=4}       '	% e0_r_r4
		'e0.r_{r=5}       '	% e0_r_r5
		'e0.r_{r=6}       '	% e0_r_r6
		'e0.r_{i=1}       '	% e0_r_i1
		'e0.r_{i=2}       '	% e0_r_i2
		'e0.r_{i=3}       '	% e0_r_i3
		'e0.r_{i=4}       '	% e0_r_i4
		'e0.r_{i=5}       '	% e0_r_i5
		];
	xAxisLabel = 't';

	% plot the variables 
	figureIndex=1; 
	xMinVal = min(X); 
	xMaxVal = max(X); 
	yMinVal = min(min(Y)); 
	yMaxVal = max(max(Y)); 
	if (plotType==0) 
		% create a plot for each state variable individually 
		for i=1:length(Y(1,:))
			if (plotControl(i)==1)
				figure(i)
				plot(X,Y(:,i),'LineWidth',linewidth)
				title('Solution of the equation system');
				xlabel(xAxisLabel);
				ylabel(yAxisLabels(i,:));
				if (axisControl==1)
					axis([xMinVal xMaxVal yMinVal yMaxVal]);
				end
				legend(yAxisLabels(i,:));
				figureIndex=figureIndex+1;
			end
		end
	elseif (plotType==1)
		% use a subplot environment
		firstTime = 1;
		numberOfPlots = sum(plotControl);
		figure(1)
		for i=1:length(Y(1,:))
			if (plotControl(i)==1)
				subplot(numberOfPlots,1,figureIndex);
				plot(X,Y(:,i),'LineWidth',linewidth)
				ylabel(yAxisLabels(i,:));
				legend(yAxisLabels(i,:));
				if (axisControl==1)
					axis([xMinVal xMaxVal yMinVal yMaxVal]);
				end
				figureIndex=figureIndex+1;
				if (firstTime) 
					title('Solution of the equation system');
					firstTime = 0;
				end 
		end 
		end
		xlabel(xAxisLabel);
	elseif (plotType==2)
		% plot in one figure
		colors = [
			'r'	% -> Red
			'g' % -> Green
			'b' % -> Blue
			'c' % -> Cyan
			'm' % -> Magenta
			'y' % -> Yellow
			'k' % -> Black
			];
		colorCtr=1;
		maxColors=7;
		figure(1)
		hold on;
		for i=1:length(Y(1,:))
			if (plotControl(i)==1)
				linespec = colors(colorCtr);
				if (colorCtr<=maxColors)
					colorCtr = colorCtr+1;
				end
				plot(X,Y(:,i),linespec,'LineWidth',linewidth);
			end
		end
		title('Solution of the equation system');
		xlabel(xAxisLabel);
		ylabel(yAxisLabels(i,:));
		legend(yAxisLabels(:,:));
		hold off;
	end


end


