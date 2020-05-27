%*********************************************************
% The namespaces have been normalized. The following
% table shows the attribuation. 
% Normalized Name --> Full Name ---> User-defined Name
% =================================== 
% e0 --> e[0]134269 --> 
%*********************************************************

%*********************************************************
% The variables are named according to the notation
% provided in the Mosaic model.
% 
% The variable names can be read as follows:
% ==========================================
% 	e0_R
% 		R: Ideal Gas Constant [J/(mol K)]
% 	 
% 	e0_HU_i##
% 		HU: Hold-up [g]
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_T
% 		T: Temperature [K]
% 	 
% 	e0_p
% 		p: Pressure [bar]
% 	 
% 	e0_n_i#
% 		n: Amount of moles [mol]
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_greek_rho_u#_i#
% 		&rho;: Density [g/l]
% 		Indices
% 			u: Units: 1,-,4=feed tanks; 5,-,8=feed pumps; 9=mixer inlet; 10=gasphase reactor; 11=liq. react.; 12=mixer dekanter; 13=splitter reactor; 14=oil phase decanter; 15=mix ph. dec.; 16=water ph. dec.; 17=X04; 18=pipe volume mix ph.; 19=pipe v. water; 20=oil splitter; 21-23=recycle pumps; 24=recycle mixer;25=product tank; 26=Membr. feed tank; 27=memb. pump; 28=memb. mod.; 29=water tank
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_V_Reactor
% 		V: Volume [l]
% 		Superscripts
% 			Reactor: Reactor related Variables
% 	 
% 	e0_c_i#_u#
% 		c: Concentration [mol/l] or DIPPR 105 Parameter
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 			u: Units: 1,-,4=feed tanks; 5,-,8=feed pumps; 9=mixer inlet; 10=gasphase reactor; 11=liq. react.; 12=mixer dekanter; 13=splitter reactor; 14=oil phase decanter; 15=mix ph. dec.; 16=water ph. dec.; 17=X04; 18=pipe volume mix ph.; 19=pipe v. water; 20=oil splitter; 21-23=recycle pumps; 24=recycle mixer;25=product tank; 26=Membr. feed tank; 27=memb. pump; 28=memb. mod.; 29=water tank
% 	 
% 	e0_greek_DeltaG_r#
% 		&Delta;G: Gibbs free enthalpy
% 		Indices
% 			r: Number of reactions: 1 = GGW; 2 = HydA; 3 = HydB; 4 = ADD
% 	 
% 	e0_K_eq_r#
% 		K: Inhibition factor
% 		Superscripts
% 			eq: Equilibrium
% 		Indices
% 			r: Number of reactions: 1 = GGW; 2 = HydA; 3 = HydB; 4 = ADD
% 	 
% 	e0_r_r#
% 		r: Reaction rate [mol/(l h)]
% 		Indices
% 			r: Number of reactions: 1 = GGW; 2 = HydA; 3 = HydB; 4 = ADD
% 	 
% 	e0_r_i#
% 		r: Reaction rate [mol/(l h)]
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_M_i#
% 		M: Molar mass [g/mol]
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_E_r#
% 		E: Activation Energy [J/mol]
% 		Indices
% 			r: Number of reactions: 1 = GGW; 2 = HydA; 3 = HydB; 4 = ADD
% 	 
% 	e0_k_ref_r#
% 		k: Preexponential factor
% 		Superscripts
% 			ref: Enhancement factor with reference temperature
% 		Indices
% 			r: Number of reactions: 1 = GGW; 2 = HydA; 3 = HydB; 4 = ADD
% 	 
% 	e0_a_i##
% 		a: DIPPR 105 Parameter
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_b_i##
% 		b: DIPPR 105 Parameter
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_c_i##
% 		c: Concentration [mol/l] or DIPPR 105 Parameter
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_d_i##
% 		d: DIPPR 105 Parameter
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_e_i#
% 		e: Vapor density polynom parameter
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 	 
% 	e0_P_GLE_i#_Sol#_SolP#
% 		P: Model Parameter
% 		Superscripts
% 			GLE: Henry Gas Liquid Equilibrium
% 		Indices
% 			i: Number of components: 1 = UnAl; 2 = UnOl; 3 = Aldols; 4 = DEUI; 5 = DEUA; 6 = H2; 7 = CO; 8 = N2, 9 = Water; 10 = Rhacac; 11 = Sx; 12 = Surfactant; 13 = DEA, 14 = Doca
% 			Sol: Number of Solvent for Gas Solubility
% 			SolP: Gas Solubility model parameter idex
% 	 
%*********************************************************

function[X,Y]=solveEquationSystem()

	% declare interval of independent variable
	X_START=0.0;
	X_END=1.0;
	X_INTERVAL=[X_START X_END];	% e0_t

	% specify the initial values for the state variables 
	Y_INIT(1) = 2.9112144;  	% e0_HU_i1  
	Y_INIT(2) = 0.0;  	% e0_HU_i2  
	Y_INIT(3) = 0.2175702;  	% e0_HU_i3  
	Y_INIT(4) = 0.0;  	% e0_HU_i4  
	Y_INIT(5) = 0.0;  	% e0_HU_i5  
	Y_INIT(6) = 11.82565;  	% e0_HU_i9  
	Y_INIT(7) = 1.29;  	% e0_HU_i13  
	Y_INIT(8) = 0.017095627;  	% e0_n_i1  
	Y_INIT(9) = 0.0;  	% e0_n_i2  
	Y_INIT(10) = 6.39E-4;  	% e0_n_i3  
	Y_INIT(11) = 0.0;  	% e0_n_i4  
	Y_INIT(12) = 0.0;  	% e0_n_i5  
	Y_INIT(13) = 0.65642333;  	% e0_n_i9  
	Y_INIT(14) = 0.003910874;  	% e0_n_i12  
	Y_INIT(15) = 0.017637407;  	% e0_n_i13  
	Y_INIT(16) = 0.05283861;  	% e0_n_i14  
	Y_INIT(17) = 0.001479374;  	% e0_n_i7  
	Y_INIT(18) = 8.13E-4;  	% e0_n_i6  
	Y_INIT(19) = 0.0016381;  	% e0_HU_i6  
	Y_INIT(20) = 0.041437257;  	% e0_HU_i7  
	Y_INIT(21) = 5.03E-5;  	% e0_n_i10  
	Y_INIT(22) = 2.02E-4;  	% e0_n_i11  
	Y_INIT(23) = 772.3244306654833;  	% e0_greek_rho_u9_i1  
	Y_INIT(24) = 776.3104270331723;  	% e0_greek_rho_u9_i2  
	Y_INIT(25) = 1090.546648999151;  	% e0_greek_rho_u9_i3  
	Y_INIT(26) = 910.3398759323168;  	% e0_greek_rho_u9_i4  
	Y_INIT(27) = 914.3602790411466;  	% e0_greek_rho_u9_i5  
	Y_INIT(28) = 961.0597510718158;  	% e0_greek_rho_u9_i9  
	Y_INIT(29) = 16547.68197801902;  	% e0_greek_rho_u9_i10  
	Y_INIT(30) = 41861.86377694738;  	% e0_greek_rho_u9_i11  
	Y_INIT(31) = 1121.3142855340068;  	% e0_greek_rho_u9_i12  
	Y_INIT(32) = 615.9736188554265;  	% e0_greek_rho_u9_i13  
	Y_INIT(33) = 690.0862705922896;  	% e0_greek_rho_u9_i14  
	Y_INIT(34) = 0.03341221;  	% e0_V_Reactor  
	Y_INIT(35) = 0.5116581;  	% e0_c_i1_u9  
	Y_INIT(36) = 0.5278731;  	% e0_c_i13_u9  
	Y_INIT(37) = 1.5814161;  	% e0_c_i14_u9  
	Y_INIT(38) = 0.0;  	% e0_c_i2_u9  
	Y_INIT(39) = 0.019119525;  	% e0_c_i3_u9  
	Y_INIT(40) = 0.0;  	% e0_c_i4_u9  
	Y_INIT(41) = 0.0;  	% e0_c_i5_u9  
	Y_INIT(42) = 0.024320394;  	% e0_c_i6_u9  
	Y_INIT(43) = 0.04427644;  	% e0_c_i7_u9  
	Y_INIT(44) = 19.64621;  	% e0_c_i9_u9  
	Y_INIT(45) = 0.001505191;  	% e0_c_i10_u9  
	Y_INIT(46) = 502.1222;  	% e0_greek_DeltaG_r1  
	Y_INIT(47) = 1.1757618;  	% e0_K_eq_r1  
	Y_INIT(48) = 1351500.2;  	% e0_r_r1  
	Y_INIT(49) = 0.0;  	% e0_r_r2  
	Y_INIT(50) = 1.5E-18;  	% e0_r_r3  
	Y_INIT(51) = 0.0;  	% e0_r_r4  
	Y_INIT(52) = -1351500.2;  	% e0_r_i1  
	Y_INIT(53) = 1.5E-18;  	% e0_r_i2  
	Y_INIT(54) = 0.0;  	% e0_r_i3  
	Y_INIT(55) = 1351500.2;  	% e0_r_i4  
	Y_INIT(56) = 0.0;  	% e0_r_i5  
	Y_INIT(57) = -1.5E-18;  	% e0_r_i6  
	Y_INIT(58) = 1351500.2;  	% e0_r_i9  
	Y_INIT(59) = -1351500.2;  	% e0_r_i13  
	Y_INIT(60) = 1.9335168620499728;  	% e0_greek_rho_u9_i6  
	Y_INIT(61) = 27.067704397594195;  	% e0_greek_rho_u9_i8  
	Y_INIT(62) = 27.51435517638899;  	% e0_greek_rho_u9_i7  

	% declare parameters 
	PARAMS(1) = 170.29;  	% e0_M_i1 
	PARAMS(2) = 172.3077;  	% e0_M_i2 
	PARAMS(3) = 340.5784;  	% e0_M_i3 
	PARAMS(4) = 226.43;  	% e0_M_i4 
	PARAMS(5) = 227.43;  	% e0_M_i5 
	PARAMS(6) = 18.01528;  	% e0_M_i9 
	PARAMS(7) = 2.24;  	% e0_HU_i12 
	PARAMS(8) = 572.762;  	% e0_M_i12 
	PARAMS(9) = 73.14;  	% e0_M_i13 
	PARAMS(10) = 9.0;  	% e0_HU_i14 
	PARAMS(11) = 170.33;  	% e0_M_i14 
	PARAMS(12) = 373.0;  	% e0_T 
	PARAMS(13) = 30.0;  	% e0_p 
	PARAMS(14) = -8209.592;  	% e0_P_GLE_i7_Sol1_SolP1 
	PARAMS(15) = 117.0208;  	% e0_P_GLE_i7_Sol1_SolP2 
	PARAMS(16) = -0.20219;  	% e0_P_GLE_i7_Sol1_SolP3 
	PARAMS(17) = -8209.592;  	% e0_P_GLE_i7_Sol2_SolP1 
	PARAMS(18) = 117.0208;  	% e0_P_GLE_i7_Sol2_SolP2 
	PARAMS(19) = -0.20219;  	% e0_P_GLE_i7_Sol2_SolP3 
	PARAMS(20) = -8209.592;  	% e0_P_GLE_i7_Sol3_SolP1 
	PARAMS(21) = 117.0208;  	% e0_P_GLE_i7_Sol3_SolP2 
	PARAMS(22) = -0.20219;  	% e0_P_GLE_i7_Sol3_SolP3 
	PARAMS(23) = -8209.592;  	% e0_P_GLE_i7_Sol4_SolP1 
	PARAMS(24) = 117.0208;  	% e0_P_GLE_i7_Sol4_SolP2 
	PARAMS(25) = -0.20219;  	% e0_P_GLE_i7_Sol4_SolP3 
	PARAMS(26) = -8209.592;  	% e0_P_GLE_i7_Sol5_SolP1 
	PARAMS(27) = 117.0208;  	% e0_P_GLE_i7_Sol5_SolP2 
	PARAMS(28) = -0.20219;  	% e0_P_GLE_i7_Sol5_SolP3 
	PARAMS(29) = 406189.431;  	% e0_P_GLE_i7_Sol9_SolP1 
	PARAMS(30) = -1243.066;  	% e0_P_GLE_i7_Sol9_SolP2 
	PARAMS(31) = 0.95534;  	% e0_P_GLE_i7_Sol9_SolP3 
	PARAMS(32) = 406189.431;  	% e0_P_GLE_i7_Sol12_SolP1 
	PARAMS(33) = -1243.066;  	% e0_P_GLE_i7_Sol12_SolP2 
	PARAMS(34) = 0.95534;  	% e0_P_GLE_i7_Sol12_SolP3 
	PARAMS(35) = 406189.431;  	% e0_P_GLE_i7_Sol13_SolP1 
	PARAMS(36) = -1243.066;  	% e0_P_GLE_i7_Sol13_SolP2 
	PARAMS(37) = 0.95534;  	% e0_P_GLE_i7_Sol13_SolP3 
	PARAMS(38) = -790.7257;  	% e0_P_GLE_i7_Sol14_SolP1 
	PARAMS(39) = 8.3594;  	% e0_P_GLE_i7_Sol14_SolP2 
	PARAMS(40) = -0.012265;  	% e0_P_GLE_i7_Sol14_SolP3 
	PARAMS(41) = 196006.2531;  	% e0_P_GLE_i6_Sol1_SolP1 
	PARAMS(42) = -4898.8245;  	% e0_P_GLE_i6_Sol1_SolP2 
	PARAMS(43) = 30.4876;  	% e0_P_GLE_i6_Sol1_SolP3 
	PARAMS(44) = 196006.2531;  	% e0_P_GLE_i6_Sol2_SolP1 
	PARAMS(45) = -4898.8245;  	% e0_P_GLE_i6_Sol2_SolP2 
	PARAMS(46) = 30.4876;  	% e0_P_GLE_i6_Sol2_SolP3 
	PARAMS(47) = 196006.2531;  	% e0_P_GLE_i6_Sol3_SolP1 
	PARAMS(48) = -4898.8245;  	% e0_P_GLE_i6_Sol3_SolP2 
	PARAMS(49) = 30.4876;  	% e0_P_GLE_i6_Sol3_SolP3 
	PARAMS(50) = 89706.8526;  	% e0_P_GLE_i6_Sol4_SolP1 
	PARAMS(51) = -2230.584;  	% e0_P_GLE_i6_Sol4_SolP2 
	PARAMS(52) = 13.9085;  	% e0_P_GLE_i6_Sol4_SolP3 
	PARAMS(53) = 89706.8526;  	% e0_P_GLE_i6_Sol5_SolP1 
	PARAMS(54) = -2230.584;  	% e0_P_GLE_i6_Sol5_SolP2 
	PARAMS(55) = 13.9085;  	% e0_P_GLE_i6_Sol5_SolP3 
	PARAMS(56) = -547307.7898;  	% e0_P_GLE_i6_Sol9_SolP1 
	PARAMS(57) = 3779.8068;  	% e0_P_GLE_i6_Sol9_SolP2 
	PARAMS(58) = -5.7113;  	% e0_P_GLE_i6_Sol9_SolP3 
	PARAMS(59) = -547307.7898;  	% e0_P_GLE_i6_Sol12_SolP1 
	PARAMS(60) = 3779.8068;  	% e0_P_GLE_i6_Sol12_SolP2 
	PARAMS(61) = -5.7113;  	% e0_P_GLE_i6_Sol12_SolP3 
	PARAMS(62) = -1086.9672;  	% e0_P_GLE_i6_Sol13_SolP1 
	PARAMS(63) = 34.3488;  	% e0_P_GLE_i6_Sol13_SolP2 
	PARAMS(64) = -0.21305;  	% e0_P_GLE_i6_Sol13_SolP3 
	PARAMS(65) = 18109.9827;  	% e0_P_GLE_i6_Sol14_SolP1 
	PARAMS(66) = -83.742;  	% e0_P_GLE_i6_Sol14_SolP2 
	PARAMS(67) = 0.10281;  	% e0_P_GLE_i6_Sol14_SolP3 
	PARAMS(68) = 2.01588;  	% e0_M_i6 
	PARAMS(69) = 28.01;  	% e0_M_i7 
	PARAMS(70) = 0.0156;  	% e0_HU_i10 
	PARAMS(71) = 310.19;  	% e0_M_i10 
	PARAMS(72) = 0.15875;  	% e0_HU_i11 
	PARAMS(73) = 784.71;  	% e0_M_i11 
	PARAMS(74) = 8.314;  	% e0_R 
	PARAMS(75) = 5429.221436274177;  	% e0_E_r1 
	PARAMS(76) = 2.881630631935921E7;  	% e0_k_ref_r1 
	PARAMS(77) = 14458.38389692845;  	% e0_E_r2 
	PARAMS(78) = 151140.1709572773;  	% e0_k_ref_r2 
	PARAMS(79) = 96644.83972725086;  	% e0_E_r3 
	PARAMS(80) = 2.708134161647571;  	% e0_k_ref_r3 
	PARAMS(81) = 33436.13818100635;  	% e0_E_r4 
	PARAMS(82) = 4231.062893490289;  	% e0_k_ref_r4 
	PARAMS(83) = 1.9385;  	% e0_a_i14 
	PARAMS(84) = 0.58748;  	% e0_b_i14 
	PARAMS(85) = 506.0108;  	% e0_c_i14 
	PARAMS(86) = 0.71269;  	% e0_d_i14 
	PARAMS(87) = 1.505;  	% e0_a_i1 
	PARAMS(88) = 0.50133;  	% e0_b_i1 
	PARAMS(89) = 598.4108;  	% e0_c_i1 
	PARAMS(90) = 0.52735;  	% e0_d_i1 
	PARAMS(91) = 0.02848;  	% e0_a_i9 
	PARAMS(92) = 0.02203;  	% e0_b_i9 
	PARAMS(93) = 447.246;  	% e0_c_i9 
	PARAMS(94) = 0.014091;  	% e0_d_i9 
	PARAMS(95) = 0.56598;  	% e0_a_i2 
	PARAMS(96) = 0.31081;  	% e0_b_i2 
	PARAMS(97) = 669.5672;  	% e0_c_i2 
	PARAMS(98) = 0.31266;  	% e0_d_i2 
	PARAMS(99) = 0.02848;  	% e0_a_i10 
	PARAMS(100) = 0.02203;  	% e0_b_i10 
	PARAMS(101) = 447.246;  	% e0_c_i10 
	PARAMS(102) = 0.014091;  	% e0_d_i10 
	PARAMS(103) = 0.25878;  	% e0_a_i5 
	PARAMS(104) = 0.21792;  	% e0_b_i5 
	PARAMS(105) = 692.8852;  	% e0_c_i5 
	PARAMS(106) = 0.28804;  	% e0_d_i5 
	PARAMS(107) = 1.9404E-10;  	% e0_a_i13 
	PARAMS(108) = 3.4576E-6;  	% e0_b_i13 
	PARAMS(109) = -757347.676;  	% e0_c_i13 
	PARAMS(110) = -108.8259;  	% e0_d_i13 
	PARAMS(111) = 0.25878;  	% e0_a_i4 
	PARAMS(112) = 0.21792;  	% e0_b_i4 
	PARAMS(113) = 692.8852;  	% e0_c_i4 
	PARAMS(114) = 0.28804;  	% e0_d_i4 
	PARAMS(115) = 0.097181;  	% e0_a_i12 
	PARAMS(116) = 0.20184;  	% e0_b_i12 
	PARAMS(117) = 1069.2551;  	% e0_c_i12 
	PARAMS(118) = 0.3072;  	% e0_d_i12 
	PARAMS(119) = 0.13876;  	% e0_a_i3 
	PARAMS(120) = 0.18201;  	% e0_b_i3 
	PARAMS(121) = 816.1091;  	% e0_c_i3 
	PARAMS(122) = 0.28091;  	% e0_d_i3 
	PARAMS(123) = 0.02848;  	% e0_a_i11 
	PARAMS(124) = 0.02203;  	% e0_b_i11 
	PARAMS(125) = 447.246;  	% e0_c_i11 
	PARAMS(126) = 0.014091;  	% e0_d_i11 
	PARAMS(127) = 0.50448;  	% e0_a_i6 
	PARAMS(128) = -0.0027972;  	% e0_b_i6 
	PARAMS(129) = 0.066184;  	% e0_c_i6 
	PARAMS(130) = 3.8561E-6;  	% e0_d_i6 
	PARAMS(131) = -9.151E-5;  	% e0_e_i6 
	PARAMS(132) = 28.0134;  	% e0_M_i8 
	PARAMS(133) = 0.55105;  	% e0_a_i8 
	PARAMS(134) = -0.0030955;  	% e0_b_i8 
	PARAMS(135) = 0.068661;  	% e0_c_i8 
	PARAMS(136) = 4.307E-6;  	% e0_d_i8 
	PARAMS(137) = -9.7341E-5;  	% e0_e_i8 
	PARAMS(138) = -0.5581;  	% e0_a_i7 
	PARAMS(139) = 0.0031211;  	% e0_b_i7 
	PARAMS(140) = 0.063578;  	% e0_c_i7 
	PARAMS(141) = -4.3159E-6;  	% e0_d_i7 
	PARAMS(142) = -8.3167E-5;  	% e0_e_i7 

	M = daeSystemMM();

	OPTIONS = odeset('Mass',M);

	[X,Y] = ode15s(@(X,Y)daeSystemLHS(X,Y,PARAMS),X_INTERVAL,Y_INIT',OPTIONS);

	displayResults(X,Y);

end


function MASS = daeSystemMM()

	MASS = zeros(62, 62); %total number of equations
	MASS(1:7, 1:7)=eye(7,7); % number of odes

end


% evaluate the differential function.
function[DYDX] = daeSystemLHS(X,Y,PARAMS)

	% read out variables  
	e0_HU_i1 = Y(1); 
	e0_HU_i2 = Y(2); 
	e0_HU_i3 = Y(3); 
	e0_HU_i4 = Y(4); 
	e0_HU_i5 = Y(5); 
	e0_HU_i9 = Y(6); 
	e0_HU_i13 = Y(7); 
	e0_n_i1 = Y(8); 
	e0_n_i2 = Y(9); 
	e0_n_i3 = Y(10); 
	e0_n_i4 = Y(11); 
	e0_n_i5 = Y(12); 
	e0_n_i9 = Y(13); 
	e0_n_i12 = Y(14); 
	e0_n_i13 = Y(15); 
	e0_n_i14 = Y(16); 
	e0_n_i7 = Y(17); 
	e0_n_i6 = Y(18); 
	e0_HU_i6 = Y(19); 
	e0_HU_i7 = Y(20); 
	e0_n_i10 = Y(21); 
	e0_n_i11 = Y(22); 
	e0_greek_rho_u9_i1 = Y(23); 
	e0_greek_rho_u9_i2 = Y(24); 
	e0_greek_rho_u9_i3 = Y(25); 
	e0_greek_rho_u9_i4 = Y(26); 
	e0_greek_rho_u9_i5 = Y(27); 
	e0_greek_rho_u9_i9 = Y(28); 
	e0_greek_rho_u9_i10 = Y(29); 
	e0_greek_rho_u9_i11 = Y(30); 
	e0_greek_rho_u9_i12 = Y(31); 
	e0_greek_rho_u9_i13 = Y(32); 
	e0_greek_rho_u9_i14 = Y(33); 
	e0_V_Reactor = Y(34); 
	e0_c_i1_u9 = Y(35); 
	e0_c_i13_u9 = Y(36); 
	e0_c_i14_u9 = Y(37); 
	e0_c_i2_u9 = Y(38); 
	e0_c_i3_u9 = Y(39); 
	e0_c_i4_u9 = Y(40); 
	e0_c_i5_u9 = Y(41); 
	e0_c_i6_u9 = Y(42); 
	e0_c_i7_u9 = Y(43); 
	e0_c_i9_u9 = Y(44); 
	e0_c_i10_u9 = Y(45); 
	e0_greek_DeltaG_r1 = Y(46); 
	e0_K_eq_r1 = Y(47); 
	e0_r_r1 = Y(48); 
	e0_r_r2 = Y(49); 
	e0_r_r3 = Y(50); 
	e0_r_r4 = Y(51); 
	e0_r_i1 = Y(52); 
	e0_r_i2 = Y(53); 
	e0_r_i3 = Y(54); 
	e0_r_i4 = Y(55); 
	e0_r_i5 = Y(56); 
	e0_r_i6 = Y(57); 
	e0_r_i9 = Y(58); 
	e0_r_i13 = Y(59); 
	e0_greek_rho_u9_i6 = Y(60); 
	e0_greek_rho_u9_i8 = Y(61); 
	e0_greek_rho_u9_i7 = Y(62); 

	% read out differential variable
	e0_t = X;

	% read out parameters 
	e0_M_i1 = PARAMS(1); 
	e0_M_i2 = PARAMS(2); 
	e0_M_i3 = PARAMS(3); 
	e0_M_i4 = PARAMS(4); 
	e0_M_i5 = PARAMS(5); 
	e0_M_i9 = PARAMS(6); 
	e0_HU_i12 = PARAMS(7); 
	e0_M_i12 = PARAMS(8); 
	e0_M_i13 = PARAMS(9); 
	e0_HU_i14 = PARAMS(10); 
	e0_M_i14 = PARAMS(11); 
	e0_T = PARAMS(12); 
	e0_p = PARAMS(13); 
	e0_P_GLE_i7_Sol1_SolP1 = PARAMS(14); 
	e0_P_GLE_i7_Sol1_SolP2 = PARAMS(15); 
	e0_P_GLE_i7_Sol1_SolP3 = PARAMS(16); 
	e0_P_GLE_i7_Sol2_SolP1 = PARAMS(17); 
	e0_P_GLE_i7_Sol2_SolP2 = PARAMS(18); 
	e0_P_GLE_i7_Sol2_SolP3 = PARAMS(19); 
	e0_P_GLE_i7_Sol3_SolP1 = PARAMS(20); 
	e0_P_GLE_i7_Sol3_SolP2 = PARAMS(21); 
	e0_P_GLE_i7_Sol3_SolP3 = PARAMS(22); 
	e0_P_GLE_i7_Sol4_SolP1 = PARAMS(23); 
	e0_P_GLE_i7_Sol4_SolP2 = PARAMS(24); 
	e0_P_GLE_i7_Sol4_SolP3 = PARAMS(25); 
	e0_P_GLE_i7_Sol5_SolP1 = PARAMS(26); 
	e0_P_GLE_i7_Sol5_SolP2 = PARAMS(27); 
	e0_P_GLE_i7_Sol5_SolP3 = PARAMS(28); 
	e0_P_GLE_i7_Sol9_SolP1 = PARAMS(29); 
	e0_P_GLE_i7_Sol9_SolP2 = PARAMS(30); 
	e0_P_GLE_i7_Sol9_SolP3 = PARAMS(31); 
	e0_P_GLE_i7_Sol12_SolP1 = PARAMS(32); 
	e0_P_GLE_i7_Sol12_SolP2 = PARAMS(33); 
	e0_P_GLE_i7_Sol12_SolP3 = PARAMS(34); 
	e0_P_GLE_i7_Sol13_SolP1 = PARAMS(35); 
	e0_P_GLE_i7_Sol13_SolP2 = PARAMS(36); 
	e0_P_GLE_i7_Sol13_SolP3 = PARAMS(37); 
	e0_P_GLE_i7_Sol14_SolP1 = PARAMS(38); 
	e0_P_GLE_i7_Sol14_SolP2 = PARAMS(39); 
	e0_P_GLE_i7_Sol14_SolP3 = PARAMS(40); 
	e0_P_GLE_i6_Sol1_SolP1 = PARAMS(41); 
	e0_P_GLE_i6_Sol1_SolP2 = PARAMS(42); 
	e0_P_GLE_i6_Sol1_SolP3 = PARAMS(43); 
	e0_P_GLE_i6_Sol2_SolP1 = PARAMS(44); 
	e0_P_GLE_i6_Sol2_SolP2 = PARAMS(45); 
	e0_P_GLE_i6_Sol2_SolP3 = PARAMS(46); 
	e0_P_GLE_i6_Sol3_SolP1 = PARAMS(47); 
	e0_P_GLE_i6_Sol3_SolP2 = PARAMS(48); 
	e0_P_GLE_i6_Sol3_SolP3 = PARAMS(49); 
	e0_P_GLE_i6_Sol4_SolP1 = PARAMS(50); 
	e0_P_GLE_i6_Sol4_SolP2 = PARAMS(51); 
	e0_P_GLE_i6_Sol4_SolP3 = PARAMS(52); 
	e0_P_GLE_i6_Sol5_SolP1 = PARAMS(53); 
	e0_P_GLE_i6_Sol5_SolP2 = PARAMS(54); 
	e0_P_GLE_i6_Sol5_SolP3 = PARAMS(55); 
	e0_P_GLE_i6_Sol9_SolP1 = PARAMS(56); 
	e0_P_GLE_i6_Sol9_SolP2 = PARAMS(57); 
	e0_P_GLE_i6_Sol9_SolP3 = PARAMS(58); 
	e0_P_GLE_i6_Sol12_SolP1 = PARAMS(59); 
	e0_P_GLE_i6_Sol12_SolP2 = PARAMS(60); 
	e0_P_GLE_i6_Sol12_SolP3 = PARAMS(61); 
	e0_P_GLE_i6_Sol13_SolP1 = PARAMS(62); 
	e0_P_GLE_i6_Sol13_SolP2 = PARAMS(63); 
	e0_P_GLE_i6_Sol13_SolP3 = PARAMS(64); 
	e0_P_GLE_i6_Sol14_SolP1 = PARAMS(65); 
	e0_P_GLE_i6_Sol14_SolP2 = PARAMS(66); 
	e0_P_GLE_i6_Sol14_SolP3 = PARAMS(67); 
	e0_M_i6 = PARAMS(68); 
	e0_M_i7 = PARAMS(69); 
	e0_HU_i10 = PARAMS(70); 
	e0_M_i10 = PARAMS(71); 
	e0_HU_i11 = PARAMS(72); 
	e0_M_i11 = PARAMS(73); 
	e0_R = PARAMS(74); 
	e0_E_r1 = PARAMS(75); 
	e0_k_ref_r1 = PARAMS(76); 
	e0_E_r2 = PARAMS(77); 
	e0_k_ref_r2 = PARAMS(78); 
	e0_E_r3 = PARAMS(79); 
	e0_k_ref_r3 = PARAMS(80); 
	e0_E_r4 = PARAMS(81); 
	e0_k_ref_r4 = PARAMS(82); 
	e0_a_i14 = PARAMS(83); 
	e0_b_i14 = PARAMS(84); 
	e0_c_i14 = PARAMS(85); 
	e0_d_i14 = PARAMS(86); 
	e0_a_i1 = PARAMS(87); 
	e0_b_i1 = PARAMS(88); 
	e0_c_i1 = PARAMS(89); 
	e0_d_i1 = PARAMS(90); 
	e0_a_i9 = PARAMS(91); 
	e0_b_i9 = PARAMS(92); 
	e0_c_i9 = PARAMS(93); 
	e0_d_i9 = PARAMS(94); 
	e0_a_i2 = PARAMS(95); 
	e0_b_i2 = PARAMS(96); 
	e0_c_i2 = PARAMS(97); 
	e0_d_i2 = PARAMS(98); 
	e0_a_i10 = PARAMS(99); 
	e0_b_i10 = PARAMS(100); 
	e0_c_i10 = PARAMS(101); 
	e0_d_i10 = PARAMS(102); 
	e0_a_i5 = PARAMS(103); 
	e0_b_i5 = PARAMS(104); 
	e0_c_i5 = PARAMS(105); 
	e0_d_i5 = PARAMS(106); 
	e0_a_i13 = PARAMS(107); 
	e0_b_i13 = PARAMS(108); 
	e0_c_i13 = PARAMS(109); 
	e0_d_i13 = PARAMS(110); 
	e0_a_i4 = PARAMS(111); 
	e0_b_i4 = PARAMS(112); 
	e0_c_i4 = PARAMS(113); 
	e0_d_i4 = PARAMS(114); 
	e0_a_i12 = PARAMS(115); 
	e0_b_i12 = PARAMS(116); 
	e0_c_i12 = PARAMS(117); 
	e0_d_i12 = PARAMS(118); 
	e0_a_i3 = PARAMS(119); 
	e0_b_i3 = PARAMS(120); 
	e0_c_i3 = PARAMS(121); 
	e0_d_i3 = PARAMS(122); 
	e0_a_i11 = PARAMS(123); 
	e0_b_i11 = PARAMS(124); 
	e0_c_i11 = PARAMS(125); 
	e0_d_i11 = PARAMS(126); 
	e0_a_i6 = PARAMS(127); 
	e0_b_i6 = PARAMS(128); 
	e0_c_i6 = PARAMS(129); 
	e0_d_i6 = PARAMS(130); 
	e0_e_i6 = PARAMS(131); 
	e0_M_i8 = PARAMS(132); 
	e0_a_i8 = PARAMS(133); 
	e0_b_i8 = PARAMS(134); 
	e0_c_i8 = PARAMS(135); 
	e0_d_i8 = PARAMS(136); 
	e0_e_i8 = PARAMS(137); 
	e0_a_i7 = PARAMS(138); 
	e0_b_i7 = PARAMS(139); 
	e0_c_i7 = PARAMS(140); 
	e0_d_i7 = PARAMS(141); 
	e0_e_i7 = PARAMS(142); 

	% evaluate the function values  
	DYDX(1) = e0_r_i1 * 60.0 * e0_V_Reactor * e0_M_i1; 
	DYDX(2) = e0_r_i2 * 60.0 * e0_V_Reactor * e0_M_i2; 
	DYDX(3) = e0_r_i3 * 60.0 * e0_V_Reactor * e0_M_i3; 
	DYDX(4) = e0_r_i4 * 60.0 * e0_V_Reactor * e0_M_i4; 
	DYDX(5) = e0_r_i5 * 60.0 * e0_V_Reactor * e0_M_i5; 
	DYDX(6) = e0_r_i9 * 60.0 * e0_V_Reactor * e0_M_i9; 
	DYDX(7) = e0_r_i13 * 60.0 * e0_V_Reactor * e0_M_i13; 
	DYDX(8) = e0_n_i1 * e0_M_i1 - (e0_HU_i1); 
	DYDX(9) = e0_n_i2 * e0_M_i2 - (e0_HU_i2); 
	DYDX(10) = e0_n_i3 * e0_M_i3 - (e0_HU_i3); 
	DYDX(11) = e0_n_i4 * e0_M_i4 - (e0_HU_i4); 
	DYDX(12) = e0_n_i5 * e0_M_i5 - (e0_HU_i5); 
	DYDX(13) = e0_n_i9 * e0_M_i9 - (e0_HU_i9); 
	DYDX(14) = e0_n_i12 * e0_M_i12 - (e0_HU_i12); 
	DYDX(15) = e0_n_i13 * e0_M_i13 - (e0_HU_i13); 
	DYDX(16) = e0_n_i14 * e0_M_i14 - (e0_HU_i14); 
	DYDX(17) = e0_n_i7 - ((e0_n_i1)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol1_SolP1 + e0_T * e0_P_GLE_i7_Sol1_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol1_SolP3))) - e0_n_i1 + (e0_n_i2)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol2_SolP1 + e0_T * e0_P_GLE_i7_Sol2_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol2_SolP3))) - e0_n_i2 + (e0_n_i3)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol3_SolP1 + e0_T * e0_P_GLE_i7_Sol3_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol3_SolP3))) - e0_n_i3 + (e0_n_i4)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol4_SolP1 + e0_T * e0_P_GLE_i7_Sol4_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol4_SolP3))) - e0_n_i4 + (e0_n_i5)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol5_SolP1 + e0_T * e0_P_GLE_i7_Sol5_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol5_SolP3))) - e0_n_i5 + (e0_n_i9)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol9_SolP1 + e0_T * e0_P_GLE_i7_Sol9_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol9_SolP3))) - e0_n_i9 + (e0_n_i12)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol12_SolP1 + e0_T * e0_P_GLE_i7_Sol12_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol12_SolP3))) - e0_n_i12 + (e0_n_i13)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol13_SolP1 + e0_T * e0_P_GLE_i7_Sol13_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol13_SolP3))) - e0_n_i13 + (e0_n_i14)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i7_Sol14_SolP1 + e0_T * e0_P_GLE_i7_Sol14_SolP2 + power((e0_T),2.0) * e0_P_GLE_i7_Sol14_SolP3))) - e0_n_i14); 
	DYDX(18) = e0_n_i6 - ((e0_n_i1)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol1_SolP1 + e0_T * e0_P_GLE_i6_Sol1_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol1_SolP3))) - e0_n_i1 + (e0_n_i2)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol2_SolP1 + e0_T * e0_P_GLE_i6_Sol2_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol2_SolP3))) - e0_n_i2 + (e0_n_i3)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol3_SolP1 + e0_T * e0_P_GLE_i6_Sol3_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol3_SolP3))) - e0_n_i3 + (e0_n_i4)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol4_SolP1 + e0_T * e0_P_GLE_i6_Sol4_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol4_SolP3))) - e0_n_i4 + (e0_n_i5)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol5_SolP1 + e0_T * e0_P_GLE_i6_Sol5_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol5_SolP3))) - e0_n_i5 + (e0_n_i9)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol9_SolP1 + e0_T * e0_P_GLE_i6_Sol9_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol9_SolP3))) - e0_n_i9 + (e0_n_i12)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol12_SolP1 + e0_T * e0_P_GLE_i6_Sol12_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol12_SolP3))) - e0_n_i12 + (e0_n_i13)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol13_SolP1 + e0_T * e0_P_GLE_i6_Sol13_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol13_SolP3))) - e0_n_i13 + (e0_n_i14)/(1.0 - (e0_p)/(2.0 * (e0_P_GLE_i6_Sol14_SolP1 + e0_T * e0_P_GLE_i6_Sol14_SolP2 + power((e0_T),2.0) * e0_P_GLE_i6_Sol14_SolP3))) - e0_n_i14); 
	DYDX(19) = e0_n_i6 * e0_M_i6 - (e0_HU_i6); 
	DYDX(20) = e0_n_i7 * e0_M_i7 - (e0_HU_i7); 
	DYDX(21) = e0_n_i10 * e0_M_i10 - (e0_HU_i10); 
	DYDX(22) = e0_n_i11 * e0_M_i11 - (e0_HU_i11); 
	DYDX(23) = e0_V_Reactor - (((e0_HU_i1)/(e0_greek_rho_u9_i1)) + ((e0_HU_i2)/(e0_greek_rho_u9_i2)) + ((e0_HU_i3)/(e0_greek_rho_u9_i3)) + ((e0_HU_i4)/(e0_greek_rho_u9_i4)) + ((e0_HU_i5)/(e0_greek_rho_u9_i5)) + ((e0_HU_i9)/(e0_greek_rho_u9_i9)) + ((e0_HU_i10)/(e0_greek_rho_u9_i10)) + ((e0_HU_i11)/(e0_greek_rho_u9_i11)) + ((e0_HU_i12)/(e0_greek_rho_u9_i12)) + ((e0_HU_i13)/(e0_greek_rho_u9_i13)) + ((e0_HU_i14)/(e0_greek_rho_u9_i14))); 
	DYDX(24) = e0_n_i1 - (e0_c_i1_u9 * e0_V_Reactor); 
	DYDX(25) = e0_n_i13 - (e0_c_i13_u9 * e0_V_Reactor); 
	DYDX(26) = e0_n_i14 - (e0_c_i14_u9 * e0_V_Reactor); 
	DYDX(27) = e0_n_i2 - (e0_c_i2_u9 * e0_V_Reactor); 
	DYDX(28) = e0_n_i3 - (e0_c_i3_u9 * e0_V_Reactor); 
	DYDX(29) = e0_n_i4 - (e0_c_i4_u9 * e0_V_Reactor); 
	DYDX(30) = e0_n_i5 - (e0_c_i5_u9 * e0_V_Reactor); 
	DYDX(31) = e0_n_i6 - (e0_c_i6_u9 * e0_V_Reactor); 
	DYDX(32) = e0_n_i7 - (e0_c_i7_u9 * e0_V_Reactor); 
	DYDX(33) = e0_n_i9 - (e0_c_i9_u9 * e0_V_Reactor); 
	DYDX(34) = e0_n_i10 - (e0_c_i10_u9 * e0_V_Reactor); 
	DYDX(35) = e0_greek_DeltaG_r1 - ( - 163320.0 + 777.96 * e0_T - 0.9082 * power((e0_T),2.0)); 
	DYDX(36) = e0_K_eq_r1 - (exp((e0_greek_DeltaG_r1)/(e0_R * e0_T))); 
	DYDX(37) = e0_r_r1 - (e0_k_ref_r1 * exp( - (e0_E_r1)/(e0_R * e0_T)) * (e0_c_i1_u9 * e0_c_i13_u9 - ((e0_c_i4_u9 * e0_c_i9_u9)/(e0_K_eq_r1)))); 
	DYDX(38) = e0_r_r2 - (e0_c_i10_u9 * e0_k_ref_r2 * exp( - (e0_E_r2)/(e0_R * e0_T)) * e0_c_i4_u9 * e0_c_i6_u9); 
	DYDX(39) = e0_r_r3 - (e0_c_i10_u9 * e0_k_ref_r3 * exp( - (e0_E_r3)/(e0_R * e0_T)) * e0_c_i1_u9 * e0_c_i6_u9); 
	DYDX(40) = e0_r_r4 - (e0_k_ref_r4 * exp( - (e0_E_r4)/(e0_R * e0_T)) * e0_c_i1_u9 * e0_c_i4_u9); 
	DYDX(41) = e0_r_i1 - ( - e0_r_r1 - e0_r_r3 - e0_r_r4); 
	DYDX(42) = e0_r_i2 - (e0_r_r3); 
	DYDX(43) = e0_r_i3 - (e0_r_r4); 
	DYDX(44) = e0_r_i4 - (e0_r_r1 - e0_r_r2 - e0_r_r4); 
	DYDX(45) = e0_r_i5 - (e0_r_r2); 
	DYDX(46) = e0_r_i6 - ( - e0_r_r2 - e0_r_r3); 
	DYDX(47) = e0_r_i9 - (e0_r_r1); 
	DYDX(48) = e0_r_i13 - ( - e0_r_r1 + e0_r_r4); 
	DYDX(49) = (e0_a_i14)/(power((e0_b_i14),1.0 + power((1.0 - (e0_T)/(e0_c_i14)),e0_d_i14))) * e0_M_i14 - (e0_greek_rho_u9_i14); 
	DYDX(50) = (e0_a_i1)/(power((e0_b_i1),1.0 + power((1.0 - (e0_T)/(e0_c_i1)),e0_d_i1))) * e0_M_i1 - (e0_greek_rho_u9_i1); 
	DYDX(51) = (e0_a_i9)/(power((e0_b_i9),1.0 + power((1.0 - (e0_T)/(e0_c_i9)),e0_d_i9))) * e0_M_i9 - (e0_greek_rho_u9_i9); 
	DYDX(52) = (e0_a_i2)/(power((e0_b_i2),1.0 + power((1.0 - (e0_T)/(e0_c_i2)),e0_d_i2))) * e0_M_i2 - (e0_greek_rho_u9_i2); 
	DYDX(53) = (e0_a_i10)/(power((e0_b_i10),1.0 + power((1.0 - (e0_T)/(e0_c_i10)),e0_d_i10))) * e0_M_i10 - (e0_greek_rho_u9_i10); 
	DYDX(54) = (e0_a_i5)/(power((e0_b_i5),1.0 + power((1.0 - (e0_T)/(e0_c_i5)),e0_d_i5))) * e0_M_i5 - (e0_greek_rho_u9_i5); 
	DYDX(55) = (e0_a_i13)/(power((e0_b_i13),1.0 + power((1.0 - (e0_T)/(e0_c_i13)),e0_d_i13))) * e0_M_i13 - (e0_greek_rho_u9_i13); 
	DYDX(56) = (e0_a_i4)/(power((e0_b_i4),1.0 + power((1.0 - (e0_T)/(e0_c_i4)),e0_d_i4))) * e0_M_i4 - (e0_greek_rho_u9_i4); 
	DYDX(57) = (e0_a_i12)/(power((e0_b_i12),1.0 + power((1.0 - (e0_T)/(e0_c_i12)),e0_d_i12))) * e0_M_i12 - (e0_greek_rho_u9_i12); 
	DYDX(58) = (e0_a_i3)/(power((e0_b_i3),1.0 + power((1.0 - (e0_T)/(e0_c_i3)),e0_d_i3))) * e0_M_i3 - (e0_greek_rho_u9_i3); 
	DYDX(59) = (e0_a_i11)/(power((e0_b_i11),1.0 + power((1.0 - (e0_T)/(e0_c_i11)),e0_d_i11))) * e0_M_i11 - (e0_greek_rho_u9_i11); 
	DYDX(60) = (e0_a_i6 + e0_b_i6 * e0_T + e0_c_i6 * e0_p + e0_d_i6 * power((e0_T),2.0) + e0_e_i6 * e0_T * e0_p) * e0_M_i6 - (e0_greek_rho_u9_i6); 
	DYDX(61) = (e0_a_i8 + e0_b_i8 * e0_T + e0_c_i8 * e0_p + e0_d_i8 * power((e0_T),2.0) + e0_e_i8 * e0_T * e0_p) * e0_M_i8 - (e0_greek_rho_u9_i8); 
	DYDX(62) = (e0_a_i7 + e0_b_i7 * e0_T + e0_c_i7 * e0_p + e0_d_i7 * power((e0_T),2.0) + e0_e_i7 * e0_T * e0_p) * e0_M_i7 - (e0_greek_rho_u9_i7); 

	DYDX=DYDX';

end

function[] = displayResults(X,Y)

	% decide for a plot type: 
	%   0 	-> Plot the variables into individual figures 
	%   1 	-> Plot into sub figures 
	%   2 	-> Plot all selected into one figure 
	%   other 	-> Do not plot 
	plotType = 1; 

	% set a line width: 
	linewidth = 1.5; 

	% define which dependent variables should be plotted
	%   1 	-> Plot.
	%   other 	-> Do not plot.
	plotControl=[ 
		1	% e0_HU_i1  
		1	% e0_n_i1  
		1	% e0_HU_i2  
		1	% e0_n_i2  
		1	% e0_HU_i3  
		1	% e0_n_i3  
		1	% e0_HU_i4  
		1	% e0_n_i4  
		1	% e0_HU_i5  
		1	% e0_n_i5  
		1	% e0_HU_i9  
		1	% e0_n_i9  
		1	% e0_n_i12  
		1	% e0_HU_i13  
		1	% e0_n_i13  
		1	% e0_n_i14  
		1	% e0_n_i7  
		1	% e0_n_i6  
		1	% e0_HU_i6  
		1	% e0_HU_i7  
		1	% e0_n_i10  
		1	% e0_n_i11  
		1	% e0_greek_rho_u9_i1  
		1	% e0_greek_rho_u9_i2  
		1	% e0_greek_rho_u9_i3  
		1	% e0_greek_rho_u9_i4  
		1	% e0_greek_rho_u9_i5  
		1	% e0_greek_rho_u9_i9  
		1	% e0_greek_rho_u9_i10  
		1	% e0_greek_rho_u9_i11  
		1	% e0_greek_rho_u9_i12  
		1	% e0_greek_rho_u9_i13  
		1	% e0_greek_rho_u9_i14  
		1	% e0_V_Reactor  
		1	% e0_c_i1_u9  
		1	% e0_c_i13_u9  
		1	% e0_c_i14_u9  
		1	% e0_c_i2_u9  
		1	% e0_c_i3_u9  
		1	% e0_c_i4_u9  
		1	% e0_c_i5_u9  
		1	% e0_c_i6_u9  
		1	% e0_c_i7_u9  
		1	% e0_c_i9_u9  
		1	% e0_c_i10_u9  
		1	% e0_greek_DeltaG_r1  
		1	% e0_K_eq_r1  
		1	% e0_r_r1  
		1	% e0_r_r2  
		1	% e0_r_r3  
		1	% e0_r_r4  
		1	% e0_r_i1  
		1	% e0_r_i2  
		1	% e0_r_i3  
		1	% e0_r_i4  
		1	% e0_r_i5  
		1	% e0_r_i6  
		1	% e0_r_i9  
		1	% e0_r_i13  
		1	% e0_greek_rho_u9_i6  
		1	% e0_greek_rho_u9_i8  
		1	% e0_greek_rho_u9_i7  
		];

	% decide wether to normalize the y axis
	%   1 	-> Normalized
	%   other 	-> Individual maximum scale
	axisControl = 1;

	%====================================================

	% labels of the dependent variables
	yAxisLabels=[
		'e0.HU_{i=1}        '	% e0_HU_i1
		'e0.n_{i=1}         '	% e0_n_i1
		'e0.HU_{i=2}        '	% e0_HU_i2
		'e0.n_{i=2}         '	% e0_n_i2
		'e0.HU_{i=3}        '	% e0_HU_i3
		'e0.n_{i=3}         '	% e0_n_i3
		'e0.HU_{i=4}        '	% e0_HU_i4
		'e0.n_{i=4}         '	% e0_n_i4
		'e0.HU_{i=5}        '	% e0_HU_i5
		'e0.n_{i=5}         '	% e0_n_i5
		'e0.HU_{i=9}        '	% e0_HU_i9
		'e0.n_{i=9}         '	% e0_n_i9
		'e0.n_{i=12}        '	% e0_n_i12
		'e0.HU_{i=13}       '	% e0_HU_i13
		'e0.n_{i=13}        '	% e0_n_i13
		'e0.n_{i=14}        '	% e0_n_i14
		'e0.n_{i=7}         '	% e0_n_i7
		'e0.n_{i=6}         '	% e0_n_i6
		'e0.HU_{i=6}        '	% e0_HU_i6
		'e0.HU_{i=7}        '	% e0_HU_i7
		'e0.n_{i=10}        '	% e0_n_i10
		'e0.n_{i=11}        '	% e0_n_i11
		'e0.&rho;_{u=9,i=1} '	% e0_greek_rho_u9_i1
		'e0.&rho;_{u=9,i=2} '	% e0_greek_rho_u9_i2
		'e0.&rho;_{u=9,i=3} '	% e0_greek_rho_u9_i3
		'e0.&rho;_{u=9,i=4} '	% e0_greek_rho_u9_i4
		'e0.&rho;_{u=9,i=5} '	% e0_greek_rho_u9_i5
		'e0.&rho;_{u=9,i=9} '	% e0_greek_rho_u9_i9
		'e0.&rho;_{u=9,i=10}'	% e0_greek_rho_u9_i10
		'e0.&rho;_{u=9,i=11}'	% e0_greek_rho_u9_i11
		'e0.&rho;_{u=9,i=12}'	% e0_greek_rho_u9_i12
		'e0.&rho;_{u=9,i=13}'	% e0_greek_rho_u9_i13
		'e0.&rho;_{u=9,i=14}'	% e0_greek_rho_u9_i14
		'e0.V^{Reactor}     '	% e0_V_Reactor
		'e0.c_{i=1,u=9}     '	% e0_c_i1_u9
		'e0.c_{i=13,u=9}    '	% e0_c_i13_u9
		'e0.c_{i=14,u=9}    '	% e0_c_i14_u9
		'e0.c_{i=2,u=9}     '	% e0_c_i2_u9
		'e0.c_{i=3,u=9}     '	% e0_c_i3_u9
		'e0.c_{i=4,u=9}     '	% e0_c_i4_u9
		'e0.c_{i=5,u=9}     '	% e0_c_i5_u9
		'e0.c_{i=6,u=9}     '	% e0_c_i6_u9
		'e0.c_{i=7,u=9}     '	% e0_c_i7_u9
		'e0.c_{i=9,u=9}     '	% e0_c_i9_u9
		'e0.c_{i=10,u=9}    '	% e0_c_i10_u9
		'e0.&Delta;G_{r=1}  '	% e0_greek_DeltaG_r1
		'e0.K_{r=1}^{eq}    '	% e0_K_eq_r1
		'e0.r_{r=1}         '	% e0_r_r1
		'e0.r_{r=2}         '	% e0_r_r2
		'e0.r_{r=3}         '	% e0_r_r3
		'e0.r_{r=4}         '	% e0_r_r4
		'e0.r_{i=1}         '	% e0_r_i1
		'e0.r_{i=2}         '	% e0_r_i2
		'e0.r_{i=3}         '	% e0_r_i3
		'e0.r_{i=4}         '	% e0_r_i4
		'e0.r_{i=5}         '	% e0_r_i5
		'e0.r_{i=6}         '	% e0_r_i6
		'e0.r_{i=9}         '	% e0_r_i9
		'e0.r_{i=13}        '	% e0_r_i13
		'e0.&rho;_{u=9,i=6} '	% e0_greek_rho_u9_i6
		'e0.&rho;_{u=9,i=8} '	% e0_greek_rho_u9_i8
		'e0.&rho;_{u=9,i=7} '	% e0_greek_rho_u9_i7
		];
	xAxisLabel = 'e0.t';

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


