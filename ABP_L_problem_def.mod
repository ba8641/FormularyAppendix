
/*
 	ABP PROBLEM DEFINITIONS (type: model file)
	#############################################################
	#   IT MUST BE DECLARED AFTER 2 files:                      #
	#                                                           #
	#   ABP_base.mod and ABP_L_H.mod                            #
	#                                                           #
	#############################################################

	by Ruben A. Proano
	rpmeie@rit.edu
	ROCHESTER INSTITUTE OF TECHNOLOGY
	
	NOTE:
	
	# 1- ABPTSS - used for the modified MINLP ABP implementation considering all bundles in B
	              
	# 2- ABPTSS_Lin - used for the parially linearized MINLP ABP implementation considering all bundles in B. 
	
	# 3- ABPFLTSS - used for solving a fully linearized ABP version that results in an MIP formulation. The linear_
	                rization is obtained via the implementation of 6 constraints resulting from 
	                a McCarthy envelop linearization together with additional variable W
		      
	# 4- ABPTCS - used for the LP problem that maximizes total customer surplus. It only considers Y as a variable
	            - options allow for the identification of iis constraints and vars
	
	# 5- ABPTPF - used for the LP problem that maximizes total profit. It only considers Y as a variable
	            - options allow for the identification of iis constraints and vars

	# 6- ABPTSS_LH_com - used in the first stage of the ABP_heuristic (combination vaccines -only).  It does not include 
		             SuperAdditivity constraint. It includes a linearized Elasticity constraint, and constraints to fix 
    		             the allocations and production decisions of each iteration
    		             Make sure loo=1 and hii=card(B)

	# 7- ABPTSS_LH_mono: used in the second stage of the ABP_heuristic (monovalent vaccines -only). It includes the
	             almost the same constraints as ABPTSS_LH_com the exception of the AntigenDemand constraint which
	             includes a deficit variable that is minimized in the objective to increase coverage
 			     Make sure loo=1 and hii=card(B)

	# 8- Environments - after each problem definition you will find the solver choice and the options for the solver. 
	                These group of instructions define the problem environment. Everytime you call a problem 
	                (e.g., solve ABPTSS) in a script, the problem will be solved with the options listed in its
	                 environment. If you drop any constraints or variables, for one problem, and you then call a new
	                 problem, such droped constraints will be automatically restaured.

	VERSION 07.05.17_RP

	07.11.17 - Added support for the problems using multiple entities.

	See additional comments at the end of the file
*/




# ABPTSS: MINLP problem formulation - bundle allocation problem
problem ABPTSS: X,Y,g,TSS, SuperAdditivity, AntigenDemand, Capacity, RecoverAnnuity, Elasticity, SegmentedPrices, ReservationPriceLim;
option solver knitro; 
option knitro_options 'ms_enable=1 ms_maxsolves=20 ma_terminate=0 mip_integral_gap_rel=1e-5 par_numthreads=6';  
option display_eps 1.0e-09; 
option solution_round 6; 

# ABPTSS_Lin: Partially linearized MINLP problem formulation - bundle allocation problem
problem ABPTSS_Lin: X,Y,g,TSS, SuperAdditivity, AntigenDemand, Capacity, RecoverAnnuity, ElasticityL, SegmentedPrices, ReservationPriceLim;
option solver knitro; 
option knitro_options 'ms_enable=0 mip_integral_gap_rel=1e-4';  
option display_eps 1.0e-09; 
option solution_round 6; 

# ABPFLTSS: Fully linearized problem formulation (MIP) -bundle allocation problem
problem ABPFLTSS: X, Y, W, g, TSS, SuperAdditivity, AntigenDemand, Capacity, ReservationPriceLim, RecoverAnnuityLMcMain, RecoverAnnuityLMc1, RecoverAnnuityLMc2,RecoverAnnuityLMc3, RecoverAnnuityLMc4, RecoverAnnuityLMc5, RecoverAnnuityLMc6, ElasticityL;
option solver cplex;
option cplex_options 'iisfind 1';
option presolve_eps 1.0e-06;

# ABPTCS: Total Customer Surplus -LP pricing problem - allows iis (aka. ATCS)
problem ABPTCS: Y,TCS, SuperAdditivity, AntigenDemand, Capacity, ElasticityNLY, ReservationPriceLim, RecoverAnnuity;
#problem ABPTCS: Y,TCS, loss, SuperAdditivity, AntigenDemand, Capacity, ElasticityNLY, SegmentedPrices, ReservationPriceLim, RecoverAnnuity, DefineLoss, GenericPrice2;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05';
option presolve 0;
option presolve_eps .01;

# ABPTPF: Total Profits - LP pricing problem - allows iis (aka ATPS)
problem ABPTPF: Y,TPF, SuperAdditivity, AntigenDemand, Capacity, ElasticityNLY, ReservationPriceLim, RecoverAnnuity;
#problem ABPTPF: Y,TPF, loss,SuperAdditivity, AntigenDemand, Capacity, ElasticityNLY, SegmentedPrices, ReservationPriceLim, RecoverAnnuity, DefineLoss, GenericPrice2;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05';
option presolve 0;
option presolve_eps .01;

# Problem definitions used for Zhang's heuristic 2012-paper
problem ABPTSS_LH_com:X,Y,g,y,TSSH_com,AntigenDemandH, CapacityH, RecoverAnnuityH,ElasticityLH,SegmentedPricesH,ReservationPriceLimH,PartialProcurementH,PartialSetUpH; 
option solver knitro;
option knitro_options 'presolve_dbg=0 ms_enable=1 ms_maxsolves=100 ma_terminate=2 mip_integral_gap_rel=1e-9 bar_feasible=3';  
option display_eps 1.0e-01; 
option solution_round 5;

# Problem definitions used for Zhang's heuristic 2012-paper
problem ABPTSS_LH_mono:X,Y,g,y,deficit,TSSH_mono,AntigenDemandH_mono, CapacityH, RecoverAnnuityH,ElasticityLH,SegmentedPricesH,ReservationPriceLimH,PartialProcurementH,PartialSetUpH; 
option solver knitro;
option knitro_options 'presolve_dbg=0 ms_enable=1 ms_maxsolves=101 ma_terminate=2 mip_integral_gap_rel=1e-9  bar_feasible=3';  
option display_eps 1.0e-01; 
option solution_round 5;

# Problem definitions used as the ABPTSS problem for multiple entities
problem ABPTSSME: X,Y, g, TSSME, SuperAdditivity, AntigenDemand, CapacityME, PriceElasticity, SegmentedPrices, ReservationPriceLimME, TargetedSale, GenericPrice, MaxBDose;
option display_eps .000000001; 
option solution_round 6;
option presolve_eps .01;
option solver cplex;
option cplex_options "iisfind=1";

# Problem definition for the TCS problem using multiple entities
problem ABPTCSME: Y, TCSME, SuperAdditivity2,  CapacityME, PriceElasticity, ReservationPriceLimME, GenericPrice2;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05';
option presolve_eps .01;
option send_statuses 0;

# Problem definition for the TPF problem using multiple entities
problem ABPTPFME: Y,TPFME, SuperAdditivity2,  CapacityME, PriceElasticity, ReservationPriceLimME, GenericPrice2;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05';
option presolve_eps .01;
option send_statuses 0;

#problem DistributionE: Z, SZ,Belonging, MinSZ, Distribution;
#problem DistributionE: Z, K, Distribution, ChooseMinDistribution, Ordering, Distribute, Distribute2;
#problem DistributionE: Z, K, Distribution, ChooseMinDistributionP, Ordering, Distribute, Distribute2;
#problem DistributionE: Z, K, s2, TargetPower, Distribute, Distribute2, Ordering, MinT2, Conditional1, Conditional2;
#problem DistributionE: Z, K, s2, TargetPower, Distribute, Distribute2, Ordering, MinT2P, Conditional1, Conditional2;
#problem DistributionE: Z, K, s2, TargetPower, Distribute, Distribute2, Ordering, MinT3, Conditional1, Conditional2;
#problem DistributionE: Z, K, s2, TargetPower, Distribute, Distribute2, Ordering, MinT3P, Conditional1, Conditional2;
#problem DistributionE: Z, K, s2, TargetPower, Distribute, Distribute2, Ordering, MinT4, Conditional3, Conditional4;
#problem DistributionE: Z, K, s2, TargetPower, Distribute, Distribute2, Ordering, MinT4P, Conditional3, Conditional4;
#problem DistributionE: Z, TargetPowerP, Distribute, Ordering;	#Population
#problem DistributionE: Z, TargetPower, Order1, Order2, Order3, Order4, Order5, Order6, Order7, Order8, Order9, Order10, Order11, Order12;
#problem DistributionE: Z, OrderA, OrderB, OrderC, OrderD;
#problem DistributionE: Z, OrderAll, OrderAll2;
#problem DistributionE: Z, TargetPowerP, Order1, Order2, Order3, Order4, Order5, Order6, Order7, Order8, Order9, Order10, Order11, Order12;
/*option solver knitro;
option knitro_options 'presolve_dbg=0 ms_enable=1 ms_maxsolves=101 ma_terminate=2 mip_integral_gap_rel=1e-9  bar_feasible=3'; */
#option solver ilogcp;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 0;
option presolve_eps .01;

problem ABPF: F, g, TSSForm, AntigenDemandF, CapacityF;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 0;
option presolve_eps .01;

problem ABPF2: X, O, TCSF, SupplyF, CostF, DemandF, RPF;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 0;
option presolve_eps .01;

problem ABPF3: X, O, TPFF, SupplyF, CostF, DemandF, RPF;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 0;
option presolve_eps .01;

problem ABPFM: X, g,TSSF3, ADF2, CapacityF2;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 1;
option presolve_eps .01;

problem ABPFM2: OM, TCSF2, CostF2, RPF2, Lazy2;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 0;
option presolve_eps .01;

problem ABPFE2: OM, TCSF2, CostF3, RPF2, Lazy2;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 0;
option presolve_eps .01;

problem ABPFM3: OM, TPFF2, CostF2, RPF2, Lazy, Lazy2, Justice;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 0;
option presolve_eps .01;

problem ABPFE3: OM, TPFF2, CostF3, RPF2, Lazy, Lazy2, Justice;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	#Check the relative gap. It should be in the magnitude of the objective 10^-4; -3, etc.*/
option presolve 0;
option presolve_eps .01;
