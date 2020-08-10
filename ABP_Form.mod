/*
	ABPTSSME, TCSME, and TPFME MODEL ELEMENTS (model file) and BASE FORMULATION
	-  based on Proano et al, Omega 40(2012):54-63
	- modular implementation
	- adapted from ABP_base.mod, by Dr.Proano.
	by Bruno Alves Maciel
	ba8641@rit.edu
	ROCHESTER INSTITUTE OF TECHNOLOGY
	
	NOTE:
	
	07.11.2017
	version I.01

	See additional comments at the end of the file
*/


/* General elements of the ABPME formulation */ \

param n_antigens >=0;
param n_bundles >=0;
param n_producers >=0;

set B:= 1..n_bundles; 
set A;
set M;
set P:= 1..n_producers;

set F within B;		#Defines formularies (from different manufacturers) that can be bought/sold.

set T{B} within M;	#Defines the target market

set A1{B} within A;
set B1{A} within B;
set B2{P} within B;
set F1{A} within F;


/* Elements used in SupperAdditivity constraint */
set Q{B} default {}; 
set QQ:= union {b in B} Q[b];
set N{QQ} within B;  

param R{F,M}; 
param l{M};
param C{B};
param CF{F};
param d{A,M};
param D{B,M};
param S{B};
param SF{F};
param k{B} within P;		#Isn't this basically B2?
param phi; 
param name{B} symbolic; 
param gni_p{M};
param theta;

param alpha;	#How much we penalize not buying the entire supply of a bundle

param minPrice{F};

param X2{f in F, m in 1..194} default 0;

param Cl{F}; 	#Keeps track of how much of the annuity costs remain to be paid.

#Parameters with the weigths of each part of the objective function
param WOri;	#Weight of the original portion of the TSS, which tries to maximize vaccine allocations taking care not to raise costs too much.
param WSat;	# Weight of the new portion of the model, which sees how much of the demand is satisfied.

param targetGNI default 1036;
param lmGNI default 4085;
param hmGNI default 12615;

param Xlow{b in B, m in M}>=0 default 0;
param Xup{b in B,m in M} default min(S[b], D[b,m]*l[m]);
param Ylow{b in B, m in M} >=0 default 0;
param Yup{b in B, m in M} default theta*R[b,m];

var K{F,P} >= 0 and <= 1 default 0;	#Fraction of the revenue of a formulary going to each manufacturer.

var X{F,M} >=0 default 0;
var Y{F,M}>=0 default 0; 
var g{F} binary default 0;

maximize TSSForm: 
WOri*(sum{b in F, m in T[b]} R[b,m]*X[b,m]- sum{b in F}CF[b]*sum{m in T[b]}X[b,m]/(SF[b]+1e-10)) - WSat*(sum{a in A, m in M: gni_p[m]<targetGNI}(d[a,m]*l[m] - sum{b in F1[a]}X[b,m]))- sum{b in B}g[b]; #Tweight-TS

maximize TCSForm: sum{b in F, m in T[b]} (R[b,m]- Y[b,m])*X[b,m];
#maximize TCSME: sum{b in B, m in T[b]: gni_p[m] < targetGNI} (R[b,m]- Y[b,m])*X[b,m] + .01*sum{b in B, m in T[b]: gni_p[m] >= targetGNI} (R[b,m]- Y[b,m])*X[b,m];	#TWeight-TC

maximize TPFME: sum{b in F, m in T[b]}Y[b,m]*X[b,m] - sum{b in F} CF[b]*g[b]*sum{m in T[b]}X[b,m]/(SF[b]+1e-10);

subject to AntigenDemand{a in A, m in M}: sum{b in F1[a]}X[b,m] <= d[a,m]*l[m];

subject to Capacity{b in F}: sum{m in M} X[b,m] <= SF[b]*g[b];

subject to RecoverAnnuityTotal{b in F}: sum{m in M} Y[b,m]*X[b,m] >= CF[b]*g[b];

#NON-LINEAR!!!
subject to RecoverAnnuityManufacturer{p in P}: sum{b in F, m in M} Y[b,m]*X[b,m]*K[b,p] >= sum{t in B: t in b and k[t]==p}C[t]*g[b];	    

subject to ReservationPriceLim{b in B, m in M}: 
	Y[b,m] <= theta*R[b,m]*g[b];

subject to CapacityME{b in F}: sum{m in T[b]} X[b,m] <= S[b]*g[b];

subject to GenericPrice{b in F, m in T[b]}: Y[b,m] >= g[b]*minPrice[b];	#Should only count when X[b,m]>0

subject to GenericPrice2{b in F, m in T[b]: X2[b,m]>0}: Y[b,m] >= g[b]*minPrice[b];	#Should only count when X[b,m]>0

subject to PriceElasticity{b in F, m in T[b]:SF[b]>0}:
	Y[b,m] >= ((SF[b]*g[b]-X[b,m])*alpha+SF[b]*g[b])/SF[b] * (C[b]/SF[b]);		#Look for better elasticity relationships?	min when X = supply; 

subject to ReservationPriceLimME{b in F, m in T[b]}: 
	Y[b,m] <= theta*R[b,m]*g[b];

subject to TargetedSale{b in F}:
	sum{m in T[b]}X[b,m] >= sum{m in M}X[b,m];

subject to MaxBDose{b in F, m in T[b]}:
	X[b,m] <= D[b,m]*l[m];
