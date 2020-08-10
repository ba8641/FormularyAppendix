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


/* General elements of the ABPME formulation */ 

#param n_entities;

set T{B} within M;	#Defines the target market
set E = 1..n_entities;

param alpha;	#How much we penalize not buying the entire supply of a bundle

param minPrice{B};

param X2{b in B, m in 1..194} default 0;

param Cl{B}; 	#Keeps track of how much of the annuity costs remain to be paid.

#Parameters with the weigths of each part of the objective function
param WOri;	#Weight of the original portion of the TSS, which tries to maximize vaccine allocations taking care not to raise costs too much.
param WSat;	# Weight of the new portion of the model, which sees how much of the demand is satisfied.

param targetGNI default 1036;
param lmGNI default 4085;
param hmGNI default 12615;

param RA{B};			#Average RP weighted by population.
param dA{A};			# Total demand for a given antigen.

param ga{B} binary default 0;
var sigma >= 1 default 1;

param tolerance default 1.05;		#Wiggle room given from the calculated bundle value

var F{B} >= 0;			#Quantity of B in F.
var FM{B,P} >= 0;
var FE{E,B} >= 0;

var O{M} >= 0;

var OM{M,P} >= 0;

var MoneyLack{P} >= 0;



maximize TSSME: 
#WOri*(sum{b in B, m in T[b]} R[b,m]*X[b,m]- sum{b in B}C[b]*sum{m in T[b]}X[b,m]/(S[b]+1e-10)) - WSat*(sum{a in A, m in M}(d[a,m]*l[m] - sum{b in B1[a]}X[b,m]))- sum{b in B}g[b]; #Added demand
WOri*(sum{b in B, m in T[b]} R[b,m]*X[b,m]- sum{b in B}C[b]*sum{m in T[b]}X[b,m]/(S[b]+1e-10)) - WSat*(sum{a in A, m in M: gni_p[m]<targetGNI}(d[a,m]*l[m] - sum{b in B1[a]}X[b,m]))- sum{b in B}g[b]; #Tweight-TS
#WOri*(sum{b in B, m in T[b]} R[b,m]*X[b,m]- sum{b in B}C[b]*sum{m in T[b]}X[b,m]/(S[b]+1e-10)) - WSat*100*(sum{a in A, m in M: gni_p[m]<targetGNI}(d[a,m]*l[m] - sum{b in B1[a]}X[b,m])) - WSat*(sum{a in A, m in M: gni_p[m]>targetGNI}(d[a,m]*l[m] - sum{b in B1[a]}X[b,m]))- sum{b in B}g[b]; #Tweight-TSG

maximize TCSME: sum{b in B, m in T[b]} (R[b,m]- Y[b,m])*X[b,m];
#maximize TCSME: sum{b in B, m in T[b]: gni_p[m] < targetGNI} (R[b,m]- Y[b,m])*X[b,m] + .01*sum{b in B, m in T[b]: gni_p[m] >= targetGNI} (R[b,m]- Y[b,m])*X[b,m];	#TWeight-TC

maximize TPFME: sum{b in B, m in T[b]}Y[b,m]*X[b,m] - sum{b in B} C[b]*g[b]*sum{m in T[b]}X[b,m]/(S[b]+1e-10);


maximize TSSForm: 
sum{b in B}(F[b]*RA[b] - C[b]*g[b]); #Tweight-TS

maximize TSSF2: 
sum{p in P}sum{b in B2[p]}(FM[b,p]*RA[b] - C[b]*g[b]); #Tweight-TS

maximize TSSF3: 
sum{m in M, b in B}(X[b,m]*R[b,m]) - sum{b in B}(C[b]*g[b]);

maximize TCSF:
sum{m in M}(sum{b in B}X[b,m]*R[b,m] - O[m]);

maximize TCSF2:
sum{m in M}(sum{b in B}X[b,m]*R[b,m] - sum{p in P}OM[m,p]);

maximize TPFF:
sum{m in M}(O[m]) - sum{b in B}(C[b]*g[b]);

maximize TPFF2:
sum{m in M, p in P}(OM[m,p]) - (sum{b in B}C[b]*g[b]);

subject to AntigenDemandF{a in A}: sum{b in B1[a]}F[b] = sigma*dA[a];


subject to ADF2{a in A, m in M}: sum{b in B1[a]}X[b,m] == d[a,m]*l[m];

subject to CapacityF{b in B}: F[b] <= S[b]*g[b];

subject to CapacityF2{b in B}: sum{m in M}X[b,m] <= S[b]*g[b];

#subject to Producer{p in P}: sum{b not in B2[p]}F[b,p] := 0;

subject to SupplyF{b in B}: sum{m in M}X[b,m] <= S[b]*g[b];

subject to SupplyF2{b in B}: sum{m in M}X[b,m] <= F[b];

subject to DemandF{m in M, a in A}: sum{b in B1[a]}X[b,m] >= d[a,m]*l[m];

subject to CostF: sum{m in M}O[m] >= sum{b in B}C[b]*g[b];

subject to CostF2{p in P}: sum{m in M}OM[m,p] >= sum{b in B2[p]}C[b]*g[b];

subject to CostF3{p in P}: sum{m in M}OM[m,p] >= sum{b in B2[p]}(C[b]*g[b]/(S[b]+1e-10) * ((S[b]*g[b] - sum{m in M}X[b,m])*alpha+S[b]*g[b])/(S[b]+1e-10));

subject to RPF {m in M}: O[m] <= sum{b in B} R[b,m]*X[b,m];

subject to RPF2 {m in M}: sum{p in P}OM[m,p] <= sum{b in B} R[b,m]*X[b,m];

subject to Justice{p in P}: sum{m in M} OM[m,p] <= (sum{b in B2[p]}C[b]*ga[b]*tolerance) / (sum{b in B}C[b]*ga[b]+1e-10) * sum{pl in P, m in M}OM[m,pl];

subject to Lazy{p in P: sum{b in B2[p]}ga[b]==0}: sum{m in M}OM[m,p] <= 0;

subject to Lazy2{p in P}: sum{m in M}OM[m,p] >= sum{m in M, b in B2[p]}X[b,m]*minPrice[b];

subject to SuperAdditivity2{b in B, m in T[b], q in Q[b]: card(Q[b])>=1}:  
	Y[b,m] + (1-g[b])*sum{t in N[q]}R[t,m] >= sum{t in N[q]}(Y[t,m]);

subject to CapacityME{b in B}: sum{m in T[b]} X[b,m] <= S[b]*g[b];

subject to GenericPrice{b in B, m in T[b]}: Y[b,m] >= g[b]*minPrice[b];	#Should only count when X[b,m]>0

subject to GenericPrice2{b in B, m in T[b]: X2[b,m]>0}: Y[b,m] >= g[b]*minPrice[b];	#Should only count when X[b,m]>0

subject to PriceElasticity{b in B, m in M:S[b]>0}:
	Y[b,m] >= ((S[b]*g[b]-X2[b,m])*alpha+S[b]*g[b])/S[b] * (C[b]/S[b]);		#Look for better elasticity relationships?	min when X = supply; 

subject to ReservationPriceLimME{b in B, m in T[b]}: 
	Y[b,m] <= theta*R[b,m]*g[b];

subject to TargetedSale{b in B}:
	sum{m in T[b]}X[b,m] >= sum{m in M}X[b,m];

subject to MaxBDose{b in B, m in M}:
	X[b,m] <= D[b,m]*l[m];
