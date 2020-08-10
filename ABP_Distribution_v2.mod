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

#Perhaps maximize the minimum purchase power of target markets?

set E = 1..num_ent;	#Entities

param s{M,M};	#Similarity between markets (by geography, for example)

param chi;	#Minimum acceptable similarity for 2 markets in an entity.

param zeta;	#Maximum acceptable gni difference wihin an entiy.


var Z{M} integer in E;
var ZM{E} integer in M;	#Shows which market the entity belongs to
#var SZ;

maximize Distribution: 
	min{e in E}(sum{m in M, b in B}(if Z[m]==e then R[b,m] else 0));

#maximize TargetPower: #(obj2)
#	min{e in E}(if sum{m in M: gni_p[m]<=targetGNI}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M, b in B}(if Z[m]==e then R[b,m] else 0) else 1e10);

maximize DistributionP: 
	min{e in E}(sum{m in M}(if Z[m]==e then l[m] else 0));

#maximize TargetPowerP: #(obj2)
#	min{e in E}(if sum{m in M: gni_p[m]<=targetGNI}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M}(if Z[m]==e then l[m] else 0) else 1e10);


#Tries maximizing low-income while minimizing high-income (obj3)
maximize TargetPower:
	min{e in E}(if sum{m in M: gni_p[m]<=targetGNI}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M, b in B: gni_p[m] <= hmGNI}(if Z[m]==e then R[b,m] else 0) else 1e10);

#Tries maximizing low-income while minimizing high-income (obj3)
maximize TargetPowerP:
	min{e in E}(if sum{m in M: gni_p[m]<=targetGNI}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M: gni_p[m] <= hmGNI}(if Z[m]==e then l[m] else 0) else 1e10);

#Maximizes the purchasing power in HIC (obj4)
#maximize TargetPower:
#	min{e in E}(if sum{m in M: gni_p[m]>=hmGNI}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M, b in B}(if Z[m]==e then R[b,m] else 0) else 1e10);

#Maximizes the purchasing power in HIC (obj4)
#maximize TargetPowerP:
#	min{e in E}(if sum{m in M: gni_p[m]>=hmGNI}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M}(if Z[m]==e then l[m] else 0) else 1e10);


#Attempts to recreate the benchmark
#maximize TargetPower:
#	max{e in E}(if sum{m in M}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M, b in B}(if Z[m]==e then R[b,m] else 0) else 1e10) - min{e in E}(if sum{m in M}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M, b in B}(if Z[m]==e then R[b,m] else 0) else 1e10);

#Attempts to recreate the benchmark
#maximize TargetPowerP:
#	max{e in E}(if sum{m in M}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M}(if Z[m]==e then l[m] else 0) else 1e10) - min{e in E}(if sum{m in M}(if Z[m]==e then 1 else 0) >= 1 then sum{m in M}(if Z[m]==e then l[m] else 0) else 1e10);


#3 if benchmark
subject to Distribute{e in E}:
	sum{m in M} (if Z[m]==e then 1 else 0) >= 1;

subject to Ordering{e in 2..4}:
	sum{m in M, b in B} (if Z[m]==e then gni_p[m] else 0) <= sum{m in M, b in B} (if Z[m]==(e-1) then gni_p[m] else 0);

subject to Order1:
	Z[1] = 1;
subject to Order2:
	Z[2] = 1;
subject to Order3:
	Z[3] = 1;
subject to Order4:
	Z[4] = 2;
subject to Order5:
	Z[5] = 2;
subject to Order6:
	Z[6] = 2;
subject to Order7:
	Z[7] = 3;
subject to Order8:
	Z[8] = 3;
subject to Order9:
	Z[9] = 3;
subject to Order10:
	Z[10] = 4;
subject to Order11:
	Z[11] = 4;
subject to Order12:
	Z[12] = 4;
#subject to Summing{e in E}:
#	m

#subject to Belonging{m in M}:  
#	sum{e in E}Z[e,m] = 1;

#SZ <= sum{m in M, b in B}(if Z[m]==e:R[b,m]);

/*subject to Similarity{e in E, m1 in M, m2 in M}:
	W[e,m1]*s[m1,m2]+W[e,m2]*s[m1,m2] >= chi;

subject to GNIRes: 
	abs(W[e,m1]*gni_p[m1]-E[e,m2]*gni_p[m2])>=zeta;	#Should only count when X[b,m]>0*/

#kEEP A GAP BETWEEN THE ip AND THE INITIAL GAP. <- There might be a different.


