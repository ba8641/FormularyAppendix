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

#set E = 1..num_ent;	#Entities

param s{M,M};	#Similarity between markets (by geography, for example)

param chi;	#Minimum acceptable similarity for 2 markets in an entity.

param zeta;	#Maximum acceptable gni difference wihin an entiy.


var Z{E,M} binary;
var K;
var ZM{E} integer in M;	#Shows which market the entity belongs to
var s2{E} binary;
#var SZ;

maximize Distribution: 
	K;

maximize TargetPower: #(obj2)
	K;

maximize DistributionP: 
	min{e in E}(sum{m in M}Z[e,m]*l[m]);

#maximize TargetPowerP: #(obj2)
#	min{e in E}(if sum{m in M: gni_p[m]<=targetGNI}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M}(if Z[e,m]==e then l[m] else 0) else 1e10);


#Tries maximizing low-income while minimizing high-income (obj3)
#maximize TargetPower:
#	min{e in E}(if sum{m in M: gni_p[m]<=targetGNI}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M, b in B: gni_p[m] <= hmGNI}(if Z[e,m]==e then R[b,m] else 0) else 1e10);

#Tries maximizing low-income while minimizing high-income (obj3)
maximize TargetPowerP:
	min{e in E}(if sum{m in M: gni_p[m]<=targetGNI}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M: gni_p[m] <= hmGNI}(if Z[e,m]==e then l[m] else 0) else 1e10);

#Maximizes the purchasing power in HIC (obj4)
#maximize TargetPower:
#	min{e in E}(if sum{m in M: gni_p[m]>=hmGNI}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M, b in B}(if Z[e,m]==e then R[b,m] else 0) else 1e10);

#Maximizes the purchasing power in HIC (obj4)
#maximize TargetPowerP:
#	min{e in E}(if sum{m in M: gni_p[m]>=hmGNI}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M}(if Z[e,m]==e then l[m] else 0) else 1e10);


#Attempts to recreate the benchmark
#maximize TargetPower:
#	max{e in E}(if sum{m in M}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M, b in B}(if Z[e,m]==e then R[b,m] else 0) else 1e10) - min{e in E}(if sum{m in M}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M, b in B}(if Z[e,m]==e then R[b,m] else 0) else 1e10);

#Attempts to recreate the benchmark
#maximize TargetPowerP:
#	max{e in E}(if sum{m in M}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M}(if Z[e,m]==e then l[m] else 0) else 1e10) - min{e in E}(if sum{m in M}(if Z[e,m]==e then 1 else 0) >= 1 then sum{m in M}(if Z[e,m]==e then l[m] else 0) else 1e10);

s.t. ChooseMinDistribution{e in E}:
	K <= (sum{m in M, b in B}Z[e,m]*R[b,m]);

s.t. ChooseMinDistributionP{e in E}:
	K <= (sum{m in M}Z[e,m]*l[m]);

s.t. MinT2{e in E}:
	K <= s2[e]* (sum{b in B, w in M}Z[e,w]*R[b,w]) + (1-s2[e])*1e10;

s.t. MinT2P{e in E}:
	K <= s2[e]* (sum{w in M}Z[e,w]*l[w]) + (1-s2[e])*1e10;
	
s.t. MinT3{e in E}:
	K <= s2[e]* (sum{b in B, w in M: gni_p[w] <= hmGNI}Z[e,w]*R[b,w]) + (1-s2[e])*1e10;

s.t. MinT3P{e in E}:
	K <= s2[e]* (sum{w in M: gni_p[w] <= hmGNI}Z[e,w]*l[w]) + (1-s2[e])*1e10;

s.t. MinT4{e in E}:
	K <= s2[e]* (sum{b in B, w in M}Z[e,w]*R[b,w]) + (1-s2[e])*1e10;

s.t. MinT4P{e in E}:
	K <= s2[e]* (sum{w in M}Z[e,w]*l[w]) + (1-s2[e])*1e10;

	
s.t. Conditional1{e in E}:
	sum{m in M: gni_p[m]<=targetGNI}Z[e,m] >= 2 - 192*(1-s2[e]);
	
s.t. Conditional2{e in E}:
	1 >= sum{m in M: gni_p[m]<=targetGNI}Z[e,m] + 1 - 192*s2[e];
	
s.t. Conditional3{e in E}:
	sum{m in M: gni_p[m]>=hmGNI}Z[e,m] >= 2 - 192*(1-s2[e]);
	
s.t. Conditional4{e in E}:
	1 >= sum{m in M: gni_p[m]>=hmGNI}Z[e,m] + 1 - 192*s2[e];

#3 if benchmark
subject to Distribute{e in E}:
	sum{m in M} Z[e,m] >= 1;
	
subject to Distribute2{m in M}:
	sum{e in E} Z[e,m] =1;

subject to Ordering{e in 2..194}:
	sum{m in M, b in B} (Z[e,m]*gni_p[m]) <= sum{m in M, b in B} (Z[e-1,m]*gni_p[m]);

subject to OrderA{m in 1..48}:
	Z[1,m] = 1;
	
subject to OrderB{m in 49..97}:
	Z[2,m] = 1;
	
subject to OrderC{m in 98..146}:
	Z[3,m] = 1;
	
subject to OrderD{m in 147..194}:
	Z[4,m] = 1;

subject to OrderAll{e in E, m in M: e ==m}:
	Z[e,m]=1;
	
subject to OrderAll2{e in E, m in M: e <>m}:
	Z[e,m]=0;

subject to Order1:
	Z[1,1] = 1;
subject to Order2:
	Z[1,2] = 1;
subject to Order3:
	Z[1,3] = 1;
subject to Order4:
	Z[2,4] = 1;
subject to Order5:
	Z[2,5] = 1;
subject to Order6:
	Z[2,6] = 1;
subject to Order7:
	Z[3,7] = 1;
subject to Order8:
	Z[3,8] = 1;
subject to Order9:
	Z[3,9] = 1;
subject to Order10:
	Z[4,10] = 1;
subject to Order11:
	Z[4,11] = 1;
subject to Order12:
	Z[4,12] = 1;
#subject to Summing{e in E}:
#	m

#subject to Belonging{m in M}:  
#	sum{e in E}Z[e,m] = 1;

#SZ <= sum{m in M, b in B}(if Z[e,m]==e:R[b,m]);

/*subject to Similarity{e in E, m1 in M, m2 in M}:
	W[e,m1]*s[m1,m2]+W[e,m2]*s[m1,m2] >= chi;

subject to GNIRes: 
	abs(W[e,m1]*gni_p[m1]-E[e,m2]*gni_p[m2])>=zeta;	#Should only count when X[b,m]>0*/

#kEEP A GAP BETWEEN THE ip AND THE INITIAL GAP. <- There might be a different.


