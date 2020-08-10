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


param RA{B};			#Average RP weighted by population.
param dA{A};			# Total demand for a given antigen.

var F{B};			#Quantity of B in F.

maximize TSSForm: 
sum{b in B}(F[b]*RA[b] - C[b]*g[b]); #Tweight-TS

subject to AntigenDemandF{a in A}: sum{b in F1[a]}F[b] >= dA[a];

subject to CapacityF{b in B}: F[b] <= S[b]*g[b];
