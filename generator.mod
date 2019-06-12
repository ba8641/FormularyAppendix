param j;

set FORM := 1..j;
set DOSE := 1..3;
set TA;
set MAN;

param f{FORM, B, DOSE};
param MAX := sum{b in B}1.1^b;

var aux{FORM} binary;
var formul{B, DOSE} binary;

var XF{FORM,M} >= 0;		#Quantity of formulary going to market m;
var YF{FORM,M} >= 0;
var gF{FORM} binary;

param RF{FORM,M};	#Reservation price of form in market.
param CF{FORM};


maximize Garbage: sum{i in FORM, b in B, do in DOSE}(formul[b,do]-f[i,b,do])^2;

maximize TSSFormul: 
sum{b in FORM, m in M} RF[b,m]*XF[b,m]- sum{b in B}C[b]*g[b] + 
(sum{a in A, m in M}(sum{form in FORM, b in B1[a], do in DOSE}XF[form,m]*do*f[form,b,do] - d[a,m]*l[m]));

s.t. Compliance: sum{b in B, do in DOSE} formul[b,do] >= 1;

#s.t. Enforce {m in 1..n_markets, b in B: D[b,m]}: formul[b] <= D[b,m];

#s.t. Difference {i in FORM}: sum{b in B}(formul[b]-f[i,b])^2 >= 1;

#s.t. Difference {i in FORM}: sum{b in B}(formul[b]-f[i,b]) * (1.1^b) >= aux[i] - MAX*(1-aux[i]);
#s.t. Difference2 {i in FORM}: sum{b in B}(formul[b]-f[i,b]) * (1.1^b) <= MAX*aux[i] -(1-aux[i]);
#s.t. Difference {i in FORM}: g not in f; 

s.t. TargetMan {b in B, do in DOSE: k[b] not in MAN}: formul[b,do] =0;

subject to ADF{a in A, m in M}:
 sum{form in FORM, b in B1[a], do in DOSE}XF[form,m]*f[form,b,do]*do <= d[a,m]*l[m];

subject to CapF{b in B}:
 sum{m in M, i in FORM, do in DOSE} XF[i,m]*f[i,b,do]*do <= S[b]*g[b];

#subject to RAF:
# sum{i in FORM, m in M} RF[i,m]*XF[i,m] >= sum{i in FORM, b in B: f[i,b] > 0}C[b]*gF[i];
 
subject to Buy{b in FORM}:
 sum{m in M}XF[b,m] <= (sum{i in B}S[i])*gF[b];
 
subject to RAFMan{p in P}:
 sum{m in M, i in FORM: sum{b in B2[p], do in DOSE}f[i,b,do] > 0} RF[i,m]*XF[i,m] >= sum{b in B2[p]}C[b]*g[b];

subject to Buy2{b in B, i in FORM: sum{do in DOSE}f[i,b,do] >0}:
 g[b] >= gF[i];
 
subject to MBDF{b in B, m in M}:
	sum{i in FORM, do in DOSE}XF[i,m]*f[i,b,do]*do <= D[b,m]*l[m];

#problem Form: g, aux, Garbage, Compliance, Enforce, Difference, Difference2, TargetMan;
problem Form: formul, aux, Garbage, Compliance, TargetMan;
#option solver ilogcp;
/*option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	
option presolve 0;
option presolve_eps .01;*/

option solver gurobi;
option gurobi_options 'iisfind=1';

#option solver gurobi;

problem TSSFORM: TSSFormul, XF, gF, g, CapF, RAFMan, Buy, Buy2, MBDF;
option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	
option presolve 0;
option presolve_eps .01;
