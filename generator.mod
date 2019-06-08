param j;

set FORM := 1..j;
set TA;
set MAN;

param f{FORM, B};
param MAX := sum{b in B}1.01^b;

var aux{FORM} binary;
var formul{B} integer >= 0;


maximize Garbage: sum{b in B} formul[b]*1.01^b;

s.t. Compliance: sum{b in B} formul[b] >= 1;

s.t. Enforce {m in 1..n_markets, b in B}: formul[b] <= D[b,m];

s.t. Difference {i in FORM}: sum{b in B}(formul[b]-f[i,b])^2 >= .1;

#s.t. Difference {i in FORM}: sum{b in B}(g[b]-f[i,b]) * (1.01^b) >= aux[i] - MAX*(1-aux[i]);
#s.t. Difference2 {i in FORM}: sum{b in B}(g[b]-f[i,b]) * (1.01^b) <= MAX*aux[i] -(1-aux[i]);
#s.t. Difference {i in FORM}: g not in f; 

s.t. TargetMan {b in B: k[b] not in MAN}: formul[b] =0;

#problem Form: g, aux, Garbage, Compliance, Enforce, Difference, Difference2, TargetMan;
problem Form: formul, aux, Garbage, Compliance, Enforce, Difference, TargetMan;
#option solver ilogcp;
/*option solver cplex;
option cplex_options 'primal iisfind=1 feasibility=1e-05 mipgap=1e-3';	
option presolve 0;
option presolve_eps .01;*/

#option solver gurobi;
#option gurobi_options 'iisfind=1';

option solver knitro;