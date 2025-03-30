/*
 * A Baseline RBC Model
 *
 * Sune Grøn Pedersen
 * Aarhus University
 *
 */

//****************************************
// House keeping
//****************************************

close all

//****************************************
// Preamble block
//****************************************

% 1.1 Endogenous variables

var
c   ${c}$   (long_name='Forbrug')
l   ${l}$   (long_name='Arbejdsudbud')
k   ${k}$   (long_name='Kapital')
in  ${l}$   (long_name='Investering')
r   ${r}$   (long_name='Rente')
w   ${w}$   (long_name='Lønrate')
y   ${y}$   (long_name='Output')
A   ${A}$   (long_name='Teknologi')
pr  ${pr}$   (long_name='Produktivitet')
;

% 1.2 Exogenous variables

varexo
    eps
;

predetermined_variables k;

% 1.3 Parameters

parameters
ALPHA
BETA
DELTA
RHOA
SIGMA
PSI
;

ALPHA   = 0.41;
BETA    = 0.99;
DELTA   = 0.047;
RHOA    = 0.89;
PSI     = 2;

% 1.4 Steady State Values

parameters
A_ss
r_ss
kol
yol
inol
col
l0
l_ss
k_ss
y_ss
inol_ss
c_ss
in_ss
w_ss
;

A_ss        = 1;
r_ss        = BETA^(-1) - 1 + DELTA;

kol         = (r_ss/ALPHA)^(1/(ALPHA-1));
yol         = kol^(ALPHA);
inol        = DELTA*kol;
col         = yol - inol;

w_ss        = (1-ALPHA)*kol^(ALPHA);

l0          = 1/3;
l_ss        = baseline_l_solver(l0,col,w_ss,PSI);

k_ss        = kol*l_ss;
y_ss        = yol*l_ss;
inol_ss     = inol*l_ss;
c_ss        = col*l_ss;
in_ss        = inol*l_ss;

model;

1/c         = PSI/(w*(1-l));
1/c         = BETA*(1-DELTA+r(+1))/c(+1);
k(+1)       = in + (1-DELTA)*k;
r           = ALPHA*(k/(A*l))^(ALPHA-1);
w           = (1-ALPHA)*A^(1-ALPHA)*(k/l)^(ALPHA);
y           = k^(ALPHA)*(A*l)^(1-ALPHA);
in          = y-c;
pr          = y/l;
log(A)      = RHOA*log(A(-1)) + eps;

end;

steady_state_model;

A           = A_ss;
r           = r_ss;
l           = ((1-ALPHA)*((1/BETA)-1+DELTA))/(PSI*((1/BETA)-1+DELTA-DELTA*ALPHA)+(1-ALPHA)*((1/BETA)-1+DELTA));
y           = y_ss;
k           = k_ss;
w           = w_ss;
c           = c_ss;
in          = in_ss;
pr          = y_ss/l_ss;

end;
model_diagnostics;
steady;

//****************************************
// Shocks
//****************************************

shocks;
    var eps;
    stderr 0.006;
end;

check;
steady;

//****************************************
// Simulation and Statistics
//****************************************

stoch_simul(order=1,drop=30,periods=192,irf=150, loglinear,hp_filter=1600,tex);










