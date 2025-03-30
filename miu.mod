/*
 * A Baseline RBC Model with Money in the Utility Function
 *
 * Sune Grøn Pedersen
 * Aarhus University
 *
 */

//****************************************
// 0. House keeping
//****************************************

close all

//****************************************
// 1. Preamble block
//****************************************

% 1.1 Endogenous variables

var

c
l
k
M
P
r
w
y
A
gam
in
prod

;

predetermined_variables k;

varexo
    epsA epsG
;

% 1.2 Parameters
parameters
ALPHA
BETA
gam_c
PSI
XI
DELTA
RHOa
SIGa
SIGg
pi
;

ALPHA   = 0.41;
BETA    = 0.99;
gam_c   = 1;
PSI     = 2;
XI      = 0.01;
DELTA   = 0.047;
RHOa    = 0.89;
SIGa    = 0.006;
SIGg    = 0.006;
pi      = 0.56;

% 1.3 SS

parameters
c_ss
l_ss
k_ss
M_ss
P_ss
r_ss
w_ss
y_ss
A_ss
gam_ss

yol
col
kol

;

A_ss        = 1;
gam_ss      = 1;
P_ss        = 1;

r_ss        = BETA^(-1) - 1 + DELTA;
w_ss        = (1-ALPHA)*(r_ss/ALPHA)^(ALPHA/(ALPHA-1));

kol     = (r_ss/ALPHA)^(1/(ALPHA-1));
yol     = (kol)^(ALPHA);
col     = yol - DELTA*kol;

l0      = 1/3;
l_ss    = miu_l_solver(l0,col,w_ss,PSI);

k_ss    = kol*l_ss;
y_ss    = yol*l_ss;
c_ss    = col*l_ss;

M_ss    = P_ss*XI*gam_ss*c_ss/(gam_ss-BETA);

//****************************************
// Model block
//****************************************

model;

1/c                         = PSI*(w*(1-l))^(-1);
1/c                         = BETA*(1-DELTA+r(+1))*(c(+1))^(-1);
1/c                         = XI*P/M+BETA*P*(c(+1)*P(+1))^(-1);
c+k(+1)-(1-DELTA)*k + M/P   = w*l + r*k + M(-1)/P + (gam-1)*M(-1)/P;
w                           = (1-ALPHA)*y/l;
r                           = ALPHA*y/k;
M                           = gam*M(-1);
y                           = k^(ALPHA)*(A*l)^(1-ALPHA);
log(A)                      = RHOa*log(A(-1)) + epsA;
log(gam)                    = (1-pi)*log(gam_c)+pi*log(gam(-1))+epsG;
in                          = y-c;
prod                        = y/l;


end;

steady_state_model;

A   = 1;
gam = 1;
P   = 1;
c   = c_ss;
l   = l_ss;
y   = y_ss;
k   = k_ss;
r   = BETA^(-1) - 1 + DELTA;
w   = (1-ALPHA)*(r_ss/ALPHA)^(ALPHA/(ALPHA-1));
M   = P*XI*gam*c/(gam-BETA);
in  = DELTA*k;
prod    = y/l;

end;
model_diagnostics;
steady;

shocks;
    var epsA;
    stderr 0.0053;
    var epsG;
    stderr 0.0007;
end;

check;
steady;

//****************************************
// Simulation and Statistics
//****************************************  

stoch_simul(order=1,drop=30,periods=192,irf=150, loglinear,hp_filter=1600,tex);

