/*
 * A Baseline RBC Model with Money in the Utility Function and Indivisible Labor
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
h
;

ALPHA   = 0.41;
BETA    = 0.99;
gam_c   = 1;
PSI     = 2.14;
XI      = 0.01;
DELTA   = 0.047;
RHOa    = 0.89;
SIGa    = 0.006;
SIGg    = 0.006;
pi      = 0.56;
h0      = 0.5;
h       = h_solver(h0,ALPHA,BETA, PSI, DELTA);

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

l0      = 0.3;
l_ss = miu_indi_l_solver(l0,col,w_ss,PSI,h);

k_ss    = kol*l_ss;
y_ss    = yol*l_ss;

c_ss    = -w_ss*h*(PSI*log(1-h))^(-1);

M_ss    = P_ss*XI*gam_ss*c_ss/(gam_ss-BETA);

//****************************************
// Model block
//****************************************

model;

1/c                         = -PSI*log(1-h)*(w*h)^(-1);
1/c                         = BETA*(1-DELTA+r(+1))*(c(+1))^(-1);
1/c                         = XI*P/M+BETA*P*(c(+1)*P(+1))^(-1);
c+in + M/P                  = w*l + r*k + M(-1)/P + (gam-1)*M(-1)/P;
w                           = (1-ALPHA)*y/l;
r                           = ALPHA*y/k;
M                           = gam*M(-1);
y                           = k^(ALPHA)*(A*l)^(1-ALPHA);
log(A)                      = RHOa*log(A(-1)) + epsA;
log(gam)                    = (1-pi)*log(gam_c)+pi*log(gam(-1))+epsG;
in                          = k(+1) - (1-DELTA)*k;
prod                        = y/l;


end;

initval;

A   = 1;
gam = 1;
P   = 1;

l   = l_ss;
y   = y_ss;
k   = k_ss;
r   = BETA^(-1) - 1 + DELTA;
w   = (1-ALPHA)*(r_ss/ALPHA)^(ALPHA/(ALPHA-1));
c   = c_ss;
M   = P*XI*gam*c/(gam-BETA);
in  = DELTA*k;
prod    = y/l;

end;
model_diagnostics;
steady(tolf = 5e-2);

shocks;
    var epsA;
    stderr 0.0053;
    var epsG;
    stderr 0.0007;
end;

check;
steady(tolf = 5e-2);

//****************************************
// Simulation and Statistics
//**************************************** 

stoch_simul(order=1,drop=30,periods=192,irf=150, loglinear,hp_filter=1600,tex);
