/*
 * A Baseline RBC Model with Calvo Price Rigidity, Money in the Utility Function 
 * and Indivisible Labor
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
y
w
r
M
P_tilde
P
P_star
PSI1
PSI2
A
gam
in
vartheta
prod
;


varexo
    epsA epsG
;

% 1.2 Parameters
parameters
ALPHA
BETA
gam_c
PSI
xi
DELTA
RHOa
SIGa
SIGg
pi
ETA
THETA
h


r_ss
w_ss
c_ss
A_ss
P_ss
P_star_ss
y_ss
k_ss
l_ss
vartheta_ss
in_ss
P_tilde_ss
M_ss
PSI1_ss
PSI2_ss
prod_ss
;


ALPHA   = 0.41;
BETA    = 0.99;
gam_c   = 1;
PSI     = 1.733;
xi      = 0.01;
DELTA   = 0.047;
RHOa    = 0.89;
SIGa    = 0.006;
SIGg    = 0.001;
pi      = 0.56;
ETA     = 11;
THETA   = 0.50;
h       = 0.5453;


% 1.3  solving for the rest of the steady states

A_ss        = 1;
gam_ss      = 1;
P_ss        = 1;
P_tilde_ss  = P_ss;
P_star_ss   = P_ss;

r_ss        = BETA^(-1) -1 + DELTA;
w_ss        = ((ETA-1)*(1-ALPHA)/ETA)^(1/(1-ALPHA)) * (ALPHA/((1-ALPHA)*r_ss))^(ALPHA/(1-ALPHA));
c_ss        = -w_ss*h*(PSI*log(1-h))^(-1);
M_ss        = c_ss*(1-BETA)^(-1)*xi*P_ss;

y_ss        = -w_ss*h/(PSI*log(1-h)*(w_ss*(r_ss*(1-ALPHA)/(w_ss*ALPHA))^(ALPHA) + 1/ETA - (DELTA-r_ss)*(r_ss*(1-ALPHA)/(w_ss*ALPHA))^(ALPHA-1)));

l_ss        = (r_ss*(1-ALPHA)/(w_ss*ALPHA))^(ALPHA)*y_ss;
k_ss        = (r_ss*(1-ALPHA)/(w_ss*ALPHA))^(ALPHA-1)*y_ss;

PSI1_ss     = P_ss^(ETA)*y_ss*(1-THETA*BETA)^(-1);
PSI2_ss     = w_ss*A_ss^(ALPHA-1)/(1-ALPHA)*(r_ss*(1-ALPHA)/(ALPHA*w_ss))^(ALPHA)*P_ss^(ETA+1)*y_ss*(1-THETA*BETA)^(-1);
vartheta_ss = y_ss/ETA;
in_ss       = DELTA*k_ss;
prod_ss     = y_ss/l_ss;

//****************************************
// Model block
//****************************************

model;

1/c                         = -w^(-1)*h^(-1)*PSI*log(1-h);
1/c                         = BETA*(1-DELTA+r(+1))*(c(+1))^(-1);
1/c                         = xi*P/M+BETA*P*(c(+1)*P(+1))^(-1);
c+k-(1-DELTA)*k(-1) + M/P   = w*l + r*k(-1) + vartheta + M(-1)/P + (gam-1)*M(-1)/P;
P_star                      = (1-THETA)*P_tilde^(-ETA)+THETA*(P(-1)/P)^(-ETA)*P_star(-1);
y/P_star                    = k(-1)^(ALPHA)*(A*l)^(1-ALPHA);
1                           = (1-THETA)*(P_tilde/P)^(1-ETA)+THETA*(P(-1)/P)^(1-ETA);
P_tilde                     = (ETA/(ETA-1))*(PSI2/PSI1);
PSI1                        = P^(ETA)*y+THETA*BETA*PSI1(+1);
PSI2                        = w*A^(ALPHA-1)/(1-ALPHA)*(r*(1-ALPHA)/(ALPHA*w))^(ALPHA)*P^(ETA+1)*y+THETA*BETA*PSI2(+1);
%l                           = A^(ALPHA-1)*(r*(1-ALPHA)/(ALPHA*w))^(ALPHA)*y;
%k                           = A^(ALPHA-1)*(r*(1-ALPHA)/(ALPHA*w))^(ALPHA-1)*y;
%r                           = w*A^(ALPHA-1)/(1-ALPHA)*(r*(1-ALPHA)/(ALPHA*w))^(ALPHA)*ALPHA*k(-1)^(ALPHA-1)*(A*l)^(1-ALPHA);
%w                           = w*A^(ALPHA-1)/(1-ALPHA)*(r*(1-ALPHA)/(ALPHA*w))^(ALPHA)*(1-ALPHA)*k(-1)^(ALPHA)*A^(1-ALPHA)*l^(-ALPHA);
l/k(-1)                     = (1-ALPHA)*r/(ALPHA*w);
y                           = c+in;
prod                        = y/l;
M                           = gam*M(-1);
vartheta                    = y - w*A^(ALPHA-1)/(1-ALPHA)*(r*(1-ALPHA)/(ALPHA*w))^(ALPHA)*y/P_star;
log(A)                      = RHOa*log(A(-1)) + epsA;
log(gam)                    = (1-pi)*log(gam_c)+pi*log(gam(-1))+epsG;

%MC  = w*A^(ALPHA-1)/(1-ALPHA)*(r*(1-ALPHA)/(ALPHA*w))^(ALPHA);

end;
model_diagnostics;
steady_state_model;

A           = A_ss;
P           = P_ss;
P_tilde     = P_tilde_ss;
P_star      = P_star_ss;
gam         = 1;
r           = r_ss;
w           = w_ss;
y           = y_ss;
c           = c_ss;
M           = M_ss;
vartheta    = vartheta_ss;
l           = l_ss;
k           = k_ss;
PSI1        = PSI1_ss;
PSI2        = PSI2_ss;
in          = in_ss;
prod        = prod_ss;

end;
steady;


shocks;
    var epsA;
    stderr 0.0053;
    var epsG;
    stderr 0.0007;
end;

write_latex_static_model;
write_latex_dynamic_model;

check;
model_diagnostics;
steady;
model_diagnostics;

//****************************************
// Simulation and Statistics
//****************************************  

stoch_simul(order=1,drop=30,periods=192,irf=150, loglinear,hp_filter=1600,tex);

