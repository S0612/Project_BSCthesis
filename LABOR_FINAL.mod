/*
 * A Baseline RBC Model with Indivisible Labor
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
    c l k r w y in A prod
;
predetermined_variables k;

% 1.3 Exogenous variables

varexo
    epsa
;

% 1.3 Parameters
parameters
    ALPHA BETA PSI DELTA RHOa SIGMA h
;

parameters
    r_ss kol_ss w_ss yol_ss iol_ss col_ss l_ss k_ss y_ss i_ss c_ss;

% 1.4 Parameter Values

ALPHA       = 0.41;
RHOa        = 0.89;     
BETA        = 0.99;
DELTA       = 0.047;
PSI         = 2;
SIGMA       = 0.006;
h0          = 0.5;
h           = h_solver(h0,ALPHA,BETA, PSI, DELTA);

% 1.5 Steady State Computation

r_ss    = BETA^(-1)+DELTA-1;

kol_ss  = (r_ss/ALPHA)^(1/(ALPHA-1));
yol_ss  = kol_ss^(ALPHA);
col_ss  = yol_ss - DELTA*kol_ss;
iol_ss  = yol_ss - col_ss;

w_ss    = (1-ALPHA)*yol_ss;

l_ss    = -(1-ALPHA)*(BETA^(-1)-1+DELTA)*((PSI/h)*log(1-h)*(BETA^(-1)-1+DELTA-DELTA*ALPHA))^(-1);

k_ss    = kol_ss*l_ss;
y_ss    = yol_ss*l_ss;
c_ss    = col_ss*l_ss;
i_ss    = iol_ss*l_ss;

//****************************************
// Model block
//****************************************

model;

    (1/(c))             = (1/(w))*(-1)*(PSI/h)*log(1-h);
    (1/(c))             = BETA*(1/(c(+1)))*(1-DELTA+(r(+1)));
    r                   = ALPHA*(y/k);
    w                   = (1-ALPHA)*((y)/(l));
    y                   = (k)^(ALPHA)*(A)^(1-ALPHA)*(l)^(1-ALPHA);
    k(+1)               = in+(1-DELTA)*k;
    in                  = y-c;
    log(A)              = RHOa*log(A(-1))+epsa; 
    prod                = y/l;
end;

//****************************************
// Steady State Computation
//****************************************

initval;

    l       = l_ss;
    r       = r_ss;
    k       = k_ss;
    w       = w_ss;
    c       = c_ss;
    in      = i_ss;
    y       = y_ss;
    prod    = y/l;
    A       = 1;

end;
steady;
resid;

//****************************************
// Shocks
//****************************************  

shocks;
    var epsa;
    stderr 0.006;
end;

check;
steady;

//****************************************
// Simulation and Statistics
//****************************************  

stoch_simul(order=1,drop=30,periods=192,irf=150, loglinear,hp_filter=1600,tex);

