/*
 * Generate data from a simple bivariate model with a PC and output gap
 */

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------


var ppi          
    y_gap       
    mpshock
    cp 
    mp
;     

varexo  eps_nu     
        eps_cp
       ;

parameters
    aalpha 
    betta    
    ggamma_b           
    pphi1
    pphi2 
    kkappa
    rho_cp
    rho_mp
    dSmp 
    ;
%----------------------------------------------------------------
% Parametrization, p. 67  and p. 113-115
%----------------------------------------------------------------
aalpha   = -1;
betta    = 0.3;
ggamma_b = 0.6;
kkappa   = 0.4;
pphi1    = 1.2;
pphi2    = -0.4;
rho_cp   = 0;
rho_mp   = 0; 
dSmp     = 0.1;

%----------------------------------------------------------------
% First Order Conditions
%----------------------------------------------------------------

model(linear); 

[name='New Keynesian Phillips Curve eq. (22)']
ppi = ggamma_b * ppi(-1) + betta*ppi(+1) + kkappa*y_gap + cp;

[name='Dynamic IS Curve eq. (23)']
y_gap = pphi1 * y_gap(-1) + pphi2 * y_gap(-2) + mp + aalpha * cp;  

[name='cost push shock']
cp      =  rho_cp * cp(-1 ) +  eps_cp;

mpshock = dSmp * eps_nu;

mp      = rho_mp * mp(-1) + mpshock;
end;


%----------------------------------------------------------------
%  define shock variances
%---------------------------------------------------------------
shocks;
    var eps_nu  = 1^2; //1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    var eps_cp  = 1^2; //unit shock to preferences
end;

%----------------------------------------------------------------
%  steady states: all 0 due to linear model
%---------------------------------------------------------------
steady;
check;

%stoch_simul(order=1,irf=30);

stoch_simul(order = 1, nocorr, nomoments, IRF = 0, periods = 1000);