clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The file provides a sample code for computing the IV estimates and AR_{a,s} based confidence bounds 
%%% of a structural equation using polynomial instruments. 
%%%
%%% The equation considered is of the form: 
%%% y_t = \gamma_b y_{t-1} + \gamma_f y_{t+1} + \lambda x_t + u_t 
%%% 
%%%  
%%% The IV estimates and AR_{a,s} based confidence bounds for \gamma_b,
%%% \gamma_f and \lambda are stored in the rows of 
%%% mRes : [ivEst , lower 90,  upper 90, lower 95, upper 95]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mData          = xlsread('.....');                                % load your data ordered as [y, x, \xi] = [outcome, forcing, instruments]

vY             = mData(:,1);                                      % current y_t
vYf            = lagmatrix(vY,-1);                                % expected y_t+1 (can replace by any expected value)
vYl            = lagmatrix(vY,1);                                 % lag y_t-1
vX             = mData(:,2);                                      % current x_t  
vI             = mData(:,3);                                      % current z_t 

iL             = 20;                                              % number of instruments 

mZ             = [vI NaN(size(vY,1),iL)];                         % construct lagged shock instrument matrix 
for j = 1:iL
	mZ(:,j+1)  = lagmatrix(vI,j);
end 

mObs           = rmmissing([vY vYl vYf vX mZ]);                   % group all observations and drop missing values 

iT             = size(mObs,1);                                    % effective sample size  

vY             = mObs(:,1) - mean(mObs(:,1),1);                   % demeaned dependent var 
mW             = mObs(:,2:4) - mean(mObs(:,2:4),1);               % demeaned endogenous vars

mZ             = mObs(:,5:5+iL);                                  % instruments  
vR             = 0:1:iL;
mZp            = [sum(mZ,2) sum(mZ .* vR,2) sum(mZ .* vR.^2,2)];  % polynomial instrument matrix 

iMaxLag        = floor(4*(iT/100)^(2/9))+1;                       % lag length   
vQS            = [1 ; zeros(iT-1,1)];                             % init quadratic spectral kernel 
for l=2:iMaxLag 
    dX         = (l-1) / (1 + iMaxLag);
    vQS(l)     = (25 / (12 * pi^2 * dX^2)) * ( sin(6*pi * dX /5)/ (6*pi * dX /5) - cos(6*pi * dX /5) );  % quadratic spectral kernel
end
mK         = toeplitz(vQS);                                       % toeplitz weight matrix  


% simple iv estimates  
vDelta         = inv(mZp'*mW) * mZp' * vY;                        

% subset confidence bounds based on AR_{a,s} 
mRes           = [vDelta zeros(3,4)];                             % storage
dGS            = 0.01;                                            % grid size 
 
vPL            = (-10:dGS:10)';                                   % grid lag 
vAR            = zeros(size(vPL,1),1);
for i = 1:size(vPL,1)  
     vAR(i)    = fSubSet2(vPL(i), vY, mW(:,1), mW(:,2:3), mZp, iT, mK);   % subset statistic  
end
mRes(1,2:5)    = [ vPL(find(vAR < chi2inv(0.9,1),1,'first')) vPL(find(vAR < chi2inv(0.9,1),1,'last')) vPL(find(vAR < chi2inv(0.95,1),1,'first')) vPL(find(vAR < chi2inv(0.95,1),1,'last'))]; % find bounds 

vPG            = (-10:dGS:10)';                                   % grid expectations 
vAR            = zeros(size(vPG,1),1);
for i = 1:size(vPG,1)
    vAR(i)     = fSubSet2(vPG(i), vY, mW(:,2), [mW(:,1) mW(:,3)], mZp, iT, mK);  % subset statistic   
end
mRes(2,2:5)    = [vPG(find(vAR < chi2inv(0.9,1) ,1,'first')) vPG(find(vAR < chi2inv(0.9,1),1,'last')) vPG(find(vAR < chi2inv(0.95,1) ,1,'first')) vPG(find(vAR < chi2inv(0.95,1),1,'last'))]; 

vPF            = (-10:dGS:10)';                                   % grid forcing variable 
vAR            = zeros(size(vPF,1),1);
for i = 1:size(vPF,1)
    vAR(i)     = fSubSet2(vPF(i), vY , mW(:,3) , mW(:,1:2), mZp, iT, mK); % subset statistic 
end
mRes(3,2:5)    = [vPF(find(vAR< chi2inv(0.9,1),1,'first')) vPF(find(vAR< chi2inv(0.9,1),1,'last')) vPF(find(vAR< chi2inv(0.95,1),1,'first')) vPF(find(vAR< chi2inv(0.95,1),1,'last'))];


