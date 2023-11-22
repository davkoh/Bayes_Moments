clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file computes the results for Table IV & Figure V 
%%% Change iForce to select the forcing variable 
%%% The estimation results for Table IV are stored in:
%%% mResR: for the   restricted parameter estimates [ cue_est , lower 90, upper 90, lower 95 , upper 95 ] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mData              = xlsread('Data_QJE.xlsx');

vPia               = mData(:,4);                                  % inflation         (4 == pix)

iForce             = 0;                                           % select the forcing variable: 0 is ygap -- 1 is ugap 

if iForce == 0
    vXa            = mData(:,5);                                  % selects ygap 
    vScale         = [-0.5 2];                                    % x-axis for figures        
end
if iForce == 1
    vXa            = mData(:,6);                                  % selects ugap 
    vScale         = [-2 0.5];                                    % x-axis for figures       
end 

vPiF               = 0.25 * (lagmatrix(vPia,-1) + lagmatrix(vPia,-2) + lagmatrix(vPia,-3) + lagmatrix(vPia,-4)); % future inflation
vPiL               = 0.25 * (lagmatrix(vPia, 1) + lagmatrix(vPia, 2) + lagmatrix(vPia, 3) + lagmatrix(vPia, 4)); % lagged inflation
	
vIV                = mData(:,8) + mData(:,10);                        % instruments (HFI) sum of FF4 and onrun10 
iL                 = 24;                                              % number of instruments 

mZ                 = [vIV NaN(size(vPia,1),iL)];                      % construct instrument matrix 
for j = 1:iL
	mZ(:,j+1)      = lagmatrix(vIV,j);
end 

mObs               = rmmissing([vPia vPiL vPiF vXa mZ]);          % group all observations and drop missing values 

iT                 = size(mObs,1);                                    % effective sample size  

vY                 = mObs(:,1) - mean(mObs(:,1),1);                   % demeaned dependent var 
mW                 = mObs(:,2:4) - mean(mObs(:,2:4),1);               % demeaned endogenous vars

vYT                = vY - mW(:,1);		                              % dependent var under \gamma_b + \gamma_f = 1								
mWT                = [(mW(:,2) - mW(:,1)) mW(:,3)];                   % endogenous vars	under \gamma_b + \gamma_f = 1

mZ                 = mObs(:,5:5+iL);                                  % instruments  
vR                 = 0:1:iL;
mZp                = [sum(mZ,2) sum(mZ .* vR,2) sum(mZ .* vR.^2,2)];  % polynomial instrument matrix 

mZl                = mObs(:,5+iL+1:end) - mean(mObs(:,5+iL+1:end),1);

iMaxLag            = floor(4*(iT/100)^(2/9))+1;                       % lag length   
vQS                = [1 ; zeros(iT-1,1)];                             % quadratic spectral kernel 
for l=2:iMaxLag
    dX             = (l-1) / (1 + iMaxLag);
    vQS(l)         = (25 / (12 * pi^2 * dX^2)) * ( sin(6*pi * dX /5)/ (6*pi * dX /5) - cos(6*pi * dX /5) );
end
mK                 = toeplitz(vQS);                                       % toeplitz weight matrix 

% baseline IV estimates 
vDelta             = inv(mZp'*mW) * mZp' * vY;                        % simple iv estimates  
vInit              = vDelta(2:3);
[vDeltaR dCueValue]  = fCUErestricted(vInit, vYT, mWT, mZp, iT, mK);  % cue estimates under restriction 

%Compute subset confidence bounds based on AR_{a,s} statistic 
mResR              = [vDeltaR zeros(2,4)];                            % storage restricted       
dGS                = 0.01;                                            % grid size 

vPG                = (-20:dGS:20)';                                   % grid inflation expectations 
vAR_res            = zeros(size(vPG,1),1); 
for i = 1:size(vPG,1)
    vAR_res(i)     = fSubSet2(vPG(i), vYT,  mWT(:,1), mWT(:,2), mZp, iT, mK);   % subset statistic under restriction    
end
mResR(1,2:5)       = [vPG(find(vAR_res < chi2inv(0.90,1),1,'first')) vPG(find(vAR_res < chi2inv(0.90,1),1,'last')) vPG(find(vAR_res < chi2inv(0.95,1),1,'first')) vPG(find(vAR_res < chi2inv(0.95,1),1,'last'))]; 

vPF                = (-20:dGS:20)';                                   % grid forcing variable 
vAR_res            = zeros(size(vPF,1),1); 
for i = 1:size(vPF,1)
    vAR_res(i)     = fSubSet2(vPF(i), vYT, mWT(:,2), mWT(:,1) , mZp, iT, mK);   % subset statistic under restriction  
end
mResR(2,2:5)       = [vPF(find(vAR_res< chi2inv(0.90,1),1,'first')) vPF(find(vAR_res< chi2inv(0.90,1),1,'last')) vPF(find(vAR_res< chi2inv(0.95,1),1,'first')) vPF(find(vAR_res< chi2inv(0.95,1),1,'last'))]; 


% joint region expected inflation and forcing variable 
vPG                = (-0.5:dGS:2)';                  % range for inflation expectation
vPL                = (vScale(1):dGS:vScale(2))';     % range for forcing variable
mARres             = zeros(size(vPG,1),size(vPL,1)); % storage AR statistic
mTestres           = zeros(size(vPG,1),size(vPL,1)); % storage test outcome
for i = 1:size(vPG,1)
    for j=1:size(vPL,1)
        vU         = vYT - mWT * [vPG(i) vPL(j)]'; 
        mPz        = mZp * inv(mZp'*mZp) * mZp' ;
        mMz        = eye(iT) - mPz;
        mARres(i,j)= vU' * mPz * vU / (iT^-1 * vU' * mMz * mK * mMz * vU);   % AR test statistic 
        
        if mARres(i,j) < chi2inv(0.95,2)         % 95% confidence
            mTestres(i,j)   = 1;
        end 
        if mARres(i,j) < chi2inv(0.90,2)         % 90% confidence 
            mTestres(i,j)   = 2;
        end 
        if mARres(i,j) < chi2inv(0.67,2)         % 67% confidence 
            mTestres(i,j)   = 3;
        end    
        
    end
end

map=[];                        % set shades of grey  
map(1,:) = [1 1 1 ];
map(2,:) = [.9 .9 .9];
map(3,:) = [.7 .7 .7];
map(4,:) = [.3 .3 .3];

figure,             % make the figure 
figure1 = figure(1);
set(0,'defaulttextInterpreter','latex')
imagesc(vPL,vPG,mTestres),colormap(figure1, map);
axis('xy')
line([0 0],[-10 10],'LineStyle','--','Color',[.5 .5 .5],'LineWidth',.5);
line([-10 10],[0 0],'LineStyle','--','Color',[.5 .5 .5],'LineWidth',.5);
hold on
h(1) =plot(vDeltaR(2),vDeltaR(1), 'o','MarkerSize',5,'Color',[.85 .33 .1],'MarkerFaceColor',[.85 .33 .1],'LineWidth',1);
%lg=legend(h,{['IV_a-',char(949)]},'Location','NorthEast','Fontsize',6)
box off 
title('Confidence region ($AR_{a}$ test)')
if iForce == 0
    xlabel('$\lambda_{\rm Y}$');
end
if iForce == 1
    xlabel('$\lambda_{\rm U}$');
end
ylabel('$\gamma_f$');
set(gca,'linewidth',.15)
set(gca, 'xtick',-6:1:6)
set(gca, 'ytick',-6:1:6,'LineWidth',0.15)
set(gca,'Layer','top')

if iForce == 0
    print('FiguresQJE/ConfSets_HFI_Y_res','-depsc','-painters');
end
if iForce == 1
    print('FiguresQJE/ConfSets_HFI_U_res','-depsc','-painters');
end




