clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This file computes the results for Table II  
%%% The data is simulated from HybridPC.mod 
%%% The estimation results for Table II are stored in: RejectionFrequencies.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dynare HybridPC.mod noclearall

iMs        = 5000;            % number of simulations
vT         = [200 500];       % sample size
vK         = [20 40];         % number of instruments

% true parameter values : \gamma_b, \gamma_f, \lambda, \sigma_m, \rho 
mPtrue     = [0.6 0.3 0.4 0.10 0.0;
              0.6 0.3 0.4 0.25 0.0;
              0.6 0.3 0.4 0.5  0.0;
              0.6 0.3 0.4 1    0.0;
              0.6 0.3 0.4 0.10 0.5;
              0.6 0.3 0.4 0.25 0.5;
              0.6 0.3 0.4 0.5  0.5;
              0.6 0.3 0.4 1    0.5];
           
% store rejection frequencies
mR_IV      = NaN(size(mPtrue,1),size(vT,2)*size(vK,2));
mR_IVa     = NaN(size(mPtrue,1),size(vT,2)*size(vK,2));
mR_AR      = NaN(size(mPtrue,1),size(vT,2)*size(vK,2));
mR_ARa     = NaN(size(mPtrue,1),size(vT,2)*size(vK,2));
mR_ARas    = NaN(size(mPtrue,1),size(vT,2)*size(vK,2));

for ss = 1:size(vT,2)
    
    iT         = vT(ss);                    % select sample size 
    iMaxLag    = floor(4*(iT/100)^(2/9))+1; % lag length
 
    vQS        = [1 ; zeros(iT-1,1)];       % quadratic spectral kernel 
    for l=2:iMaxLag
        dX     = l / (1 + iMaxLag);
        vQS(l) = (25 / (12 * pi^2 * dX^2)) * ( sin(6*pi * dX /5)/ (6*pi * dX /5) - cos(6*pi * dX /5) );
    end
    mK         = toeplitz(vQS);           % toeplitz weight matrix 
    
    for ll = 1:size(vK,2)
        for pp = 1:size(mPtrue,1)
                        
            iK         = vK(ll);          % select number of structural shocks 
            vPtrue     = mPtrue(pp,1:3)'; % true parameters to test  
     
            vR_IV      = zeros(1,iMs);    % storage for current case 
            vR_IVa     = zeros(1,iMs);
            vR_AR      = zeros(1,iMs);
            vR_ARa     = zeros(1,iMs);
            vR_ARas    = zeros(1,iMs);
    
            set_param_value('dSmp',mPtrue(pp,4));    % set the variance of the structural shock for dynare 
            set_param_value('rho_cp',mPtrue(pp,5));  % set the serial correlation for the cost-push shock for dynare 
    
            for kk = 1:iMs
                stoch_simul(var_list_);                              % generates data from basic model using dynare 
  
                vY     = ppi(iK:iT+iK-1);                                % dependent variable  
                mX     = [ppi(iK-1:iT+iK-2) ppi(iK+1:iT+iK) y_gap(iK:iT+iK-1)]; % endogenous variables  
                mZ     = NaN(iT,iK) ;                             % construct lagged shocks     
                for kkl = 1:iK
                    mZ(:,kkl) = mpshock(kkl:iT+kkl-1);
                end
                mZ     = fliplr(mZ);
            
                vR     = 0:1:iK-1;
                mI     = [sum(mZ,2) sum(mZ .* vR,2) sum(mZ .* vR.^2,2)];  % polynomial instrument matrix 
                mIIinv = inv(mI'*mI);
            
                % IV
                mPz       = mZ * inv(mZ' * mZ) * mZ';
                vBeta     = inv(mX'*mPz*mX) * mX'*mPz*vY;
                vRes      = vY - mX * vBeta;
                mCov      = NeweyWest(vRes,mZ,-1, 0);    
                mOmega    = inv(mX' * mPz * mX) * mX' * mZ * inv(mZ' * mZ) * mCov * inv(mZ'*mZ) * mZ' * mX * inv(mX' * mPz * mX);
                vR_IV(kk) = ( (vBeta-vPtrue)' / mOmega ) * (vBeta-vPtrue)  > chi2inv(0.95,iK);
        
                % IV_a
                mPi       = mI * mIIinv * mI'; 
                vBeta     = inv(mX'*mPi*mX) * mX'*mPi*vY;
                vRes      = vY - mX * vBeta;
                mCov      = NeweyWest(vRes,ones(iT,1),-1, 0);
                mOmega    = (iT^-1 * mCov) * inv(mX' * mPi * mX);
                vR_IVa(kk)= ( (vBeta-vPtrue)' / mOmega ) * (vBeta-vPtrue)  > chi2inv(0.95,3);   
 
                % AR
                vYtilde           = vY - mX * vPtrue;	
                vTheta            = inv(mZ' * mZ) * mZ' * vYtilde;
                vRes              = vYtilde - mZ * vTheta;
                mCov              = NeweyWest(vRes,mZ,-1, 0);
                dWaldStat         = vYtilde' * mZ * inv(mCov) * mZ' * vYtilde;               
                vR_AR(kk)         = dWaldStat > chi2inv(0.95,iK);  
                
                % AR_a
                vYtilde           = vY - mX * vPtrue;	
                [EstCov,se,vBeta] = hac(ones(iT,1),vYtilde,'Intercept',false,'weights','QS','bandwidth',floor(4*(iT/100)^(2/9))+1,'display','off');
                dWaldStat         = vYtilde' * mI * inv( mI'* mI * iT * EstCov) * mI' * vYtilde;               
                vR_ARa(kk)        = dWaldStat > chi2inv(0.95,3);      
                
                % AR_as
                dARStat           = fSubSet2(vPtrue(3), vY, mX(:,3), mX(:,1:2), mI, iT, mK); 
                vR_ARas(kk)       = dARStat > chi2inv(0.95,1); 
            end
            mR_IV(pp,(ss-1)*size(vK,2)+ll)    = mean(vR_IV,2);
            mR_IVa(pp,(ss-1)*size(vK,2)+ll)   = mean(vR_IVa,2);
            mR_AR(pp,(ss-1)*size(vK,2)+ll)    = mean(vR_AR,2);
            mR_ARa(pp,(ss-1)*size(vK,2)+ll)   = mean(vR_ARa,2);
            mR_ARas(pp,(ss-1)*size(vK,2)+ll)  = mean(vR_ARas,2);
           
        end
    end 
end

save('RejectionFrequencies.mat','mR_IV','mR_IVa','mR_AR','mR_ARa','mR_ARas');
