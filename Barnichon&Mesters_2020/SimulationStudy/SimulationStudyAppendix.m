clear all
close all
clc


iMs        = 5000;                                  % number of simulations for each design 
vT         = [200 500];                             % sample sizes
vK         = [5 10 20 40 80];                       % number of polynomial shock instruments
vSE        = [1];                                   % degree of engogeneity 


mParS      = [1.4 -0.5 1.0 1.4 -0.5];               % true model parameters 
           
vParSigma  = [0.1 0.25 0.5 1];                      % variance structural shocks 

mError     = [0.0 0.0 0.0 0.0 1.0  0.0 0.0 0.0  ;   % error term coefficients 
              0.0 0.0 0.0 0.0 0.0  1.0 0.0 0.0  ;
              0.0 0.0 0.0 0.0 0.05 0.0 0.9 0.05 ;
              0.5 0.5 0.5 0.5 1.0  0.0 0.0 0.0  ;
              0.5 0.5 0.5 0.5 0.0  1.0 0.0 0.0  ;
              0.5 0.5 0.5 0.5 0.05 0.0 0.9 0.05 ]; 
              

% store rejection frequencies
mR_ARa      = NaN( size(vK,2),size(vParSigma,2),size(vT,2),size(vSE,2),size(mParS,1),size(mError,1));
mR_ARas     = NaN( size(vK,2),size(vParSigma,2),size(vT,2),size(vSE,2),size(mParS,1),size(mError,1));

for ss = 1:size(vT,2)
    iT         = vT(ss);                       % select sample size 
    iMaxLag    = floor(4*(iT/100)^(2/9))+1;    % lag length   
    vQS        = [1 ; zeros(iT-1,1)];          % quadratic spectral kernel 
    for l=2:iMaxLag +1 
        dX     = l / (1 + iMaxLag);
        vQS(l) = (25 / (12 * pi^2 * dX^2)) * ( sin(6*pi * dX /5)/ (6*pi * dX /5) - cos(6*pi * dX /5) );
    end
    mK         = toeplitz(vQS);                % toeplitz weight matrix 
    for ll = 1:size(vK,2)
        iK         = vK(ll);                   % select number of structural shocks  
        for sse = 1:size(vSE,2)                % degree of endogeneity 
            for ppt = 1:size(mParS,1)          % true structural parameters 
                for ppss = 1:size(vParSigma,2) % strength structural shocks 
                    for ppe = 1:size(mError,1) % design structural shocks 
                        vR_ARa     = zeros(1,iMs);  
                        vR_ARas    = zeros(1,iMs);
    
                        for kk = 1:iMs  % simulations 
                            vPtrue       = mParS(ppt,:)';
                            
                            [vYt,vX,vXi,vU] = SimulateData(mError(ppe,:), vParSigma(ppss), mParS(ppt,:), vSE(sse), iT + iK + 1);
               
                            vY  = vYt(iK+2:end);
                            mW  = [vYt(iK+1:end-1) vYt(iK:end-2) vX(iK+2:end)];
                                      
                            mZ     = NaN(iT,iK) ;                             % construct lagged shocks     
                            for kkl = 1:iK
                                mZ(:,kkl) = vXi(kkl:iT+kkl-1);
                            end
                            mZ                   = fliplr(mZ);

                            vR                   = 0:1:iK-1;
                            mI                   = [sum(mZ,2) sum(mZ .* vR,2) sum(mZ .* vR.^2,2)];  % polynomial instrument matrix 
               
                            vYtilde              = vY - mW * vPtrue(1:3);	
                            dWaldStat            = vYtilde' * mI * inv( mI'* mI) * mI' * vYtilde / (iT^-1 * vYtilde'* mK * vYtilde);               
                            vR_ARa(kk)           = dWaldStat > chi2inv(0.95,3);     

                            dARStat              = fSubSet2(vPtrue(3), vY, mW(:,3), mW(:,1:2), mI, iT, mK); 
                            vR_ARas(kk)          = dARStat > chi2inv(0.95,1);                     
                        end
                    mR_ARa(ll , ppss, ss, sse, ppt, ppe)    = mean(vR_ARa,2);
                    mR_ARas(ll, ppss, ss, sse, ppt, ppe)    = mean(vR_ARas,2);
                    end
                end
            end
        end
    end 
end

save('RejectionFrequenciesLARGE.mat','mR_ARa','mR_ARas');