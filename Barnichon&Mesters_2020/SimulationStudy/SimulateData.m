function [vY,vX,vXi,vU] = SimulateData(vParE, dParS, vParM, dS, iT)   % true parameter values vParE : % a_u b_u a_e b_e d_u0 d_u1 d_u2 d_u3 d_e1 d_e2 
                                                                        
vXi      = zeros(iT,1); 
vSigmaXi = zeros(iT,1); 
vNuXi    = zeros(iT,1); 
vRxi     = normrnd(0,1,[iT,1]);

vU       = zeros(iT,1); 
vSigmaU  = zeros(iT,1); 
vNuU     = zeros(iT,1); 
vRu      = normrnd(0,1,[iT,1]); 

vY       = zeros(iT,1);
vX       = zeros(iT,1);
for t=1:iT
   if t==1    
       vXi(t)       = dParS * vRxi(t); 
       
       vSigmaU(t)   = vParE(5) + vParE(6) * vXi(t)^2;
       vNuU(t)      = sqrt(vSigmaU(t)) * vRu(t); 
       vU(t)        = vNuU(t); 
       
       vX(t)        = vXi(t) + vU(t);
       vY(t)        = vParM(3) * vX(t) + vU(t);
   else
       vXi(t)       = vParE(3) * vXi(t-1) + vParE(4) * vNuXi(t-1) + dParS * vRxi(t);   
       
       vSigmaU(t)   = vParE(5) + vParE(6) * vXi(t)^2 + vParE(7) * vSigmaU(t-1) + vParE(8) * vNuU(t-1)^2;  
       vNuU(t)      = sqrt(vSigmaU(t)) * vRu(t); 
       vU(t)        = vParE(1) * vU(t-1) + vParE(2) * vNuU(t-1) + vNuU(t);   
       
       if t==2 
            vX(t)   = vParM(4) * vX(t-1) + vXi(t) + dS * vU(t); 
            vY(t)   = vParM(1) * vY(t-1) + vParM(3) *vX(t) + vU(t);  
       else
            vX(t)   = vParM(4) * vX(t-1) + vParM(5) * vX(t-2) + vXi(t) + dS * vU(t); 
            vY(t)   = vParM(1) * vY(t-1) + vParM(2) * vY(t-2) + vParM(3) * vX(t) + vU(t);  
       end
       
   end
end
end