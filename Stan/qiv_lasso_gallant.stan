functions {

// Helper Function that creates columnwise means
  vector col_sums(matrix Z) {
    
    vector[cols(Z)] s ;		
    
    for (j in 1:cols(Z)){
      s[j] = sum(col(Z, j));
    }
    
    return s;
  }
  

// Helper Function for the quantile loss
vector q_loss(real q, vector u, real a){
  vector[rows(u)] temp;
  
  for (i in  1:rows(u)){
    temp[i] = (1 - ( 1/(1+exp(-a*u[i]))));
  }
  
        return temp;
    }
  
// Sample Moments
  vector moment(vector pars, real tau, matrix Z, vector y, matrix X, real a, real alpha){
    vector[rows(X)] tloss;
    matrix[rows(Z),cols(Z)] temp;
   // vector[rows(pars)] ret_val;
    
    tloss = q_loss(tau, y-X*pars - alpha,a);
    
    temp = diag_pre_multiply(rep_vector(tau,rows(Z)) - tloss,Z)./rows(X);
    //ret_val = col_sums(temp); 
    //  return ret_val;
    return col_sums(temp); 
      }
      
// Weight Matrix
  matrix weight_m(real tau,vector pars, vector y, matrix Z, matrix X, real a, real alpha){
    vector[rows(y)] tloss;
    matrix[cols(Z),cols(Z)] W;
    matrix[cols(Z),cols(Z)] identity;
    matrix[rows(X),cols(Z)] temp;
    
    
    identity = diag_matrix(rep_vector(1,cols(Z))); 
    
    tloss = q_loss(tau, y-X*pars - alpha,a);
    
    temp = diag_pre_multiply(rep_vector(tau,rows(Z)) - tloss,Z);
    
    return ((temp)'*(temp)./rows(Z))\identity;//inverse((temp)'*(temp)./rows(Z));//((temp)'*(temp)./rows(Z))\identity;
    
  }
  
// Objective Function
 real ObjFn(vector pars, vector y, matrix Z,matrix X, real tau, real a, real alpha){
   vector[cols(Z)] gbar;
   matrix[cols(Z),cols(Z)] W;
   
   gbar = moment(pars,tau,Z,y,X,a, alpha);
   W = weight_m(tau,pars,y,Z,X,a,alpha);
   return -(gbar'*W*gbar)*rows(Z)/2;
 }

// Create Log Posterior Kernel
real kernel(vector pars, vector y, matrix Z,matrix X, real tau, vector mu, real a, real alpha){
  vector[cols(pars)] lprior_temp;
  real lprior;
  

  
  return ObjFn(pars,y,Z,X,tau,a,alpha);// + lprior;
}
}

data {
    int N;                 // Number of observation
    int K;                 // Number of predictors
    int L;                 // Number of Instruments
    real tau;               // Quantile to be estimated
    vector[K] mu;          // Prior on quantile coefficients
    vector[N] y;           // Response variable
    matrix[N,K] X;        // Design matrix
    matrix[N,L] Z;        // Instrument matrix
    //real<lower=0> lasso_df;  // prior degrees of freedom
    //real<lower=0> lasso_scale;  // prior scale
    real<lower=0> lambda;
    vector[K] w_hat;
    real a;
}

parameters {
    vector[K] beta;
    real alpha;
    //real<lower=20,upper=100> a;
    //real<lower=0> lasso_lambda;
    //vector<lower=0>[K] w_hat;
}

model {
    alpha ~ normal(0,5);
  beta ~ double_exponential(mu,lambda * w_hat);
  target += kernel(beta,y,Z,X,tau,mu,a,alpha);
}

