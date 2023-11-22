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
   
   gbar = moment(pars,tau,Z,y,X,a,alpha);
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
  real <lower=0> scale_global ; // scale for the half -t prior for tau
  real <lower=1> nu_global ; // degrees of freedom for the half -t priors for tau
  real <lower=1> nu_local ; // degrees of freedom for the half - t priors for lambdas
  real <lower=0> slab_scale ; // slab scale for the regularized horseshoe
  real <lower=0> slab_df ; // slab degrees of freedom for the regularized horseshoe
  real a;
  real <lower=0> lambda_gallant;
}

parameters {
  //vector[K] beta;
 // real<lower=20,upper=100> a;
  //vector<lower=0>[K] lambda;
  //real<lower=0> tau1;
  real alpha;
  vector[K] z;
real <lower=0> aux1_global ;
real <lower=0> aux2_global ;
vector <lower=0>[K] aux1_local ;
vector <lower=0>[K] aux2_local ;
real <lower=0> caux ;
}

transformed parameters {
real <lower=0> tau1; // global shrinkage parameter
vector<lower=0>[K] lambda ; // local shrinkage parameter
vector<lower=0>[K] lambda_tilde ; // ' truncated ' local shrinkage parameter
real <lower=0> c; // slab scale
vector[K] beta ; // regression coefficients
lambda = aux1_local .*sqrt(aux2_local);
tau1 = aux1_global * sqrt ( aux2_global ) * scale_global;
c =slab_scale*sqrt(caux);
lambda_tilde =sqrt(c^2*square(lambda) ./ (c^2+tau1^2*square(lambda)));
beta =z .*lambda_tilde*tau1*lambda_gallant;
}

model {
  z ~ normal(0,1); // half -t priors for lambdas and tau , and inverse - gamma for c ^2
aux1_local ~ normal(0,1);
aux2_local ~  inv_gamma(0.5*nu_local,0.5*nu_local);
aux1_global ~ normal(0,1);
aux2_global ~ inv_gamma(0.5*nu_global,0.5*nu_global);
caux ~ inv_gamma(0.5*slab_df,0.5*slab_df);

  alpha ~ normal(0,20); // variance was at 5 before

//  a ~ uniform(20,100);
  target += kernel(beta,y,Z,X,tau,mu,a,alpha);
}

