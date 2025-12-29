//
// This Stan program fits an N-state hidden Markov model
// with state-dependent distributions for step length (>0)
// and turning angle (-pi, pi)
// 
// Parameters to be estimated: 
//   location(mu)/scale(sigma) for step length distributions
//   location(kappa)/concentration for turning angle
//   entries of transition probability matrix (gamma)
//

data {
  int<lower=0> Nstates;
  int<lower=0> Tlen;
  array[Tlen] int track_index;
  
  vector[Tlen] steplength;
  vector[Tlen] angle;
  
  int nCovs;
  matrix[Tlen, nCovs+1] covs;

}

parameters {
  vector<lower=0>[Nstates] mu;
  vector<lower=0>[Nstates] sigma;
  vector<lower=0, upper=1>[Nstates] mixp;
  
  vector[Nstates] xangle;
  vector[Nstates] yangle;
  
  matrix[Nstates*(Nstates-1), nCovs+1] beta;
  simplex[Nstates] initial_dist;

}

transformed parameters{
  vector[Nstates] shape = mu .* mu ./(sigma .* sigma);
  vector[Nstates] rate = mu ./ (sigma .* sigma);
  vector[Nstates] loc;
  vector[Nstates] kappa;
  for(n in 1:Nstates){
    loc[n] = atan2(yangle[n], xangle[n]);
    kappa[n] = sqrt(xangle[n]*xangle[n] + yangle[n]*yangle[n]); 
  }
}

model {
  vector[Nstates] lalpha;
  vector[Nstates] lalpha_temp;
  array[Tlen] matrix[Nstates, Nstates] tpm;
  array[Tlen] matrix[Nstates, Nstates] log_tpm;
  array[Tlen] matrix[Nstates, Nstates] log_tpm_tr;
  
  // prior distributions
  mu ~ normal(0.15, 0.05);
  sigma ~ student_t(3, 0, 1);
  xangle ~ normal(0, 1);
  yangle ~ normal(0, 0.5);
  mixp ~ beta(10, 1);
  for(n in 1:(Nstates*(Nstates -1))){
    beta[n] ~ normal(0, 1);
  }


  // derive array of (log-)transition probabilities
  for(t in 1:Tlen){
    int betarow = 1;
    for(i in 1:Nstates){
      for(j in 1:Nstates){
        if(i == j){
          tpm[t,i,j] = 1;
        } else {
          tpm[t,i,j] = exp(beta[betarow]*to_vector(covs[t]));
          betarow = betarow + 1;
        }
      }
    }
    
    for(n in 1:Nstates){
      log_tpm[t][n] = log(tpm[t][n]/sum(tpm[t][n]));
    }
  }
  
  // transpose
  for(t in 1:Tlen){
    for(i in 1:Nstates){
      for(j in 1:Nstates){
        log_tpm_tr[t,j,i] = log_tpm[t,i,j];
      }
    }
  }

  // likelilhood computation
  for(n in 1:Nstates){
    lalpha[n] = log(initial_dist[n]);
    if(steplength[1] > 0)
      lalpha[n] = lalpha[n] + log(mixp[n]) +
        gamma_lpdf(steplength[1]| shape[n], rate[n]);
    if(steplength[1]==0)
      lalpha[n] = lalpha[n] + log(1-mixp[n]); 
    if(angle[1] > (-pi())) 
      lalpha[n] = lalpha[n] + von_mises_lpdf(angle[1] | loc[n], kappa[n]);
  }
  
  for(t in 2:Tlen){
    if(track_index[t]!=track_index[t-1]){
      for(n in 1:Nstates){
        lalpha_temp[n] = log(initial_dist[n]);
        if(steplength[t]>0)
          lalpha_temp[n] = lalpha_temp[n] + log(mixp[n]) +
            gamma_lpdf(steplength[t]| shape[n], rate[n]); 
        if(steplength[t]==0)
          lalpha_temp[n] = lalpha_temp[n] + log(1-mixp[n]); 
        if(angle[t] > (-pi()))  
          lalpha_temp[n] = lalpha_temp[n] + 
            von_mises_lpdf(angle[t]| loc[n], kappa[n]);
      }
    } else {
      for(n in 1:Nstates){
        lalpha_temp[n] = log_sum_exp(to_vector(log_tpm_tr[t,n]) + lalpha); 
        if(steplength[t]>0)
          lalpha_temp[n] = lalpha_temp[n] + log(mixp[n]) +
            gamma_lpdf(steplength[t]| shape[n], rate[n]); 
        if(steplength[t]==0)
          lalpha_temp[n] = lalpha_temp[n] + log(1-mixp[n]); 
        if(angle[t] > (-pi()))  
          lalpha_temp[n] = lalpha_temp[n] + 
            von_mises_lpdf(angle[t]| loc[n], kappa[n]);
      }
    }  
    lalpha = lalpha_temp;
    if(t == Tlen || track_index[t+1]!=track_index[t])
      target+=log_sum_exp(lalpha);
  }
}

