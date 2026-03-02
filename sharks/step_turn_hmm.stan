//
// This Stan program fits an N-state hidden Markov model
// with state-dependent distributions for step length (>=0)
// and turning angle (-pi, pi)
// 
// Parameters to be estimated: 
//   location(mu)/scale(sigma)/zero point mass(mixp) for 
//     step length distributions
//   xangle/yangle for turning angle distributions
//   entries of transition probability matrix (gamma)
//   entries of initial state distribution (initial_dist)
//

data {
  int<lower=0> Nstates;
  int<lower=0> Tlen;
  array[Tlen] int track_index;
  
  vector[Tlen] steplength;
  vector[Tlen] angle;
  
  int<lower=0> Ntracks;
}

parameters {
  //positive_ordered[Nstates] mu;
  vector<lower=0>[Nstates] mu;
  vector<lower=0>[Nstates] sigma;
  vector<lower=0, upper=1>[Nstates] mixp;
  
  vector[Nstates] xangle;
  vector[Nstates] yangle;
  
  array[Nstates] simplex[Nstates] tpm;
  simplex[Nstates] initial_dist;
}

transformed parameters{
  vector[Nstates] shape = mu .* mu ./(sigma .* sigma);
  vector[Nstates] rate = mu ./ (sigma .* sigma);
  matrix[Nstates, Nstates] tpm_matrix;
  vector[Nstates] loc;
  vector[Nstates] kappa;
  for(n in 1:Nstates){
    loc[n] = atan2(yangle[n], xangle[n]);
    kappa[n] = sqrt(xangle[n]*xangle[n] + yangle[n]*yangle[n]); 
  }
  for(i in 1:Nstates)
    for(j in 1:Nstates)
      tpm_matrix[i,j] = tpm[i,j];
}

model {
  vector[Nstates] lalpha;
  vector[Nstates] lalpha_temp;
  array[Nstates] vector[Nstates] log_tpm_tr;
  
  // prior distributions
  mu ~ normal(0.15, 0.05);
  sigma ~ student_t(3, 0, 1);
  xangle ~ normal(0, 1);
  yangle ~ normal(0, 0.5);
  mixp ~ beta(10, 1);

  // log/transpose tpm
  for(i in 1:Nstates)
    for(j in 1:Nstates) 
      log_tpm_tr[j,i] = log(tpm[i,j]);

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
        lalpha_temp[n] = log_sum_exp(to_vector(log_tpm_tr[n]) + lalpha); 
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


generated quantities{
  array[Tlen] vector[Nstates] state_probs;
  array[Tlen] int state_sequence;

  array[Tlen] vector[Nstates] lalpha_mat;
  array[Tlen] vector[Nstates] lbeta_mat;
  array[Nstates] vector[Nstates] log_tpm_tr;
  array[Nstates] vector[Nstates] log_tpm;
  vector[Tlen] pseudo_residuals;
  vector[Nstates] ps_temp;
  array[Ntracks] real llike;
  
  vector[Nstates] lbeta_temp;
  vector[Nstates] ffbs_prob_unnorm;
  vector[Nstates] ffbs_prob_norm;
  
  // log/transpose tpm
  for(i in 1:Nstates){
    for(j in 1:Nstates){ 
      log_tpm_tr[j,i] = log(tpm[i,j]);
      log_tpm[i,j] = log(tpm[i,j]);
    }
  }
  
  // forward-backward algorithm for local state decoding
  for(n in 1:Nstates){
    lalpha_mat[1,n] = log(initial_dist[n]);
    if(steplength[1] > 0)
      lalpha_mat[1,n] = lalpha_mat[1,n] + log(mixp[n]) +
        gamma_lpdf(steplength[1]| shape[n], rate[n]);
    if(steplength[1]==0)
      lalpha_mat[1,n] = lalpha_mat[1,n] + log(1-mixp[n]); 
    if(angle[1] > (-pi())) 
      lalpha_mat[1,n] = lalpha_mat[1,n] + von_mises_lpdf(angle[1] | loc[n], kappa[n]);
  }
  
  for(t in 2:Tlen){
    if(track_index[t]!=track_index[t-1]){
      for(n in 1:Nstates){
        lalpha_mat[t,n] = log(initial_dist[n]);
        if(steplength[t]>0)
          lalpha_mat[t,n] = lalpha_mat[t,n] + log(mixp[n]) +
            gamma_lpdf(steplength[t]| shape[n], rate[n]); 
        if(steplength[t]==0)
          lalpha_mat[t,n] = lalpha_mat[t,n] + log(1-mixp[n]); 
        if(angle[t] > (-pi()))  
          lalpha_mat[t,n] = lalpha_mat[t,n] + 
            von_mises_lpdf(angle[t]| loc[n], kappa[n]);
      }
    } else {      
      for(n in 1:Nstates){
        lalpha_mat[t,n] = log_sum_exp(to_vector(log_tpm_tr[n]) + lalpha_mat[t-1]); 
        if(steplength[t]>0)
          lalpha_mat[t,n] = lalpha_mat[t,n] + log(mixp[n]) +
            gamma_lpdf(steplength[t]| shape[n], rate[n]); 
        if(steplength[t]==0)
          lalpha_mat[t,n] = lalpha_mat[t,n] + log(1-mixp[n]); 
        if(angle[t] > (-pi()))  
          lalpha_mat[t,n] = lalpha_mat[t,n] + 
            von_mises_lpdf(angle[t]| loc[n], kappa[n]);
      }
    }
    if(t == Tlen || track_index[t+1]!=track_index[t])
      llike[track_index[t]]=log_sum_exp(lalpha_mat[t]);
  }

  lbeta_mat[Tlen] = rep_vector(0, Nstates);

  for(tt in 1:(Tlen-1)){
    if(track_index[Tlen - tt+1]!=track_index[Tlen - tt]){
      lbeta_mat[Tlen-tt] = rep_vector(0, Nstates);
    } else{
      for(n in 1:Nstates){
        for(m in 1:Nstates){
          lbeta_temp[m] = log_tpm[n,m]; 
          if(steplength[Tlen-tt+1]>0)
            lbeta_temp[m] = lbeta_temp[m] + log(mixp[m]) +
              gamma_lpdf(steplength[Tlen-tt+1]| shape[m], rate[m]); 
          if(steplength[Tlen-tt+1]==0)
            lbeta_temp[m] = lbeta_temp[m] + log(1-mixp[m]); 
          if(angle[Tlen-tt+1] > (-pi()))  
            lbeta_temp[m] = lbeta_temp[m] + 
              von_mises_lpdf(angle[Tlen-tt+1]| loc[m], kappa[m]);
        }
        lbeta_mat[Tlen-tt,n] = log_sum_exp(lbeta_temp + lbeta_mat[Tlen-tt+1]); 
      }
    }
  }
  
  for(t in 1:Tlen){
    for(n in 1:Nstates){
      state_probs[t,n] = exp(lalpha_mat[t,n] + lbeta_mat[t,n] - llike[track_index[t]]);
    }
  }
  
  // forward-filtering backward-sampling for state sequence samples
  state_sequence[Tlen] = categorical_rng(state_probs[Tlen]);  
  for(t in 2:Tlen){
    if(track_index[Tlen - t+1]!=track_index[Tlen - t + 2]){
      state_sequence[Tlen - t + 1] = categorical_rng(state_probs[Tlen - t + 1]);
    }else{
      ffbs_prob_unnorm = exp(log_tpm_tr[state_sequence[Tlen-t +2]] + lalpha_mat[Tlen-t +1]);
      ffbs_prob_norm = ffbs_prob_unnorm/sum(ffbs_prob_unnorm);
      state_sequence[Tlen - t + 1] = categorical_rng(ffbs_prob_norm);
    }
  }
  
  // forecast pseudo-residuals
  for(tt in 1:Tlen){
    if(tt == 1 || track_index[tt] != track_index[tt-1] || steplength[tt] < 0){
      pseudo_residuals[tt] =  0;      
    } else {
      for(n in 1:Nstates){
          ps_temp[n] = log_sum_exp(to_vector(log_tpm_tr[n]) + lalpha_mat[tt-1]);
        if(steplength[tt]>0)
          ps_temp[n] = ps_temp[n] + bernoulli_lcdf(1 | mixp[n]) + 
            gamma_lcdf(steplength[tt]| shape[n], rate[n]); 
        if(steplength[tt]==0)
          ps_temp[n] = ps_temp[n] + bernoulli_lcdf(0 | mixp[n]);  
        if(angle[tt] > (-pi()))  
          ps_temp[n] = ps_temp[n] + 
              von_mises_lcdf(angle[tt]| loc[n], kappa[n]);
 
      }
        pseudo_residuals[tt] = inv_Phi(exp(log_sum_exp(ps_temp) - 
                                       log_sum_exp(lalpha_mat[tt-1])));
    }
  }
  
}
