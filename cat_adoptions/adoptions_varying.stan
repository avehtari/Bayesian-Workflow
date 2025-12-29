// cats vary in their adoption probabilities
data{
  int N;
  array[N] int adopted; // 1/0 indicator
  array[N] int days;    // days until event
  array[N] int color;   // 1=black, 2=other
}
parameters{
  // average adoptions
  vector<lower=0,upper=1>[2] p;
  vector<lower=0>[2] theta; // dispersion
  array[N] vector<lower=0,upper=1>[2] q; // cat specific probabilities
}
model{
  p ~ beta(1,10);
  theta ~ exponential(1);
  for (i in 1:N) {
    real P = 0;
    for (j in 1:2)
      q[i,j] ~ beta(p[j]*theta[j], (1-p[j])*theta[j]);
    P = q[i,color[i]];
    if (adopted[i]==1) {
      target += log((1-P)^(days[i]-1) * P);
    } else {
      target += log((1-P)^days[i]);
    }
  }
}
