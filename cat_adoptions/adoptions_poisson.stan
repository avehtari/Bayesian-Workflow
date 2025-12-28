// poisson version
data{
  int N;
  array[N] int adopted; // 1/0 indicator
  vector[N] days;    // days until event
  array[N] int color;   // 1=black, 2=other
}
parameters{
  vector<lower=0>[2] lambda;
}
model{
  lambda ~ exponential(10.0);
  adopted ~ poisson(lambda[color] .* days);
}
