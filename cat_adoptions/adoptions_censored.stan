// all events, including censored
data{
  int N;
  array[N] int adopted; // 1/0 indicator
  array[N] int days;    // days until event
  array[N] int color;   // 1=black, 2=other
}
parameters{
  vector<lower=0,upper=1>[2] p;
}
model{
  p ~ beta(1,10);
  for (i in 1:N) {
    real P = p[color[i]];
    if (adopted[i]==1) {
      target += log((1-P)^(days[i]-1) * P);
    } else {
      target += log((1-P)^days[i]);
    }
  }
}
