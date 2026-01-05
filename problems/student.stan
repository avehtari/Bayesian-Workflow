// student-student model
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
}
model {
  mu ~ student_t(4, 0, 100);
  y ~ student_t(4, mu, 1);
}
