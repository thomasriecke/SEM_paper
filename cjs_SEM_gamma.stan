
functions {
  // this function calculates the probability of never being seen again given an individual's last encounter and it's survival probability  
  matrix prob_uncaptured(int nI, int nT, real p, vector phi) {
    matrix[nI, nT] chi;
    
    for (i in 1 : nI) {
      chi[i, nT] = 1.0;
      for (t in 1 : (nT - 1)) {

        int t_curr = nT - t;
        int t_next = t_curr + 1;
        
        /*
          int t_curr;
        int t_next;
        
        t_curr = nT - t;
        t_next = t_curr + 1;
        
        t_curr = nT - t;
        t_next = t_curr + 1;
        
        */

        chi[i, t_curr] = (1 - phi[i]) + phi[i] * (1 - p) * chi[i, t_next];
      }
    }
    return chi;
  }
}

data {
  int nI;
  int nT;
  array[nI] int f;
  array[nI] int l;
  array[nI, nT, 2] int<lower=-10> y; // Capture-history, -9 used to fill cells that should not be sampled, throwing trap if indexing is incorrect
}


parameters {
  vector[nI] z;
  real<lower=0> sigma;
  real lambda;
  vector[2] mu;
  real<lower = 0, upper = 1> p;
}

transformed parameters {
  vector[nI] phi;
  vector[nI] gamma;
  vector[nI] h;
  vector[nI] epsilon;
  matrix<lower=0, upper=1>[nI, nT] chi;
  epsilon = z * sigma; 

  h = exp(mu[1] + epsilon * lambda);
  phi = exp(-h);
  gamma = exp(mu[2] + epsilon);  
  
  chi = prob_uncaptured(nI, nT, p, phi);
}

model {

  mu ~ normal(-0.5,2);
  sigma ~ gamma(1,1);
  z ~ normal(0, 1);
  lambda ~ normal(0, 2.5);
  
  // Likelihood
  for (i in 1 : nI) {
    for (t in (f[i] + 1) : l[i]) {
      1 ~ bernoulli(phi[i]);
      y[i, t, 1] ~ bernoulli(p);
    }
    1 ~ bernoulli(chi[i, l[i]]);

    for (t in f[i]:nT){
      if (y[i,t,1] == 1){
        y[i,t,2] ~ poisson(gamma[i]);
      }
    }

  }
}

generated quantities {
  real cov;
  cov = lambda * sigma * sigma;
}
