data {
  
  int ns; //number of studies in the network
  int nvar[ns]; // number of params in ith study
  int P; // Number of covariate interaction variables
  int G; //Number of treatments
  real prior_sd; //prior variance for params
  
  vector[nvar[1]] y1;
  matrix[nvar[1], nvar[1]] Sigmahat1;
  matrix[nvar[1], ((P+1)*(G-1))] X1;
  
  vector[nvar[2]] y2;
  matrix[nvar[2], nvar[2]] Sigmahat2;
  matrix[nvar[2],((P+1)*(G-1))] X2;
  
  vector[nvar[3]] y3;
  matrix[nvar[3], nvar[3]] Sigmahat3;
  matrix[nvar[3], ((P+1)*(G-1))] X3;
  
}

parameters {
  
  vector[((P+1)*(G-1))] phi;
  
}

transformed parameters {
  
  vector[nvar[1]] delta1;
  vector[nvar[2]] delta2;
  vector[nvar[3]] delta3;
  
  delta1 = (X1*phi);
  delta2 = (X2*phi);
  delta3 = (X3*phi);
  
}


model {
  
  y1 ~ multi_normal(delta1, Sigmahat1);
  // delta1 = X1*phi;
  
  y2 ~ multi_normal(delta2, Sigmahat2);
  // delta2 = X2*phi;
  
  y3 ~ multi_normal(delta3, Sigmahat3);
  // delta3 = X3*phi;
  
  for(i in 1:((P+1)*(G-1))) {
    
    phi[i] ~ normal(0, prior_sd);
    
  }
  
}
