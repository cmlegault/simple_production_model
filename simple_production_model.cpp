#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(obs_B);
  DATA_VECTOR(obs_Y);

  // Parameters
  PARAMETER_VECTOR(logTMB_B);
  PARAMETER_VECTOR(logitTMB_F);
  PARAMETER(logTMB_K);
  PARAMETER(logTMB_r);
  PARAMETER(logsigmaF);   // process error in F (random walk)
  PARAMETER(logsigmaB1);  // process error in B
  PARAMETER(logsigmaB2);  // observation error in B (assumes q=1 for survey)
  PARAMETER(logsigmaY);   // observation error in F
  
  // Derived quantities
  int nY = obs_B.size();  // number of observations
  Type TMB_K=exp(logTMB_K);
  Type TMB_r=exp(logTMB_r);
  Type sigmaF=exp(logsigmaF);
  Type sigmaB1=exp(logsigmaB1);
  Type sigmaB2=exp(logsigmaB2);
  Type sigmaY=exp(logsigmaY);
  vector<Type> TMB_B=exp(logTMB_B);
  vector<Type> TMB_F=invlogit(logitTMB_F);
  vector<Type> pred_Y(nY-1);
  vector<Type> pred_B(nY);
  pred_Y.setZero();
  pred_B.setZero();
  
  Type ans=0; // likelihood
  
  // compute time series 
  for(int y=1; y<nY; y++){
    pred_Y(y-1) = TMB_B(y-1) * TMB_F(y-1);
	pred_B(y) = TMB_B(y-1) + TMB_r * TMB_B(y-1) * (1.0 - TMB_B(y-1) / TMB_K) - pred_Y(y-1);
  }
	
  // process error in B
  for(int y=1; y<nY; y++){
    ans -= dnorm(log(TMB_B(y)), log(pred_B(y)), sigmaB1, true); 
  }
  
  // random walk in F
  for(int y=1; y<(nY-1); y++){
	ans -= dnorm(logitTMB_F(y), logitTMB_F(y-1), sigmaF, true);
  }
  
  // catch likelihood
  for(int y=0; y<(nY-1); y++){
    ans -= dnorm(obs_Y(y), pred_Y(y), pred_Y(y)*sigmaY, true); 
  }
  
  // likelihood for biomass observations (survey with q=1)
  for(int y=1; y<nY; y++){
    ans -= dnorm(obs_B(y), TMB_B(y), TMB_B(y)*sigmaB2, true);
  }
  
  // penalty for B0 not equal to K
  ans -= dnorm(logTMB_B(0), logTMB_K, Type(0.01), true);
  
  REPORT(pred_Y);
  REPORT(pred_B);
  
  ADREPORT(sigmaF);
  ADREPORT(sigmaB1);
  ADREPORT(sigmaB2);
  ADREPORT(sigmaY);
  ADREPORT(TMB_K);
  ADREPORT(TMB_r);
  ADREPORT(TMB_B);
  ADREPORT(TMB_F);
  
  return ans;
}
