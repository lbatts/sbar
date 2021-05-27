// CSA model as implemented within the NOAA Fisheries Toolbox with some slight adaptations
//
//
// Author: LB
//Date: 21/06/2019

/// @file csa.hpp
// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

using namespace density;

template<class Type>
Type csa(objective_function<Type>* obj) {

  // data input
  DATA_VECTOR(obs_catch);
  DATA_MATRIX(obs_ind);
  DATA_IMATRIX(indices_class);
  DATA_VECTOR(indices_ts);
  DATA_MATRIX(sr);


  //process data


  // parameters
  PARAMETER_VECTOR(logitqhat);
  PARAMETER(logphat1);
  PARAMETER_VECTOR(logrhat);
  PARAMETER_VECTOR(logf_calc);
  //PARAMETER_VECTOR(logc_calc);
  //if undivided survey then estimated rec selectivity is a parameters as well
  PARAMETER_VECTOR(logitsrx);
  //PARAMETER_VECTOR(logitrecsel);
  PARAMETER_VECTOR(log_surveycv);
  PARAMETER(log_catchcv);
  PARAMETER(lognmort);
  //PARAMETER_VECTOR(logit_recsel);

  //process params
  vector <Type> qhat = invlogit(logitqhat);
  Type phat1 = exp(logphat1);
  vector <Type> rhat = exp(logrhat);
  //vector <Type> c_calc = exp(logc_calc);
  vector <Type> f_calc = exp(logf_calc);
  vector <Type> srx = invlogit(logitsrx); ///fix this parameter if no undivided surveys
  Type nmort = exp(lognmort);

  //vector <Type> recsel = invlogit(logit_recsel);
  //vector <Type> recsel = exp(logit_recsel);

  vector <Type> surveycv = exp(log_surveycv);
  Type catchcv = exp(log_catchcv);


  int nyears = obs_catch.size();

  vector <Type> phat(nyears);
  phat[0] = phat1;
  Type m; Type n; Type f; Type z; Type cx; Type nx;
  vector <Type> z_calc(nyears);
  vector <Type> c_calc(nyears);

  //std::cout << "nyears" << std::endl << nyears << std::endl << std::endl;

  // Basic population dynamics
  for(int i=0;i<nyears;i++)
    {

    //obs_catch[i] = obs_land[i] + obs_disc[i];
    m = nmort;
    n = phat[i] + rhat[i];
    f = f_calc[i];
    z = f + m;
    cx = f * n * (1.0 - exp(-z)) / z;
    nx = n * exp(-z);
    z_calc[i] = z;
    c_calc[i] = cx;
    if(i < (nyears-1)){
    phat[i+1] = nx;
    }
  }



  //Calculation of Negative Log Likelihood
    Type negloglikelihood = 0.0;
  //Type nll = 0.0;

  // Calculating Catch Fit

  Type sdcatch = sqrt(log(1.0 + pow(catchcv,2))); //calculates catch sd from cv, for lognormal data
  vector <Type> resid_catch(nyears);
  Type lx;
  Type lambda_catch = 1.0;

  //std::cout << "phat" << std::endl << phat << std::endl << std::endl;
  //std::cout << "rhat" << std::endl << rhat << std::endl << std::endl;

  for (int i=0; i < nyears; i++)
  {

    if (obs_catch[i] > 0.0 && c_calc[i] > 0.0)

    {

      resid_catch[i] = log(obs_catch[i]) - log(c_calc[i]);
      lx = log(sdcatch) +  0.5 * pow((resid_catch[i] / sdcatch),2);
      SIMULATE {
        obs_catch[i] = exp(rnorm(log(c_calc[i]), sdcatch));  // Simulate response
      }

      //nll -= dnorm(log(obs_catch[i]), log(c_calc[i]),sdcatch,true);
      negloglikelihood += lx * lambda_catch;

    }
  }

  //std::cout << "nll" << std::endl << nll << std::endl << std::endl;
// Calculate residuals for each observation (ktype[j] = 1 - recruit,2 = post-recruit, 3 = undivided survey)


  int no_ind = obs_ind.rows();


  //std::cout << "no ind" << std::endl << no_ind << std::endl << std::endl;

  matrix <Type> logpred_survey(no_ind,nyears);
  matrix <Type> resid_log(no_ind,nyears);
  matrix <Type> pop_num(3,nyears);
  matrix <Type> sdsurv(no_ind,nyears);
  vector <Type> bhat(nyears);
  bhat = phat+rhat;
  pop_num.row(0) = rhat;
  pop_num.row(1) = phat;
  pop_num.row(2) = bhat;
  //sr.row(0) = srx;
  int surv;
  Type ts;
  Type cv;
  Type q;
  int tp;

  //std::cout << "pred_survey size" << std::endl << pred_survey.size() << std::endl << std::endl;
  //std::cout << "no_ind" << std::endl << no_ind << std::endl << std::endl;
  //std::cout << "nyears" << std::endl << nyears << std::endl << std::endl;


  for(int i=0; i<no_ind; i++){
    for (int j=0; j < nyears; j++){

      surv = indices_class(i,0);
      tp  = indices_class(i,1);
      ts = indices_ts[(surv-1)];
      cv = surveycv[surv-1];       ///replace with surveycv[surv-1] if you want to calculate survey specific cv rather than indices specific
      q = qhat[(surv-1)];

       //RECRUITS      recsel[(surv-1)]*
       if(tp ==1){
      logpred_survey(i,j) = log(q * pop_num(0,j) * sr((surv-1),j) * exp(-ts*z_calc[j]));
      sdsurv(i,j) = sqrt(log(1.0 + pow(cv,2)));
       }

       //POST-RECRUITS
       if(tp ==2){
      logpred_survey(i,j) = log(q * pop_num(1,j) * exp(-ts*z_calc[j]));
      sdsurv(i,j) = sqrt(log(1.0 + pow(cv,2)));
       }

       //UNDIVIDED
       if(tp ==3){
      logpred_survey(i,j) = log(q * (pop_num(0,j)*srx[i] + pop_num(1,j)) * exp(-ts*z_calc[j]));
      sdsurv(i,j) = sqrt(log(1.0 + pow(cv,2)));
       }




    }
  }


    //std::cout << "pred survey" << std::endl << pred_survey << std::endl << std::endl;


  vector <Type> lambda_survey(qhat.size(),1);//this is for survey weighting
  matrix <Type> l_calc(no_ind,nyears);
  matrix <Type> l_calc_wt(no_ind,nyears);
  //Type sdsurvtemp;

  for (int i=0; i < no_ind; i++){
    for (int j=0; j < nyears; j++){

      if(obs_ind(i,j) >0){
      resid_log(i,j) = log(obs_ind(i,j)) - logpred_survey(i,j);

        //sdsurvtemp = sqrt(log(1.0 + pow(surveycv[i],2)));

      l_calc(i,j) = log(sdsurv(i,j)) +  0.5 * pow((resid_log(i,j) / sdsurv(i,j)),2);
      
      SIMULATE {
        obs_ind(i,j) = exp(rnorm(logpred_survey(i,j), sdsurv(i,j)));  // Simulate response
      }

    l_calc_wt(i,j) = l_calc(i,j) * 1;

    //nll -= dnorm(log(indices_obs(i,j)), log(pred_survey(i,j)),sdsurv(i,j),true);

    negloglikelihood += l_calc_wt(i,j);
          }
      }
  }

  //std::cout << "sd" << std::endl << sdsurv << std::endl << std::endl;

  vector <Type> lnphat = log(phat);
  vector <Type> lnbhat = log(bhat);
  vector <Type> lnc = log(c_calc);
  SIMULATE{
    REPORT(obs_ind);          // Report the simulation
    REPORT(obs_catch);
  }
  ADREPORT(phat);
  ADREPORT(rhat);
  ADREPORT(bhat);
  ADREPORT(lnphat);
  ADREPORT(lnbhat);
  ADREPORT(lnc);
  ADREPORT(c_calc);
  ADREPORT(logpred_survey);
  ADREPORT(sdsurv);
  ADREPORT(f_calc);
  ADREPORT(phat1);
  ADREPORT(nmort)
    ADREPORT(qhat);;
  //ADREPORT(pred_survey.row(0));
  //ADREPORT(pred_survey.row(1));
  //ADREPORT(pred_survey.row(2));
  //ADREPORT(pred_survey.row(3));

  //std::cout << "nll" << std::endl << nll << std::endl << std::endl;
  return negloglikelihood;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
