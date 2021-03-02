// Two stage Delay Difference model as derived in Schnute 1987
//
//Model key is a density function related to number of fish within any given size range
//
// Author: LB
//Date: 26/02/2018

/// @file schnute_new_V2.hpp
// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
using namespace density;

template<class Type>
Type schnute_new_V2(objective_function<Type>* obj) {


  // data input
  DATA_INTEGER(indices_class);
  DATA_VECTOR(indices_ts);
  DATA_VECTOR(obs_catch);
  DATA_MATRIX(obs_ind);
  DATA_MATRIX(mean_wts); //a matrix dim (3,no_years) of mean weights for recruits (Y), previosuly exploited population (Z) and whole recruited population (X),
  //DATA_IVECTOR(rec_ages);
  DATA_SCALAR(nu);//fraction of total mortality that occurs before spawning
  //DATA_VECTOR(omega);
  DATA_INTEGER(SRcode);
  DATA_VECTOR(spawn_prop);



  //std::cout << "rec ages length" << std::endl << rec_ages_length << std::endl << std::endl;
  //parameters

  //timing parameters
  PARAMETER_VECTOR(logrec_param);

  PARAMETER(logB0);
  PARAMETER_VECTOR(logq);//proportionality parameter for index
  PARAMETER(logW);
  PARAMETER(logrho);
  PARAMETER_VECTOR(logf_calc);
  PARAMETER(logcatch_sigma);
  PARAMETER(logitsigma);//fraction of fish opoulation that survives natural mortality
  //PARAMETER(q_sg)//proportionality parameter for index
  //PARAMETER(rec_a); //parameter for ricker
  //PARAMETER(rec_b);//parameter for ricker
  //PARAMETER_VECTOR(rec_prop);//parameters that determine the proportion of each age that contributes to recruitment
  //PARAMETER(growth_W);
  //PARAMETER(growth_rho);
  PARAMETER_VECTOR(logindex_sigma);


  //process params
  Type sigma=invlogit(logitsigma);
  vector <Type> index_sigma = exp(logindex_sigma);
  vector <Type> qhat = exp(logq);
  vector <Type> rec_param = exp(logrec_param);

  Type B0 = exp(logB0);

  Type rho = exp(logrho);
  Type W = exp(logW);
  vector <Type> f_calc = exp(logf_calc);
  Type catch_sigma = exp(logcatch_sigma);
  //process data input
  //int rec_ages_length = rec_ages.size();
  //vector <Type> rec_prop_all(rec_ages_length);
  //if(rec_ages_length>1){

  //rec_prop_all.head(rec_ages_length-1) = rec_prop;
  //rec_prop_all.tail(1) = 1-rec_prop.sum();
  //}else{

  //rec_prop_all = rec_prop;
  //}

  //  std::cout << "rec_prop_all" << std::endl << rec_prop_all << std::endl << std::endl;

  int no_years = obs_catch.size();
  int no_ind = obs_ind.rows();

  //std::cout << "no years" << std::endl << no_years << std::endl << std::endl;
  vector <Type> omega (no_years);

  for(int i = 0; i <no_years;i++){

    omega[i] = (mean_wts(0,i)/mean_wts(2,i))*((mean_wts(1,i)-mean_wts(2,i))/(mean_wts(1,i) - mean_wts(0,i)));
  }


  vector <Type> biomass_alt(no_years);
  vector <Type> biomass(no_years);
  vector <Type> ssb(no_years);
  vector <Type> rec_bio(no_years);
  vector <Type> rec_no(no_years);
  vector <Type> post_rec(no_years);

  vector <Type> N(no_years);
  vector <Type> PR(no_years);
   vector <Type> C(no_years);


  Type X_yearon = 0.0;
  vector <Type> xprog(no_years);
  matrix <Type> logpred_survey(no_ind,no_years);
  Type theta;
  Type q;
  int tp;
  vector <Type> catch_pred(no_years);
  vector <Type> mu(no_years);
  Type m = -log(sigma);
  Type f;
  tp  = indices_class;


  f = f_calc[0];

  biomass[0] = B0;
  catch_pred[0] = f * biomass[0] * (1.0 - exp(-(f+m))) / (f+m);
  mu[0] = (f*(1 - exp(-m)) - ((m*exp(-m))*(1-exp(-f))))/ (f*(1 - exp(-m))*(1-exp(-m-f)));


  ssb[0] = spawn_prop[0]*((biomass[0]*(1-nu*(1-sigma))) - ((nu*(1 - mu[0]*(1-sigma)))*catch_pred[0]));

  theta = indices_ts[0];
  q = qhat[0];
  //biomass_alt[0] = (obs_ind(0,0) + q*((theta*(1 - mu[0]*(1-sigma)))*catch_pred[0]))/q*(1-theta*(1-sigma));
  Type rec_a = rec_param[0];
  Type rec_b = rec_param[1];


  for(int i=0;i<(no_years-1);i++){

    if(SRcode==0){//ricker SR model
      warning("This SRcode requires a index of recruit biomass only");
      rec_bio[i+1] = (rec_param[i+1]);
    }else{
      if(SRcode==1){//ricker SR model
        Type rec_a = rec_param[0];
        Type rec_b = rec_param[1];
        rec_no[i+1] = ((rec_a*ssb[i])*exp(-rec_b*ssb[i]));
        rec_bio[i+1] = rec_no[i+1]*mean_wts(0,i+1);
      }else{
        if(SRcode==2){//bevholt SR model

          rec_no[i+1] = (rec_a*ssb[i]/(rec_b +ssb[i]));
          rec_bio[i+1] = rec_no[i+1]*mean_wts(0,i+1);
        }else{
          error("SRcode not recognised");
        }
      }
    }


    X_yearon = W + rho*mean_wts(2,i);
    xprog[i] = (X_yearon/mean_wts(2,i));

      //RECRUITS
      if(tp ==1){

        biomass[i+1] = rec_bio[i+1]/omega[i+1];
        post_rec[i+1] = biomass[i+1] - rec_bio[i+1]; //xprog[i]*((sigma*(biomass[i] - (mu[i]*catch_pred[i]))) - ((1-mu[i])*catch_pred[i]));

        f = f_calc[i+1];
        catch_pred[i+1] = f * biomass[i+1] * (1.0 - exp(-(f+m))) / (f+m);
        mu[i+1] = (f*(1 - exp(-m)) - ((m*exp(-m))*(1-exp(-f))))/ (f*(1 - exp(-m))*(1-exp(-m-f)));
        ssb[i+1] = spawn_prop[i+1]*((biomass[i+1]*(1-nu*(1-sigma))) - ((nu*(1 - mu[i+1]*(1-sigma)))*catch_pred[i+1]));
        //biomass_alt[i+1] = (obs_ind(0,i+1) + q*((theta*(1 - mu[i+1]*(1-sigma)))*catch_pred[i+1]))/q*(1-theta*(1-sigma));

      }

      //POST-RECRUITS
      if(tp ==2){
        post_rec[i+1] = xprog[i]*((sigma*(biomass[i] - (mu[i]*catch_pred[i]))) - ((1-mu[i])*catch_pred[i]));

        biomass[i+1] = post_rec[i+1]/(1-omega[i+1]);
        rec_bio[i+1] = biomass[i+1] - post_rec[i+1];

        rec_no[i+1] = rec_bio[i+1] / mean_wts(0,i+1);

        f = f_calc[i+1];
        catch_pred[i+1] = f * biomass[i+1] * (1.0 - exp(-(f+m))) / (f+m);
        mu[i+1] = (f*(1 - exp(-m)) - ((m*exp(-m))*(1-exp(-f))))/ (f*(1 - exp(-m))*(1-exp(-m-f)));
        ssb[i+1] = spawn_prop[i+1]*((biomass[i+1]*(1-nu*(1-sigma))) - ((nu*(1 - mu[i+1]*(1-sigma)))*catch_pred[i+1]));

        //biomass_alt[i+1] = (obs_ind(0,i+1) + q*((theta*(1 - mu[i+1]*(1-sigma)))*catch_pred[i+1]))/q*(1-theta*(1-sigma));
      }

      //UNDIVIDED
      if(tp ==3){
        post_rec[i+1] = xprog[i]*((sigma*(biomass[i] - (mu[i]*catch_pred[i]))) - ((1-mu[i])*catch_pred[i]));

        biomass[i+1] = rec_bio[i+1] + post_rec[i+1];
        f = f_calc[i+1];
        catch_pred[i+1] = f * biomass[i+1] * (1.0 - exp(-(f+m))) / (f+m);
        mu[i+1] = (f*(1 - exp(-m)) - ((m*exp(-m))*(1-exp(-f))))/ (f*(1 - exp(-m))*(1-exp(-m-f)));
        ssb[i+1] = spawn_prop[i+1]*((biomass[i+1]*(1-nu*(1-sigma))) - ((nu*(1 - mu[i+1]*(1-sigma)))*catch_pred[i+1]));
        //biomass_alt[i+1] = (obs_ind(0,i+1) + q*((theta*(1 - mu[i+1]*(1-sigma)))*catch_pred[i+1]))/q*(1-theta*(1-sigma));

      }
    }

  rec_no[0] = rec_bio[0] / mean_wts(0,0);

  for(int i=0;i<(no_years);i++){

    N[i] = biomass[i] / mean_wts(2,i);
    PR[i] = post_rec[i] / mean_wts(1,i);
    C[i] = catch_pred[i] / mean_wts(2,i);


          for(int j=0; j<no_ind; j++){

    tp  = indices_class;
    theta = indices_ts[j];
    q = qhat[j];

    logpred_survey(j,i) = log(q*(((1-theta*(1-sigma))*biomass[i]) - ((theta*(1 - mu[i]*(1-sigma)))*catch_pred[i])));
    biomass_alt[i] = (exp(logpred_survey(0,i)) + q*((theta*(1 - mu[i]*(1-sigma)))*catch_pred[i]))/(q*(1-theta*(1-sigma)));
  }
}

  //std::cout << "post_rec" << std::endl << post_rec << std::endl << std::endl;
  //std::cout << "rec bio" << std::endl << rec_bio << std::endl << std::endl;
  //std::cout << "rec no" << std::endl << rec_no << std::endl << std::endl;
  //std::cout << "x prog" << std::endl << xprog << std::endl << std::endl;

  ////std::cout << "prev_ex" << std::endl << prev_ex << std::endl << std::endl;
  //std::cout << "biomass" << std::endl << biomass << std::endl << std::endl;
  //std::cout << "pred_survey" << std::endl << pred_survey << std::endl << std::endl;
  //std::cout << "catch" << std::endl << catch_pred << std::endl << std::endl;
  //std::cout << "omega" << std::endl << omega << std::endl << std::endl;

  //std::cout << "ssb" << std::endl << ssb << std::endl << std::endl;

  Type nll = 0.0;


  for (int i=0; i < no_years; i++){
    for (int j=0; j < no_ind; j++){


      if(obs_ind(j,i) >0){

        nll -= dnorm(log(obs_ind(j,i)), logpred_survey(j,i),index_sigma[j],true);

      }
    }
  }

  for (int i=0; i < no_years; i++){

    nll -= dnorm(log(obs_catch[i]), log(catch_pred[i]), catch_sigma, true);

  }


  vector <Type> lnb = log(biomass);
  vector <Type> lnpr = log(post_rec);
  vector <Type> lnr = log(rec_bio);
  //matrix <Type> lni = log(pred_survey);
  vector <Type> lnc = log(catch_pred);

  vector <Type> lnN = log(N);
  vector <Type> lnPR = log(PR);
  vector <Type> lnR = log(rec_no);
  vector <Type> lnC = log(C);

  ADREPORT(lnb);
  ADREPORT(lnpr);
  ADREPORT(lnr);
  ADREPORT(lnc);

  ADREPORT(lnN);
  ADREPORT(lnPR);
  ADREPORT(lnR);
  ADREPORT(lnC);

  ADREPORT(biomass);
  ADREPORT(N);
  //ADREPORT(ssb);
  ADREPORT(post_rec);
  ADREPORT(PR);
  //ADREPORT(rec);
  ADREPORT(rec_no);
  ADREPORT(rec_bio);
  ADREPORT(logpred_survey);
  ADREPORT(omega);
  ADREPORT(xprog);
  ADREPORT(catch_pred);
  ADREPORT(C);
  ADREPORT(mu);
  ADREPORT(f_calc);
  ADREPORT(biomass_alt);
  ADREPORT(sigma);
  ADREPORT(qhat);
    ADREPORT(B0);
  ADREPORT(rec_param);
  ADREPORT(catch_sigma);
  ADREPORT(index_sigma);






  return nll;

}


#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
