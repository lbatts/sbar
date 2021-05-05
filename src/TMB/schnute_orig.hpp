// Two stage Delay Difference model as derived in Schnute 1987
//
//Model key is a density function related to number of fish within any given size range
// This particular version is faithful to the Schnute description and is a process error only model with a lognormal autoregressive process
//Author LB
//Date: 26/03/2021

/// @file schnute_orig.hpp
// **DON'T** #include <TMB.hpp> as it is not include-guarded

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
//
using namespace density;
//
//
template<class Type>
Type schnute_orig(objective_function<Type>* obj) {


  // data input
  DATA_INTEGER(indices_class);
  DATA_VECTOR(indices_ts);
  DATA_VECTOR(obs_catch);
  DATA_MATRIX(obs_ind);
  DATA_MATRIX(mean_wts); //a matrix dim (3,no_years) of mean weights for recruits (Y), previously exploited population (Z) and whole recruited population (X),
  //DATA_IVECTOR(rec_ages);
  DATA_SCALAR(nu);//fraction of total mortality that occurs before spawning
  DATA_SCALAR(mu);//fraction of catch taken prior to natural mortality
  DATA_INTEGER(SRcode);
  DATA_VECTOR(spawn_prop);
  DATA_VECTOR(l_calc_wt);



  //std::cout << "rec ages length" << std::endl << rec_ages_length << std::endl << std::endl;
  //parameters

  //timing parameters
  PARAMETER_VECTOR(logrec_param);

  PARAMETER_VECTOR(logitq);//proportionality parameter for index
  PARAMETER(logW);
  PARAMETER(logrho);
  PARAMETER(logitsigma);//fraction of fish opoulation that survives natural mortality
  //PARAMETER_VECTOR(rec_prop);//parameters that determine the proportion of each age that contributes to recruitment
  PARAMETER_VECTOR(logindex_sigma);


  //process params
  Type sigma=invlogit(logitsigma);
  vector <Type> index_sigma = exp(logindex_sigma);
  vector <Type> qhat = invlogit(logitq);
  vector <Type> rec_param = exp(logrec_param);

  Type rho = exp(logrho);
  Type W = exp(logW);
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


  //vector <Type> biomass_alt(no_years);
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
  tp  = indices_class;
  matrix <Type> logB0(no_ind,no_years);
  matrix <Type> B0_ctrl(no_ind,no_years);
  vector <Type> logB0_y(no_ind);
  vector <Type> B0_ctrl2(no_ind);
  //Type B0;


  //biomass_alt[0] = (obs_ind(0,0) + q*((theta*(1 - mu[0]*(1-sigma)))*catch_pred[0]))/q*(1-theta*(1-sigma));
  Type rec_a = rec_param[0];
  Type rec_b = rec_param[1];
  tp  = indices_class;

  //std::cout << "B0" << std::endl << exp(logB0(0,0)) << std::endl << std::endl;
  //std::cout << "B0sum" << std::endl << B0.sum() << std::endl << std::endl;
  //std::cout << "B02sum" << std::endl << B0_ctrl.sum()  << std::endl << std::endl;
  for(int i=0;i<(no_years-1);i++){

    for(int j=0; j<no_ind; j++){

      theta = indices_ts[j];
      q = qhat[j];
      logB0(j,i) = 0;
      B0_ctrl(j,i) = 0;
      if(obs_ind(j,i)>0){
        logB0(j,i) = log((obs_ind(j,i) + q*((theta*(1 - mu*(1-sigma)))*obs_catch[i]))/q*(1-theta*(1-sigma)));
        B0_ctrl(j,i)= 1;
        logB0_y = logB0.col(i);
      }
    }

    B0_ctrl2 = B0_ctrl.col(i);
    biomass[i] = exp(logB0_y.sum()/ B0_ctrl2.sum());
    ssb[i] = spawn_prop[i]*((biomass[i]*(1-nu*(1-sigma))) - ((nu*(1 - mu*(1-sigma)))*obs_catch[i]));

    N[i] = biomass[i] / mean_wts(2,i);

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

    for(int j=0; j<no_ind; j++){
      theta = indices_ts[j];
      q = qhat[j];
      //RECRUITS
      if(tp ==1){


        post_rec[i+1] = (rec_bio[i+1]/omega[i+1]) - rec_bio[i+1];
        logpred_survey(j,i+1) = log(q*(((1-theta*(1-sigma))*(rec_bio[i+1]/omega[i+1])) - ((theta*(1 - mu*(1-sigma)))*obs_catch[i+1])));

      }

      //POST-RECRUITS
      if(tp ==2){
        post_rec[i+1] = xprog[i]*((sigma*(biomass[i] - (mu*obs_catch[i]))) - ((1-mu)*obs_catch[i]));

        rec_bio[i+1] = (post_rec[i+1]/(1-omega[i+1])) - post_rec[i+1];
        rec_no[i+1] = rec_bio[i+1] / mean_wts(0,i+1);

        logpred_survey(j,i+1) = log(q*(((1-theta*(1-sigma))*(post_rec[i+1]/(1-omega[i+1]))) - ((theta*(1 - mu*(1-sigma)))*obs_catch[i+1])));

      }

      //UNDIVIDED
      if(tp ==3){
        post_rec[i+1] = xprog[i]*((sigma*(biomass[i] - (mu*obs_catch[i]))) - ((1-mu)*obs_catch[i]));

        //biomass[i+1] = rec_bio[i+1] + post_rec[i+1];
        logpred_survey(j,i+1) = log(q*(((1-theta*(1-sigma))*(rec_bio[i+1] + post_rec[i+1])) - ((theta*(1 - mu*(1-sigma)))*obs_catch[i+1])));
      }

    }
    PR[i+1] = post_rec[i+1] / mean_wts(1,i+1);
    C[i+1] = obs_catch[i+1] / mean_wts(2,i+1);


  }



  for(int j=0; j<no_ind; j++){

    theta = indices_ts[j];
    q = qhat[j];
    logB0(j,(no_years-1)) = 0;
    B0_ctrl(j,(no_years-1)) = 0;
    if(obs_ind(j,(no_years-1))>0){
      logB0(j,(no_years-1)) = log((obs_ind(j,(no_years-1)) + q*((theta*(1 - mu*(1-sigma)))*obs_catch[(no_years-1)]))/q*(1-theta*(1-sigma)));
      B0_ctrl(j,(no_years-1))= 1;
      logB0_y = logB0.col((no_years-1));
    }
  }
  B0_ctrl2 = B0_ctrl.col((no_years-1));

  biomass[(no_years-1)] = exp(logB0_y.sum()/ B0_ctrl2.sum());
  ssb[(no_years-1)] = spawn_prop[(no_years-1)]*((biomass[(no_years-1)]*(1-nu*(1-sigma))) - ((nu*(1 - mu*(1-sigma)))*obs_catch[(no_years-1)]));

  N[(no_years-1)] = biomass[(no_years-1)] / mean_wts(2,(no_years-1));
  C[0] = obs_catch[0] / mean_wts(2,0);


  std::cout << "post_rec" << std::endl << post_rec << std::endl << std::endl;
  //std::cout << "rec bio" << std::endl << rec_bio << std::endl << std::endl;
  //std::cout << "rec no" << std::endl << rec_no << std::endl << std::endl;
  //std::cout << "x prog" << std::endl << xprog << std::endl << std::endl;

  ////std::cout << "prev_ex" << std::endl << prev_ex << std::endl << std::endl;
  //std::cout << "biomass" << std::endl << biomass(0) << std::endl << std::endl;
  std::cout << "pred_survey" << std::endl << logpred_survey.col(1) << std::endl << std::endl;
  //std::cout << "catch" << std::endl << obs_catch << std::endl << std::endl;
  //std::cout << "omega" << std::endl << omega << std::endl << std::endl;

  //std::cout << "ssb" << std::endl << ssb << std::endl << std::endl;

  Type nll = 0.0;


  for (int i=1; i < no_years; i++){
    for (int j=0; j < no_ind; j++){


      if(obs_ind(j,i) >0){

        nll -= l_calc_wt[j] * (dnorm(log(obs_ind(j,i)), logpred_survey(j,i),index_sigma[j],true));

      }
    }
  }


  vector <Type> lnb = log(biomass);
  vector <Type> lnpr = log(post_rec);
  vector <Type> lnr = log(rec_bio);
  //matrix <Type> lni = log(pred_survey);


  vector <Type> lnN = log(N);
  vector <Type> lnPR = log(PR);
  vector <Type> lnR = log(rec_no);
  vector <Type> lnC = log(C);

  ADREPORT(lnb);
  ADREPORT(lnpr);
  ADREPORT(lnr);


  ADREPORT(lnN);
  ADREPORT(lnPR);
  ADREPORT(lnR);
  ADREPORT(lnC);

  ADREPORT(biomass);
  ADREPORT(N);
  ADREPORT(ssb);
  ADREPORT(post_rec);
  ADREPORT(PR);
  //ADREPORT(rec);
  ADREPORT(rec_no);
  ADREPORT(rec_bio);
  ADREPORT(logpred_survey);
  ADREPORT(omega);
  ADREPORT(xprog);
  ADREPORT(C);
  ADREPORT(sigma);
  ADREPORT(qhat);
  ADREPORT(logB0);
  ADREPORT(B0_ctrl);
  ADREPORT(rec_param);
  ADREPORT(index_sigma);


  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
