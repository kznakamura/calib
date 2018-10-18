void swap(float& time1, float& timeerr1, float& time2, float& timeerr2, float& alpha, float& alphaerr, float& beta, float& betaerr){
  float time_tmp, timeerr_tmp, ab_tmp, aberr_tmp;
  time_tmp = time1;
  time1 = time2;
  time2 = time_tmp;

  timeerr_tmp = timeerr1;
  timeerr1 = timeerr2;
  timeerr2 = timeerr_tmp;

  ab_tmp = alpha;
  alpha = beta;
  beta =ab_tmp;

  aberr_tmp = alphaerr;
  alphaerr = betaerr;
  betaerr =aberr_tmp;
}


void mcr_recovery(string target, string setup, string darkdate){


  const int T_NREADCH = 6;
  const int MPPC_AMP = 5;
  const float CELL_PITCH = 50; //um unit
  const float LED_WIDTH = 1000; //ns unit
  const int SMALL_NUM = 5;
  const float V_OP = 55.0; //about condition1: VUV3=55.0V, 50um=66V, 25um=69V(setup02_25um:67.8V)
  const int SMALL_EVAL =1; //0:default, 1:error evaluation
  const TString SETUP_FILE=target+"_setup.dat";
  const TString BREAKDOWN_FILE=target+"_breakdown.dat";
  const int RECFITRANGE_MIN = 0;
  //const int READCYCLE = 0; //0: default comment out 20180301~ 
  const int CYCLE = 20;
  const float GERROR_MIN = 0.98, GERROR_MAX = 1.02;
  const string SAVEFILE = "all"; //save: org,two,modiv,twomodiv,all / NOTsave: else
  const bool PDFSAVE = 1; //0: not save pdf, 1: save canvas pdf
  const bool DRAW_CANVAS = 1; //0: not draw each canvas, 1: draw each canvas
  const string ANADATE = "_ana181010";
  
  TString SETUP_PATH="/home/kznakamura/DAQ/setup/"+SETUP_FILE;  
  TString BREAKDOWN_PATH="/home/kznakamura/DAQ/setup/"+BREAKDOWN_FILE;  
  
  
  float cell_pixel = pow(3000,2)/pow(CELL_PITCH,2);

  string date, condition, darknum;
  istringstream datestream(darkdate);
  getline(datestream,date,'_');
  getline(datestream,condition,'_');
  getline(datestream,darknum,'_');


  TString darkfile="analysed_dark_"+date+"_"+condition+"_"+darknum+".root";
  TString integfile="sum_integ_"+date+"_"+condition+".root";
  cout<<"input file: "<< darkfile <<", "<< integfile <<endl;

  TString recoveryfile="recovery_"+date+"_"+condition+".root";

  string dummysetup;
  int ch_numtmp, serial_numtmp;
  int ch_num[T_NREADCH], serial_num[T_NREADCH];
  int indexsetup=1;
  
  cout<<SETUP_PATH<<endl;
  ifstream finsetup(SETUP_PATH);
  if(!finsetup){
    cout<<"#### can not find "<<SETUP_FILE<<" ####"<<endl;
    return 1;
  }else{
    cout<<SETUP_FILE<<" was found!!"<<endl;
    while(finsetup >> dummysetup) if(dummysetup == setup) break;
    while(finsetup >> ch_numtmp >> serial_numtmp){
      ch_num[indexsetup]=ch_numtmp;
      serial_num[indexsetup]=serial_numtmp;
      cout<<setup<<": ch="<<ch_num[indexsetup]<<" serial="<<serial_num[indexsetup]<<endl;
      indexsetup++;
      if(indexsetup==T_NREADCH) break;
    }
  }

    
  string dummybreak;
  int file_tmp;
  float breakvolt_tmp;
  float deltavolt[T_NREADCH]={};
  int indexbreak=1;
  bool defaultbreakflag=0;
      
  ifstream finbreak(BREAKDOWN_PATH);
  if(!finbreak){
    cout<<"##### Can not open "<<BREAKDOWN_FILE<<" #####"<<endl;
    cout<<"##### default delta_volt(=3V) is used #####"<<endl;
    while(indexbreak<T_NREADCH){
      deltavolt[indexbreak]=3.;
      defaultbreakflag=1;
      indexbreak++;
    }
  }else{
    cout<<BREAKDOWN_FILE<<" was found!!"<<endl;
    while(finbreak >> dummybreak) if(dummybreak == setup) break;
    while(finbreak >> file_tmp >> breakvolt_tmp){
      deltavolt[indexbreak]=V_OP-breakvolt_tmp;
      if(deltavolt[indexbreak]>10){
	cout<<"#### Delta volt is too BIG ####"<<endl<<"### Something is wrong about Vop ###"<<endl;
	cout<<"ch="<<indexbreak<<" deltavolt="<<deltavolt[indexbreak]<<endl;
	//return 0;
      }else if(deltavolt[indexbreak]<0){
	cout<<"#### Delta volt is NEGATIVE ####"<<endl<<"### Something is wrong about Vop ###"<<endl;
	cout<<"ch="<<indexbreak<<" deltavolt="<<deltavolt[indexbreak]<<endl;
      }
      indexbreak++;
      if(indexbreak==T_NREADCH) break;
    }
  }

  
  
  TFile *readfile1 = new TFile(darkfile);//Form("analysed_%s",inputfile1));
  TTree *tree_hybrid = (TTree*)readfile1->Get("hybrid");
  TFile *readfile2 = new TFile(integfile);//Form("sum_%s",inputfile2));
  TTree *tree_integral = (TTree*)readfile2->Get("integral");
  TTree *tree_integralerr = (TTree*)readfile2->Get("integralerr");
    
  /*
  int cycle = (int)(tree_integral->GetEntries());
  if(READCYCLE==0){
    const int CYCLE = cycle;
  }else{
    const int CYCLE = READCYCLE;
  }
  */
  
  
  float integral_mean[T_NREADCH][CYCLE], integral_err[T_NREADCH][CYCLE];
  
  for(int cy=0; cy<CYCLE; cy++){
    tree_integral->GetEntry(cy);
    tree_integralerr->GetEntry(cy);
    for(int ch=0 ; ch<T_NREADCH ; ch++){
	integral_mean[ch][cy] = tree_integral->GetLeaf(Form("ch%d",ch))->GetValue(0);
	integral_err[ch][cy] = tree_integralerr->GetLeaf(Form("ch%d",ch))->GetValue(0);
    }
  }

  TFile *newfile = new TFile(recoveryfile,"recreate");
  TTree *tree_photonnum = new TTree("photonnum","number of photon");
  TTree *tree_reoriginal = new TTree("reoriginal","original");
  TTree *tree_retwoparam = new TTree("retwoparam","two param");
  TTree *tree_remodiv = new TTree("remodiv","modified Vop");
  TTree *tree_remoditwo = new TTree("remoditwo","modified Vop + two param");
  

  TGraphErrors *g1[T_NREADCH], *g2[T_NREADCH], *g_recovery1[T_NREADCH], *g_recovery2[T_NREADCH], *g_recovery3[T_NREADCH], *g_recovery4[T_NREADCH];
  TGraph *g_error1[T_NREADCH], *g_error2[T_NREADCH], *g_error3[T_NREADCH], *g_error4[T_NREADCH];
  TF1 *f_smallfit[T_NREADCH], *f_recoveryfit1[T_NREADCH], *f_recoveryfit2[T_NREADCH], *f_recoveryfit3[T_NREADCH], *f_recoveryfit4[T_NREADCH];
  float smallfit_slp[T_NREADCH], smallfit_cpt[T_NREADCH];
  float eff_gain[T_NREADCH]={};
  float ntrue[T_NREADCH][CYCLE]={}, nobs[T_NREADCH][CYCLE]={};
  float ntrue_error[T_NREADCH][CYCLE]={}, nobs_error[T_NREADCH][CYCLE]={};
  float recovery1_time[T_NREADCH]={},recovery1_timeerr[T_NREADCH]={};
  float recovery2_time1[T_NREADCH]={},recovery2_time2[T_NREADCH]={},recovery2_timeerr1[T_NREADCH]={},recovery2_timeerr2[T_NREADCH]={}, 
        recovery2_alpha[T_NREADCH]={}, recovery2_beta[T_NREADCH]={}, recovery2_alphaerr[T_NREADCH]={}, recovery2_betaerr[T_NREADCH]={};
  float recovery3_time[T_NREADCH]={},recovery3_timeerr[T_NREADCH]={};
  float recovery4_time1[T_NREADCH]={},recovery4_time2[T_NREADCH]={},recovery4_timeerr1[T_NREADCH]={},recovery4_timeerr2[T_NREADCH]={}, 
        recovery4_alpha[T_NREADCH]={}, recovery4_beta[T_NREADCH]={},recovery4_alphaerr[T_NREADCH]={}, recovery4_betaerr[T_NREADCH]={};
  float recovery1_diff[T_NREADCH][CYCLE]={}, recovery2_diff[T_NREADCH][CYCLE]={}, recovery3_diff[T_NREADCH][CYCLE]={}, recovery4_diff[T_NREADCH][CYCLE]={};
  float recovery1_chisq[T_NREADCH]={},recovery2_chisq[T_NREADCH]={},recovery3_chisq[T_NREADCH]={},recovery4_chisq[T_NREADCH]={};
  int b_ch;
  int b_cycle[CYCLE]={};
  float b_ntrue[CYCLE]={}, b_nobs[CYCLE]={};
  float b_ntrueerr[CYCLE]={}, b_nobserr[CYCLE]={};

  tree_photonnum -> Branch("ch", &b_ch);
  tree_photonnum -> Branch("cycle", b_cycle, Form("b_cycle[%d]/I",CYCLE));
  tree_photonnum -> Branch("ntrue", b_ntrue, Form("b_ntrue[%d]/F",CYCLE));
  tree_photonnum -> Branch("nobs", b_nobs, Form("b_nobs[%d]/F",CYCLE));
  tree_photonnum -> Branch("ntrueerr", b_ntrueerr, Form("b_ntrueerr[%d]/F",CYCLE));
  tree_photonnum -> Branch("nobserr", b_nobserr, Form("b_nobserr[%d]/F",CYCLE));

  int channel[T_NREADCH]={};
  for(int ch=0; ch<T_NREADCH; ch++){
    channel[ch]=ch;
  }

  tree_reoriginal -> Branch("ch", channel, Form("channel[%d]/I",T_NREADCH));
  tree_reoriginal -> Branch("rectime", recovery1_time, Form("recovery1_time[%d]/F",T_NREADCH));
  tree_reoriginal -> Branch("recerror", recovery1_timeerr, Form("recovery1_timeerr[%d]/F",T_NREADCH));
  tree_reoriginal -> Branch("chisq", recovery1_chisq, Form("recovery1_chisq[%d]/F",T_NREADCH));

  tree_retwoparam -> Branch("ch", channel, Form("channel[%d]/I",T_NREADCH));
  tree_retwoparam -> Branch("alpha", recovery2_alpha, Form("recovery2_alpha[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("alphaerr", recovery2_alphaerr, Form("recovery2_alphaerr[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("rectime1", recovery2_time1, Form("recovery2_time1[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("recerror1", recovery2_timeerr1, Form("recovery2_timeerr1[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("beta", recovery2_beta, Form("recovery2_beta[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("betaerr", recovery2_betaerr, Form("recovery2_betaerr[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("rectime2", recovery2_time2, Form("recovery2_time2[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("recerror2", recovery2_timeerr2, Form("recovery2_timeerr2[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("chisq", recovery2_chisq, Form("recovery2_chisq[%d]/F",T_NREADCH));

  tree_remodiv -> Branch("ch", channel, Form("channel[%d]/I",T_NREADCH));
  tree_remodiv -> Branch("rectime", recovery3_time, Form("recovery3_time[%d]/F",T_NREADCH));
  tree_remodiv -> Branch("recerror", recovery3_timeerr, Form("recovery3_timeerr[%d]/F",T_NREADCH));
  tree_remodiv -> Branch("chisq", recovery3_chisq, Form("recovery3_chisq[%d]/F",T_NREADCH));
  
  tree_remoditwo -> Branch("ch", channel, Form("channel[%d]/I",T_NREADCH));
  tree_remoditwo -> Branch("alpha", recovery4_alpha, Form("recovery4_alpha[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("alphaerr", recovery4_alphaerr, Form("recovery4_alphaerr[%d]/F",T_NREADCH));
  tree_remoditwo -> Branch("rectime1", recovery4_time1, Form("recovery4_time1[%d]/F",T_NREADCH));
  tree_remoditwo -> Branch("recerror1", recovery4_timeerr1, Form("recovery4_timeerr1[%d]/F",T_NREADCH));
  tree_remoditwo -> Branch("beta", recovery4_beta, Form("recovery4_beta[%d]/F",T_NREADCH));
  tree_retwoparam -> Branch("betaerr", recovery4_betaerr, Form("recovery4_betaerr[%d]/F",T_NREADCH));
  tree_remoditwo -> Branch("rectime2", recovery4_time2, Form("recovery4_time2[%d]/F",T_NREADCH));
  tree_remoditwo -> Branch("recerror2", recovery4_timeerr2, Form("recovery4_timeerr2[%d]/F",T_NREADCH));
  tree_remoditwo -> Branch("chisq", recovery4_chisq, Form("recovery4_chisq[%d]/F",T_NREADCH));

  //error calculation
  float integral_errsum[T_NREADCH][SMALL_NUM]={};
  float integral_delta[T_NREADCH]={};
  float integral_slp[T_NREADCH]={}, integral_cpt[T_NREADCH]={};
  float integral_slperr[T_NREADCH]={}, integral_cpterr[T_NREADCH]={};
  float sum_w[T_NREADCH]={}, sum_wx2[T_NREADCH]={}, sum_wx[T_NREADCH]={}, sum_wy[T_NREADCH]={}, sum_wxy[T_NREADCH]={};

  for(int ch=1; ch<T_NREADCH; ch++){
    for(int cy=0; cy<SMALL_NUM; cy++){
      integral_errsum[ch][cy]=1/(pow(integral_err[ch][cy],2)+pow(integral_err[0][cy],2));
    }
    for(int cy=0; cy<SMALL_NUM; cy++){
      sum_w[ch]+=integral_errsum[ch][cy];
      sum_wx2[ch]+=integral_errsum[ch][cy]*pow(integral_mean[0][cy],2);
      sum_wx[ch]+=integral_errsum[ch][cy]*integral_mean[0][cy];
      sum_wy[ch]+=integral_errsum[ch][cy]*integral_mean[ch][cy];
      sum_wxy[ch]+=integral_errsum[ch][cy]*integral_mean[0][cy]*integral_mean[ch][cy];      
    }
    integral_delta[ch]=sum_w[ch]*sum_wx2[ch]-pow(sum_wx[ch],2);
    integral_slp[ch]=(sum_w[ch]*sum_wxy[ch]-sum_wx[ch]*sum_wy[ch])/integral_delta[ch];
    integral_cpt[ch]=(sum_wx2[ch]*sum_wy[ch]-sum_wx[ch]*sum_wxy[ch])/integral_delta[ch];
    integral_slperr[ch]=sqrt(sum_w[ch]/integral_delta[ch]);
    integral_cpterr[ch]=sqrt(sum_wx2[ch]/integral_delta[ch]);
    //    cout<<"ch="<<ch<<" katamuki="<<integral_slp[ch]<<"+/-"<<integral_slperr[ch]<<endl;
  }
  
    
  for(int ch=1; ch<T_NREADCH; ch++){

    // tree_integral->Draw(Form("ch%d:ch0",ch),"ch0<700","goff");
    //g1[ch] = new TGraph(tree_integral->GetSelectedRows(),tree_integral->GetV2(),tree_integral->GetV1());
    g1[ch] = new TGraphErrors(SMALL_NUM,integral_mean[0],integral_mean[ch],integral_err[0],integral_err[ch]);
    f_smallfit[ch] = new TF1(Form("f_smallfit_%d",ch),"[0]*x+[1]", 0, 3000);

    switch(SMALL_EVAL){
    case 0:
    g1[ch]->Fit(Form("f_smallfit_%d",ch),"Q","",0,integral_mean[ch][SMALL_NUM-1]);
    smallfit_slp[ch] = f_smallfit[ch] -> GetParameter(0);
    smallfit_cpt[ch] = f_smallfit[ch] -> GetParameter(1);
    // f[ch] -> SetParameters(0,100);
    tree_hybrid -> GetEntry(0);
    eff_gain[ch] = tree_hybrid -> GetLeaf("hybrid_effgain") -> GetValue(ch);
    for(int cy=0; cy<CYCLE; cy++){
      //tree_integral -> GetEntry(cy);
      ntrue[ch][cy] = integral_mean[0][cy]*smallfit_slp[ch]*0.2/(MPPC_AMP*eff_gain[ch]*1.602)*pow(10,7);
      nobs[ch][cy]  = integral_mean[ch][cy]*0.2/(MPPC_AMP*eff_gain[ch]*1.602)*pow(10,7);
      //tree_integralerr -> GetEntry(cy);
      ntrue_error[ch][cy] = integral_err[0][cy]*smallfit_slp[ch]*0.2/(MPPC_AMP*eff_gain[ch]*1.602)*pow(10,7);
      nobs_error[ch][cy]  = integral_err[ch][cy]*0.2/(MPPC_AMP*eff_gain[ch]*1.602)*pow(10,7);
    }
    break;
    case 1:
      f_smallfit[ch]->SetParameters(integral_slp[ch],integral_cpt[ch]);
      tree_hybrid -> GetEntry(0);
      eff_gain[ch] = tree_hybrid -> GetLeaf("hybrid_effgain") -> GetValue(ch);
      //eff_gain[ch] = tree_hybrid -> GetLeaf("hybrid_opgain") -> GetValue(ch);  //test of one photon gain 170107
    for(int cy=0; cy<CYCLE; cy++){
      // tree_integral -> GetEntry(cy);
      ntrue[ch][cy] = integral_mean[0][cy]*integral_slp[ch]*0.2/(MPPC_AMP*eff_gain[ch]*1.602)*pow(10,7);
      nobs[ch][cy]  = integral_mean[ch][cy]*0.2/(MPPC_AMP*eff_gain[ch]*1.602)*pow(10,7);
      //     tree_integralerr -> GetEntry(cy);
      ntrue_error[ch][cy] = integral_mean[0][cy]*integral_slperr[ch]*0.2/(MPPC_AMP*eff_gain[ch]*1.602)*pow(10,7);
      nobs_error[ch][cy]  = 0;//integral_err[ch][cy]*0.2/(MPPC_AMP*eff_gain[ch]*1.602)*pow(10,7);
    } 
      break;
    }
    cout << "\n### channel="<<ch<<", eff_gain="<<eff_gain[ch]<<", deltavolt="<<deltavolt[ch]<<"V ###"<<endl;  
      
    g2[ch] = new TGraphErrors(SMALL_NUM,ntrue[ch],nobs[ch],ntrue_error[ch],nobs_error[ch]);
    g_recovery1[ch] = new TGraphErrors(CYCLE,ntrue[ch],nobs[ch],ntrue_error[ch],nobs_error[ch]);
    g_recovery2[ch] = new TGraphErrors(CYCLE,ntrue[ch],nobs[ch],ntrue_error[ch],nobs_error[ch]);
    g_recovery3[ch] = new TGraphErrors(CYCLE,ntrue[ch],nobs[ch],ntrue_error[ch],nobs_error[ch]);
    g_recovery4[ch] = new TGraphErrors(CYCLE,ntrue[ch],nobs[ch],ntrue_error[ch],nobs_error[ch]);
    //g_recovery1[ch] -> SetMarkerStyle(20);
    //    g_recovery1[ch] ->Fit("pol1","","",0,1000);
    f_recoveryfit1[ch] = new TF1(Form("f_recoveryfit1_%d",ch),"x/(1+[0]*x)", 0, 30000);
    f_recoveryfit1[ch] -> SetParameter(0,120/cell_pixel/LED_WIDTH);
    g_recovery1[ch] -> Fit(Form("f_recoveryfit1_%d",ch),"QW","",RECFITRANGE_MIN,ntrue[ch][CYCLE-1]);
    recovery1_time[ch] = (f_recoveryfit1[ch] -> GetParameter(0))*cell_pixel*LED_WIDTH;
    recovery1_timeerr[ch] = (f_recoveryfit1[ch] -> GetParError(0))*cell_pixel*LED_WIDTH;
    recovery1_chisq[ch] = f_recoveryfit1[ch] -> GetChisquare();
    cout<<left<<setw(22)<<"original"<<"tau="<<setprecision(3)<<recovery1_time[ch] <<"ns"<<endl;


    f_recoveryfit2[ch] = new TF1(Form("f_recoveryfit2_%d",ch),"[2]*x/(1+[0]*x)+[3]*x/(1+[1]*x)", 0, 30000);
    //f_recoveryfit2[ch] = new TF1(Form("f_recoveryfit2_%d",ch),"[2]*x/(1+[0]*x)+(1-[2])*x/(1+[1]*x)", 0, 30000);
    f_recoveryfit2[ch] -> SetParameters(50/cell_pixel/LED_WIDTH,100/cell_pixel/LED_WIDTH,0.5,0.5);
    //f_recoveryfit2[ch] -> SetParameters(50/cell_pixel/LED_WIDTH,100/cell_pixel/LED_WIDTH,0.5);
    f_recoveryfit2[ch] -> SetParLimits(0,0,10000/cell_pixel/LED_WIDTH);
    f_recoveryfit2[ch] -> SetParLimits(1,0,10000/cell_pixel/LED_WIDTH);
    f_recoveryfit2[ch] -> SetParLimits(2,0,1);
    g_recovery2[ch] -> Fit(Form("f_recoveryfit2_%d",ch),"QW","",RECFITRANGE_MIN,ntrue[ch][CYCLE-1]);
    recovery2_time1[ch] = (f_recoveryfit2[ch] -> GetParameter(0))*cell_pixel*LED_WIDTH;
    recovery2_timeerr1[ch] = (f_recoveryfit2[ch] -> GetParError(0))*cell_pixel*LED_WIDTH;
    recovery2_time2[ch] = (f_recoveryfit2[ch] -> GetParameter(1))*cell_pixel*LED_WIDTH;
    recovery2_timeerr2[ch] = (f_recoveryfit2[ch] -> GetParError(1))*cell_pixel*LED_WIDTH;
    recovery2_alpha[ch] = f_recoveryfit2[ch]->GetParameter(2);
    recovery2_alphaerr[ch] = f_recoveryfit2[ch]->GetParError(2);
    recovery2_beta[ch] =  f_recoveryfit2[ch]->GetParameter(3);
    // recovery2_beta[ch] =  1-recovery2_alpha[ch];
    recovery2_betaerr[ch] = f_recoveryfit2[ch]->GetParError(3);
    recovery2_chisq[ch] = f_recoveryfit2[ch] -> GetChisquare();
    //recovery2_betaerr[ch] = f_recoveryfit2[ch]->GetParError(2);
    if(recovery2_time1[ch]>recovery2_time2[ch]){
      swap(recovery2_time1[ch], recovery2_timeerr1[ch], recovery2_time2[ch], recovery2_timeerr2[ch], recovery2_alpha[ch], recovery2_alphaerr[ch], recovery2_beta[ch], recovery2_betaerr[ch]);
      }

    cout<<left<<setw(22)<<"2param"<<"alpha="<<setprecision(2)<<recovery2_alpha[ch]<<" tau1="<<setprecision(3)<<recovery2_time1[ch] <<"ns, beta="<<setprecision(2)<<recovery2_beta[ch]<<" tau2="<<setprecision(3)<<recovery2_time2[ch] <<"ns"<<endl;

    //f_recoveryfit3[ch] = new TF1(Form("f_recoveryfit3_%d",ch),"x/((1+[0]*x)*(1+[1]*x))", 0, 30000); old modiVop ~20171213
    f_recoveryfit3[ch] = new TF1(Form("f_recoveryfit3_%d",ch),"x/((1+[0]*x)*(1+[1]*x/(1+[0]*x)))", 0, 30000);
    f_recoveryfit3[ch] -> FixParameter(1,50*eff_gain[ch]*1.602*pow(10,-10)/(deltavolt[ch]*LED_WIDTH));
    f_recoveryfit3[ch] -> SetParameter(0,120/cell_pixel/LED_WIDTH);
    //    f_recoveryfit3[ch] -> SetParameter(0,120/cell_pixel/LED_WIDTH);
    g_recovery3[ch] -> Fit(Form("f_recoveryfit3_%d",ch),"QW","",RECFITRANGE_MIN,ntrue[ch][CYCLE-1]);
    recovery3_time[ch] = (f_recoveryfit3[ch] -> GetParameter(0))*cell_pixel*LED_WIDTH;
    recovery3_timeerr[ch] = (f_recoveryfit3[ch] -> GetParError(0))*cell_pixel*LED_WIDTH;
    recovery3_chisq[ch] = f_recoveryfit3[ch] -> GetChisquare();
    cout<<left<<setw(22)<<"modified Vop"<<"tau="<<setprecision(3)<<recovery3_time[ch] <<"ns"<<endl;

    //    f_recoveryfit4[ch] = new TF1(Form("f_recoveryfit4_%d",ch),"[2]*x/((1+[0]*x)*(1+[4]*x))+[3]*x/((1+[1]*x)*(1+[5]*x))", 0, 30000); old modiVop ~20171213
    f_recoveryfit4[ch] = new TF1(Form("f_recoveryfit4_%d",ch),"[2]*x/((1+[0]*x)*(1+[4]*x/(1+[0]*x)))+[3]*x/((1+[1]*x)*(1+[5]*x/(1+[1]*x)))", 0, 30000); //alpha+beta =! 1
    //f_recoveryfit4[ch] = new TF1(Form("f_recoveryfit4_%d",ch),"[2]*x/((1+[0]*x)*(1+[4]*x/(1+[0]*x)))+(1-[2])*x/((1+[1]*x)*(1+[5]*x/(1+[1]*x)))", 0, 30000); //alpha+beta = 1
    f_recoveryfit4[ch] -> SetParameters(50/cell_pixel/LED_WIDTH,100/cell_pixel/LED_WIDTH,0.5,0.5);
    //f_recoveryfit4[ch] -> SetParameters(50/cell_pixel/LED_WIDTH,100/cell_pixel/LED_WIDTH,0.5);
    //f_recoveryfit4[ch] -> SetParLimits(0,5/cell_pixel/LED_WIDTH,100/cell_pixel/LED_WIDTH);
    f_recoveryfit4[ch] -> SetParLimits(0,0,10000/cell_pixel/LED_WIDTH);
    f_recoveryfit4[ch] -> SetParLimits(1,0,10000/cell_pixel/LED_WIDTH);
    f_recoveryfit4[ch] -> SetParLimits(2,0,1);
    f_recoveryfit4[ch] -> FixParameter(4,50*eff_gain[ch]*1.602*pow(10,-10)/(deltavolt[ch]*LED_WIDTH));
    f_recoveryfit4[ch] -> FixParameter(5,50*eff_gain[ch]*1.602*pow(10,-10)/(deltavolt[ch]*LED_WIDTH));
   
    g_recovery4[ch] -> Fit(Form("f_recoveryfit4_%d",ch),"QW","",RECFITRANGE_MIN,ntrue[ch][CYCLE-1]);
    recovery4_time1[ch] = (f_recoveryfit4[ch] -> GetParameter(0))*cell_pixel*LED_WIDTH;
    recovery4_timeerr1[ch] = (f_recoveryfit4[ch] -> GetParError(0))*cell_pixel*LED_WIDTH;
    recovery4_time2[ch] = (f_recoveryfit4[ch] -> GetParameter(1))*cell_pixel*LED_WIDTH;
    recovery4_timeerr2[ch] = (f_recoveryfit4[ch] -> GetParError(1))*cell_pixel*LED_WIDTH;
    recovery4_alpha[ch] = f_recoveryfit4[ch]->GetParameter(2);
    recovery4_alphaerr[ch] = f_recoveryfit4[ch]->GetParError(2);
    recovery4_beta[ch] =  f_recoveryfit4[ch]->GetParameter(3);
    recovery4_betaerr[ch] = f_recoveryfit4[ch]->GetParError(3);
    recovery4_chisq[ch] = f_recoveryfit4[ch] -> GetChisquare();
    //recovery4_beta[ch] =  1-recovery4_alpha[ch];
    //recovery4_betaerr[ch] = f_recoveryfit4[ch]->GetParError(2);

    if(recovery4_time1[ch]>recovery4_time2[ch]){
      swap(recovery4_time1[ch], recovery4_timeerr1[ch], recovery4_time2[ch], recovery4_timeerr2[ch], recovery4_alpha[ch], recovery4_alphaerr[ch], recovery4_beta[ch], recovery4_betaerr[ch]);
    }
    cout<<left<<setw(22)<<"2param + modified Vop"<<"alpha="<<setprecision(2)<<recovery4_alpha[ch]<<" tau1="<<setprecision(3)<<recovery4_time1[ch] <<"ns, beta="<<setprecision(2)<<recovery4_beta[ch]<<" tau2="<<setprecision(3)<<recovery4_time2[ch] <<"ns"<<endl;    

    for(int cy=0; cy<CYCLE; cy++){
      recovery1_diff[ch][cy]=nobs[ch][cy]/f_recoveryfit1[ch]->Eval(ntrue[ch][cy]);
      recovery2_diff[ch][cy]=nobs[ch][cy]/f_recoveryfit2[ch]->Eval(ntrue[ch][cy]);
      recovery3_diff[ch][cy]=nobs[ch][cy]/f_recoveryfit3[ch]->Eval(ntrue[ch][cy]);
      recovery4_diff[ch][cy]=nobs[ch][cy]/f_recoveryfit4[ch]->Eval(ntrue[ch][cy]);
    }
    g_error1[ch] = new TGraph(CYCLE,ntrue[ch],recovery1_diff[ch]);
    g_error2[ch] = new TGraph(CYCLE,ntrue[ch],recovery2_diff[ch]);
    g_error3[ch] = new TGraph(CYCLE,ntrue[ch],recovery3_diff[ch]);
    g_error4[ch] = new TGraph(CYCLE,ntrue[ch],recovery4_diff[ch]);
  }

  for(int ch=1; ch<T_NREADCH; ch++){
    b_ch=ch;
    for(int cy=0; cy<CYCLE; cy++){
      b_cycle[cy]=cy;
      b_ntrue[cy]=ntrue[ch][cy];
      b_nobs[cy]=nobs[ch][cy];
      b_ntrueerr[cy]=ntrue_error[ch][cy];
      b_nobserr[cy]=nobs_error[ch][cy];
    }
    tree_photonnum->Fill();
  }
  tree_reoriginal->Fill();
  tree_retwoparam->Fill();
  tree_remodiv->Fill();
  tree_remoditwo->Fill();

  tree_photonnum->Write();
  tree_reoriginal->Write();
  tree_retwoparam->Write();
  tree_remodiv->Write();
  tree_remoditwo->Write();
  
  if(SAVEFILE=="all"){
    string filetmp;
    for(int i=0; i<4; i++){
      if(i==0) filetmp="org";
      if(i==1) filetmp="two";
      if(i==2) filetmp="modiv";
      if(i==3) filetmp="twomodiv";

      TString RECOVERY_FILE=target+"_recovery_cnd"+condition+"_"+filetmp+ANADATE+".dat";
      TString RECOVERY_PATH="/home/kznakamura/DAQ/analyzed/"+RECOVERY_FILE;  

	ofstream foutrecovery(RECOVERY_PATH, ios_base::app);
	if(!foutrecovery){
	  cout<<"### Can not open "<<RECOVERY_FILE<<" ###"<<endl;
	  return 1;
	}
	cout<<endl<<"### Analized data is filled in "<< RECOVERY_FILE<<" ###"<<endl;
	foutrecovery<<endl;
	foutrecovery<<"### "<<setup<<" ###"<<endl;
	if(defaultbreakflag){
	  foutrecovery<<"### default delta_volt(=3V) is used for modified Vop analysys ###"<<endl;
	}
	if(filetmp=="org"){
	  cout<<"### analized with original recovery time ###"<<endl;
	  foutrecovery<<"### analized with original recovery time ###"<<endl;
	}else if(filetmp=="two"){
	  cout<<"### analized with two param recovery time ###"<<endl;
	  foutrecovery<<"### analized with two param recovery time ###"<<endl;
	}else if(filetmp=="modiv"){
	  cout<<"### analized with modified Vop recovery time ###"<<endl;
	  foutrecovery<<"### analized with modified Vop recovery time ###"<<endl;
	}else if(filetmp=="twomodiv"){
	  cout<<"###analized with two param + modified Vop recovery time ###"<<endl;
	  foutrecovery<<"### analized with two param + modified Vop recovery time ###"<<endl;
	}
	

	string targetch=target+"ch";
	if(filetmp=="org" || filetmp=="modiv"){
	  foutrecovery<<left<<setw(15)<<targetch<<setw(15)<<"Serial"<<setw(15)<<"effgain"<<setw(15)<<"tau(ns)"<<setw(15)<<"tauerr(ns)"<<setw(15)<<"chisq"<<endl;
	  if(filetmp=="org"){
	    for(int ch=1; ch<T_NREADCH; ch++){ 
	      foutrecovery<<left<<setw(15)<<ch_num[ch]<<setw(15)<<serial_num[ch]<<setw(15)<<eff_gain[ch]<<setw(15)<<recovery1_time[ch]<<setw(15)<<recovery1_timeerr[ch]<<setw(15)<<recovery1_chisq[ch]<<endl;
	    }
	  }else if(filetmp=="modiv"){
	    for(int ch=1; ch<T_NREADCH; ch++){ 
	      foutrecovery<<left<<setw(15)<<ch_num[ch]<<setw(15)<<serial_num[ch]<<setw(15)<<eff_gain[ch]<<setw(15)<<recovery3_time[ch]<<setw(15)<<recovery3_timeerr[ch]<<setw(15)<<recovery3_chisq[ch]<<endl;
	    }
	  }
	}else if(filetmp=="two" || filetmp=="twomodiv"){
	  foutrecovery<<left<<setw(15)<<targetch<<setw(15)<<"Serial"<<setw(15)<<"effgain"<<setw(15)<<"alpha"<<setw(15)<<"alphaerr"<<setw(15)<<"beta"<<setw(15)<<"betaerr"<<setw(15)<<"tau1(ns)"<<setw(15)<<"tau1err(ns)"<<setw(15)<<"tau2(ns)"<<setw(15)<<"tau2err(ns)"<<setw(15)<<"chisq"<<endl;
      if(filetmp=="two"){
	for(int ch=1; ch<T_NREADCH; ch++){ 
	  foutrecovery<<left<<setw(15)<<ch_num[ch]<<setw(15)<<serial_num[ch]<<setw(15)<<eff_gain[ch]<<setw(15)<<recovery2_alpha[ch]<<setw(15)<<recovery2_alphaerr[ch]<<setw(15)<<recovery2_beta[ch]<<setw(15)<<recovery2_betaerr[ch]<<setw(15)<<recovery2_time1[ch]<<setw(15)<<recovery2_timeerr2[ch]<<setw(15)<<recovery2_time2[ch]<<setw(15)<<recovery2_timeerr2[ch]<<setw(15)<<recovery2_chisq[ch]<<endl;
	}
      }else if(filetmp=="twomodiv"){
	for(int ch=1; ch<T_NREADCH; ch++){ 
	  foutrecovery<<left<<setw(15)<<ch_num[ch]<<setw(15)<<serial_num[ch]<<setw(15)<<eff_gain[ch]<<setw(15)<<recovery4_alpha[ch]<<setw(15)<<recovery4_alphaerr[ch]<<setw(15)<<recovery4_beta[ch]<<setw(15)<<recovery4_betaerr[ch]<<setw(15)<<recovery4_time1[ch]<<setw(15)<<recovery4_timeerr2[ch]<<setw(15)<<recovery4_time2[ch]<<setw(15)<<recovery4_timeerr2[ch]<<setw(15)<<recovery4_chisq[ch]<<endl;
	}
      }
	}	
	foutrecovery.close();
    }
  }else if(SAVEFILE=="org" || SAVEFILE=="two" || SAVEFILE=="modiv" || SAVEFILE=="twomodiv"){

      TString RECOVERY_FILE=target+"_recovery_cnd"+condition+"_"+SAVEFILE+ANADATE+".dat";
      TString RECOVERY_PATH="/home/kznakamura/DAQ/analyzed/"+RECOVERY_FILE;  

      ofstream foutrecovery(RECOVERY_PATH, ios_base::app);
      if(!foutrecovery){
	cout<<"### Can not open "<<RECOVERY_FILE<<" ###"<<endl;
	return 1;
      }
	cout<<endl<<"### Analized data is filled in "<< RECOVERY_FILE<<" ###"<<endl;
	foutrecovery<<endl;
	foutrecovery<<"### "<<setup<<" ###"<<endl;
	if(defaultbreakflag){
	  foutrecovery<<"### default delta_volt(=3V) is used for modified Vop analysys ###"<<endl;
	}
	if(SAVEFILE=="org"){
	  cout<<"### analized with original recovery time ###"<<endl;
	  foutrecovery<<"### analized with original recovery time ###"<<endl;
	}else if(SAVEFILE=="two"){
	  cout<<"### analized with two param recovery time ###"<<endl;
	  foutrecovery<<"### analized with two param recovery time ###"<<endl;
	}else if(SAVEFILE=="modiv"){
	  cout<<"### analized with modified Vop recovery time ###"<<endl;
	  foutrecovery<<"### analized with modified Vop recovery time ###"<<endl;
	}else if(SAVEFILE=="twomodiv"){
	  cout<<"###analized with two param + modified Vop recovery time ###"<<endl;
	  foutrecovery<<"### analized with two param + modified Vop recovery time ###"<<endl;
	}
	
	string targetch=target+"ch";
	if(SAVEFILE=="org" || SAVEFILE=="modiv"){
	  foutrecovery<<left<<setw(15)<<targetch<<setw(15)<<"Serial"<<setw(15)<<"effgain"<<setw(15)<<"tau(ns)"<<setw(15)<<"tauerr(ns)"<<setw(15)<<"chisq"<<endl;
	  if(SAVEFILE=="org"){
	    for(int ch=1; ch<T_NREADCH; ch++){ 
	      foutrecovery<<left<<setw(15)<<ch_num[ch]<<setw(15)<<serial_num[ch]<<setw(15)<<eff_gain[ch]<<setw(15)<<recovery1_time[ch]<<setw(15)<<recovery1_timeerr[ch]<<setw(15)<<recovery1_chisq[ch]<<endl;
	    }
	  }else if(SAVEFILE=="modiv"){
	    for(int ch=1; ch<T_NREADCH; ch++){ 
	      foutrecovery<<left<<setw(15)<<ch_num[ch]<<setw(15)<<serial_num[ch]<<setw(15)<<eff_gain[ch]<<setw(15)<<recovery3_time[ch]<<setw(15)<<recovery3_timeerr[ch]<<setw(15)<<recovery3_chisq[ch]<<endl;
	    }
	  }
	}else if(SAVEFILE=="two" || SAVEFILE=="twomodiv"){
	  foutrecovery<<left<<setw(15)<<targetch<<setw(15)<<"Serial"<<setw(15)<<"effgain"<<setw(15)<<"alpha"<<setw(15)<<"alphaerr"<<setw(15)<<"beta"<<setw(15)<<"betaerr"<<setw(15)<<"tau1(ns)"<<setw(15)<<"tau1err(ns)"<<setw(15)<<"tau2(ns)"<<setw(15)<<"tau2err(ns)"<<setw(15)<<"chisq"<<endl;
	  if(SAVEFILE=="two"){
	    for(int ch=1; ch<T_NREADCH; ch++){ 
	      foutrecovery<<left<<setw(15)<<ch_num[ch]<<setw(15)<<serial_num[ch]<<setw(15)<<eff_gain[ch]<<setw(15)<<recovery2_alpha[ch]<<setw(15)<<recovery2_alphaerr[ch]<<setw(15)<<recovery2_beta[ch]<<setw(15)<<recovery2_betaerr[ch]<<setw(15)<<recovery2_time1[ch]<<setw(15)<<recovery2_timeerr2[ch]<<setw(15)<<recovery2_time2[ch]<<setw(15)<<recovery2_timeerr2[ch]<<setw(15)<<recovery2_chisq[ch]<<endl;
	    }
	  }else if(SAVEFILE=="twomodiv"){
	    for(int ch=1; ch<T_NREADCH; ch++){ 
	      foutrecovery<<left<<setw(15)<<ch_num[ch]<<setw(15)<<serial_num[ch]<<setw(15)<<eff_gain[ch]<<setw(15)<<recovery4_alpha[ch]<<setw(15)<<recovery4_alphaerr[ch]<<setw(15)<<recovery4_beta[ch]<<setw(15)<<recovery4_betaerr[ch]<<setw(15)<<recovery4_time1[ch]<<setw(15)<<recovery4_timeerr2[ch]<<setw(15)<<recovery4_time2[ch]<<setw(15)<<recovery4_timeerr2[ch]<<setw(15)<<recovery4_chisq[ch]<<endl;
	    }
	  }
	}
	
	foutrecovery.close();
  }else{
    cout<<endl<<"### Analized data is NOT saved ###"<<endl;
    } 
  
     
      //Cnavas
  TCanvas *c1[T_NREADCH], *c2[T_NREADCH];
  TLine *l[T_NREADCH];

  for(int ch=1; ch<T_NREADCH; ch++){
    c1[ch] = new TCanvas(Form("c1_ch%d", ch),Form("c1_ch%d",ch),200*(ch-1),0,800,400);
    c1[ch]->Divide(2,1,0.01,0.01);
    c1[ch]->cd(1);  
    g1[ch] ->SetMarkerStyle(20);
    g1[ch]->SetTitle(Form("small fit (serial:%d);PMT ADC count;MPPC ADC count",serial_num[ch]));
    g1[ch] -> SetMinimum(0);
    g1[ch]->GetYaxis()->SetTitleOffset(1.5);
    g1[ch] -> Draw("ap");
    if(SMALL_EVAL==1){
      f_smallfit[ch]->Draw("same");
    }

    c1[ch] -> cd(2);
    g2[ch] ->SetMarkerStyle(20);
    g2[ch]->SetTitle(Form("small recovery original (serial:%d);Nref (#mus^{-1});Nobs (#mus^{-1})",serial_num[ch]));
    g2[ch] -> SetMinimum(0);
    g2[ch]->GetYaxis()->SetTitleOffset(1.2);
    g2[ch] -> Draw("ap");
    f_recoveryfit1[ch] -> Draw("same");

    
    c2[ch] = new TCanvas(Form("c2_ch%d", ch), Form("c2_ch%d",ch),200*(ch-1),200,1600,600);
    c2[ch]->Divide(4,2,0.01,0.01);
    c2[ch]->cd(1);
    g_recovery1[ch] -> SetMaximum(nobs[ch][CYCLE-1]+2000);
    g_recovery1[ch] -> SetMarkerStyle(20);
    g_recovery1[ch] -> GetXaxis() -> SetRangeUser(0,ntrue[ch][CYCLE-1]+5000);
    g_recovery1[ch]->SetTitle(Form("recovery original (serial:%d);Nref (#mus^{-1});Nobs (#mus^{-1})",serial_num[ch]));
    g_recovery1[ch]->GetYaxis()->SetTitleOffset(1.05);
    g_recovery1[ch]->GetYaxis()->SetTitleSize(0.05);
    g_recovery1[ch]->GetXaxis()->SetTitleSize(0.05);
    g_recovery1[ch] -> Draw("ap");
    c2[ch]->cd(2);
    g_recovery2[ch] -> SetMaximum(nobs[ch][CYCLE-1]+2000);
    g_recovery2[ch] ->SetMarkerStyle(20);
    g_recovery2[ch] -> GetXaxis() -> SetRangeUser(0,ntrue[ch][CYCLE-1]+5000);
    g_recovery2[ch]->SetTitle(Form("recovery 2param (serial:%d);Nref (#mus^{-1});Nobs (#mus^{-1})",serial_num[ch]));
    g_recovery2[ch]->GetYaxis()->SetTitleOffset(1.05);
    g_recovery2[ch]->GetYaxis()->SetTitleSize(0.05);
    g_recovery2[ch]->GetXaxis()->SetTitleSize(0.05);
    g_recovery2[ch] -> Draw("ap");
    c2[ch]->cd(3);
    g_recovery3[ch] -> SetMaximum(nobs[ch][CYCLE-1]+2000);
    g_recovery3[ch] ->SetMarkerStyle(20);
    g_recovery3[ch] -> GetXaxis() -> SetRangeUser(0,ntrue[ch][CYCLE-1]+5000);
    g_recovery3[ch]->SetTitle(Form("recovery modified Vop (serial:%d);Nref (#mus^{-1});Nobs (#mus^{-1})",serial_num[ch]));
    g_recovery3[ch]->GetYaxis()->SetTitleOffset(1.05);
    g_recovery3[ch]->GetYaxis()->SetTitleSize(0.05);
    g_recovery3[ch]->GetXaxis()->SetTitleSize(0.05);
    g_recovery3[ch] -> Draw("ap");
    c2[ch]->cd(4);
    g_recovery4[ch] -> SetMaximum(nobs[ch][CYCLE-1]+2000);
    g_recovery4[ch] ->SetMarkerStyle(20);
    g_recovery4[ch] -> GetXaxis() -> SetRangeUser(0,ntrue[ch][CYCLE-1]+5000);
    g_recovery4[ch]->SetTitle(Form("recovery 2param + modified Vop (serial:%d);Nref (#mus^{-1});Nobs (#mus^{-1})",serial_num[ch]));
    g_recovery4[ch]->GetYaxis()->SetTitleOffset(1.05);
    g_recovery4[ch]->GetYaxis()->SetTitleSize(0.05);
    g_recovery4[ch]->GetXaxis()->SetTitleSize(0.05);
    g_recovery4[ch] -> Draw("ap");
    
    c2[ch]->cd(5);
    g_error1[ch] -> SetMarkerStyle(20);
    g_error1[ch] -> GetXaxis() -> SetRangeUser(0,ntrue[ch][CYCLE-1]+1000);
    g_error1[ch] -> SetMinimum(GERROR_MIN);
    g_error1[ch] -> SetMaximum(GERROR_MAX);
    g_error1[ch]->SetTitle(Form("recovery original error (serial:%d);Nref (#mus^{-1});Data/Fit",serial_num[ch]));
    g_error1[ch]->GetYaxis()->SetTitleOffset(1.05);
    g_error1[ch]->GetYaxis()->SetTitleSize(0.05);
    g_error1[ch]->GetXaxis()->SetTitleSize(0.05);
    g_error1[ch] -> Draw("ap");
    l[ch] = new TLine(0.,1.,ntrue[ch][CYCLE-1]+1000.,1.);
    l[ch] -> SetLineWidth(1);
    l[ch] -> SetLineStyle(2);
    l[ch] -> Draw("same");

    c2[ch]->cd(6);
    g_error2[ch] -> SetMarkerStyle(20);
    g_error2[ch] -> GetXaxis() -> SetRangeUser(0,ntrue[ch][CYCLE-1]+1000);
    g_error2[ch] -> SetMinimum(GERROR_MIN);
    g_error2[ch] -> SetMaximum(GERROR_MAX);
    g_error2[ch]->SetTitle(Form("recovery 2param error (serial:%d);Nref (#mus^{-1});Data/Fit",serial_num[ch]));
    g_error2[ch]->GetYaxis()->SetTitleOffset(1.05);
    g_error2[ch]->GetYaxis()->SetTitleSize(0.05);
    g_error2[ch]->GetXaxis()->SetTitleSize(0.05); 
    g_error2[ch] -> Draw("ap");
    l[ch] -> Draw("same");

    c2[ch]->cd(7);
    g_error3[ch] -> SetMarkerStyle(20);
    g_error3[ch] -> GetXaxis() -> SetRangeUser(0,ntrue[ch][CYCLE-1]+1000);
    g_error3[ch] -> SetMinimum(GERROR_MIN);
    g_error3[ch] -> SetMaximum(GERROR_MAX);
    g_error3[ch]->SetTitle(Form("recovery modified Vop error (serial:%d);Nref (#mus^{-1});Data/Fit",serial_num[ch]));
    g_error3[ch]->GetYaxis()->SetTitleOffset(1.05);
    g_error3[ch]->GetYaxis()->SetTitleSize(0.05);
    g_error3[ch]->GetXaxis()->SetTitleSize(0.05);
    g_error3[ch] -> Draw("ap");
    l[ch] -> Draw("same");

    c2[ch]->cd(8);
    g_error4[ch] -> SetMarkerStyle(20);
    g_error4[ch] -> GetXaxis() -> SetRangeUser(0,ntrue[ch][CYCLE-1]+1000);
    g_error4[ch] -> SetMinimum(GERROR_MIN);
    g_error4[ch] -> SetMaximum(GERROR_MAX);
    g_error4[ch]->SetTitle(Form("recovery 2param + modified Vop error (serial:%d);Nref (#mus^{-1});Data/Fit",serial_num[ch]));
    g_error4[ch]->GetYaxis()->SetTitleOffset(1.05);
    g_error4[ch]->GetYaxis()->SetTitleSize(0.05);
    g_error4[ch]->GetXaxis()->SetTitleSize(0.05);
    g_error4[ch] -> Draw("ap");
    l[ch] -> Draw("same");
  }

  TCanvas *c1_smallfit[T_NREADCH], *c1_smallrec[T_NREADCH];
  TCanvas *c2_orgerr, *c2_twoerr;

  
  if(DRAW_CANVAS){
    c2_orgerr = new TCanvas(Form("orgerr_%s",setup.c_str()),Form("orgerr_%s",setup.c_str()),0,0,(T_NREADCH-1)*400,300);
    c2_orgerr -> Divide(T_NREADCH-1,1,0.01,0.01);
    c2_twoerr = new TCanvas(Form("twoerr_%s",setup.c_str()),Form("twogerr_%s",setup.c_str()),0,300,(T_NREADCH-1)*400,300);
    c2_twoerr -> Divide(T_NREADCH-1,1,0.01,0.01);    

    for(int ch=1; ch<T_NREADCH; ch++){
      c2_orgerr -> cd(ch);
      g_error1[ch] -> Draw("ap");
      l[ch] -> Draw("same");
   
      c2_twoerr -> cd(ch);
      g_error2[ch] -> Draw("ap");
      l[ch] -> Draw("same");
    }

    
  }

  system("mkdir -p canvas/recovery");
  for(int ch=1; ch<T_NREADCH; ch++){
    c1[ch] -> Write();
    c2[ch] -> Write();
    if(PDFSAVE){
      c1[ch] -> Print(Form("canvas/recovery/small_ch%d.pdf",ch));
      c2[ch] -> Print(Form("canvas/recovery/rec_ch%d.pdf",ch));
      if(DRAW_CANVAS){	
	system("mkdir -p canvas/recovery/residual");
	c2_orgerr -> Print(Form("canvas/recovery/residual/recorgerr_%s.pdf",setup.c_str()));
	c2_twoerr -> Print(Form("canvas/recovery/residual/rectwoerr_%s.pdf",setup.c_str()));
      }
    }
  }
    newfile->Close();  
}

