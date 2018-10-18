#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <TString.h>
#include <iomanip>
//#include <vector>
#include "Root_h.h"

using namespace std;

const double THRESHOLD_FACTOR=5.0;
const double THRESHOLD=10.0; //mV unit
const int T_NSAMPLE=2048;
const int T_NREADCH=6;
const int T_NEVENT=5000;
const float AMP=100.;
const int CLOCK_WIDTH=50;
const int WIN_NUM=10;
const double WIN_WIDTH=1.2; //us units
const bool SKIP_EVENT=0; 
  

int main(int argc, char* argv[]){
  if(argc!=2){
    cout<<"usage : ./ana_dark (inputfile name)"<<endl;
    exit(1);
  }
  char* inputfile = argv[1];
  cout<<"input file : "<<Form("%s",inputfile)<<endl;
  //  TApplication app( "app", &argc, argv );
  
  TTree *t_wave;//, *t_time;
  TFile *readfile_ana = new TFile(inputfile);
  t_wave = (TTree*)readfile_ana->Get("rawwave");
  //t_time = (TTree*)readfile_ana->Get("timearray");


  int range_thre_bin,          range_thre_min,          range_thre_max;
  int range_thre_fit;
  int range_threlow_bin,      range_threlow_min,      range_threlow_max;
  int range_threlow_fit;
  int range_window_bin,       range_window_min,        range_window_max;
  int range_window_fit;
  int range_hybrid_bin,        range_hybrid_min,        range_hybrid_max;
  int range_hybrid_fit;
  int range_rise_bin,     range_rise_min,     range_rise_max;
  int range_baselinesigma_bin, range_baselinesigma_min, range_baselinesigma_max;
  range_thre_bin     = 2500; range_thre_min     = 0; range_thre_max = 2499;
  range_thre_fit     = 40;
  range_threlow_bin = 200;  range_threlow_min = 0; range_threlow_max = 199;
  range_threlow_fit     = 40;
  range_window_bin = 2600; range_window_min  = -100; range_window_max = 2499;
  range_window_fit     = 80;
  range_hybrid_bin      = 2600; range_hybrid_min      = -100; range_hybrid_max  = 2499;
  range_hybrid_fit     = 80;
  range_rise_bin     = 10; range_rise_min     = 0; range_rise_max = 100;
  range_baselinesigma_bin= 20; range_baselinesigma_min =0; range_baselinesigma_max = 10;


   // read parameter.dat file
  string tmpc;
  fstream file_param("parameter.dat",std::ios::in);
  if(!file_param){
    cout<<" parameter.dat not found !"<<endl;
  }else{
    cout<<" Found parameter.dat!"<<endl;
    while(file_param >> tmpc){
      if(tmpc=="range_thre"          ) file_param >> range_thre_bin     >> range_thre_min     >> range_thre_max;
      if(tmpc=="range_threfit"       ) file_param >> range_thre_fit;
      if(tmpc=="range_threlow"       ) file_param >> range_threlow_bin >> range_threlow_min >> range_threlow_max;
      if(tmpc=="range_threlowfit"       ) file_param >> range_threlow_fit;
      if(tmpc=="range_window"        ) file_param >> range_window_bin >> range_window_min >> range_window_max;
      if(tmpc=="range_windowfit"       ) file_param >> range_window_fit;
      if(tmpc=="range_hybrid"        ) file_param >> range_hybrid_bin      >> range_hybrid_min      >> range_hybrid_max;
      if(tmpc=="range_hybridfit"       ) file_param >> range_hybrid_fit;
      if(tmpc=="range_rise"          ) file_param >> range_rise_bin     >> range_rise_min     >> range_rise_max;
      if(tmpc=="range_baselinesigma" ) file_param >> range_baselinesigma_bin  >> range_baselinesigma_min  >> range_baselinesigma_max;
     }
  }
  

  //Histogram    
  TH1F* h_thre[T_NREADCH];
  TH1F* h_threlow[T_NREADCH];
  TH1F* h_window[T_NREADCH];
  TH1F* h_hybrid[T_NREADCH];
  TH1F* h_baselinesigma[T_NREADCH];
  TH2F* h_riseinteg[T_NREADCH];
  TGraph* g_avephoton[T_NREADCH];
  for ( int ch=0 ; ch<T_NREADCH ; ch++ ) {
    h_thre[ch] = new TH1F(Form("h_thre_ch%d",ch),Form("h_thre_ch%d",ch),range_thre_bin,range_thre_min,range_thre_max);
    h_threlow[ch] = new TH1F(Form("h_threlow_ch%d",ch),Form("h_threlow_ch%d",ch),range_threlow_bin,range_threlow_min,range_threlow_max);
    h_window[ch] = new TH1F(Form("h_window_ch%d",ch),Form("h_window_ch%d",ch),range_window_bin,range_window_min,range_window_max);
    h_hybrid[ch] = new TH1F(Form("h_hybrid_ch%d",ch),Form("h_hybrid_ch%d",ch),range_hybrid_bin,range_hybrid_min,range_hybrid_max);
    h_baselinesigma[ch] = new TH1F(Form("h_baselinesigma_ch%d",ch),Form("h_baselinesigma_ch%d",ch),range_baselinesigma_bin,range_baselinesigma_min,range_baselinesigma_max);
    h_riseinteg[ch] = new TH2F(Form("h_riseinteg_ch%d",ch),Form("h_riseinteg_ch%d",ch),range_threlow_bin,range_threlow_min,range_threlow_max/2,range_rise_bin,range_rise_min,range_rise_max);
  }

  
  // float thre_mean[T_NREADCH];
  TFile *new_file = new TFile(Form("analysed_%s",inputfile),"recreate");
  //TTree *tree_hybrid = new TTree("hybrid_integral","hybrid integral");
  TTree *tree_threlow = new TTree("threlow","threlow info");
  TTree *tree_hybrid = new TTree("hybrid","hybrid info");
  TTree *tree_threentry = new TTree("threentry","threshold entry");
  TTree *tree_winentry = new TTree("winentry","window entry");
  TTree *tree_interval = new TTree("interval","threshold interval");

  float t_fadc[T_NREADCH][T_NSAMPLE];

  
  TH1F* h_fadc[T_NREADCH];
  TF1* f_baseline[T_NREADCH];
  float baseline[T_NREADCH];
  float baseline_sigma[T_NREADCH];

  int skip_event=0;

  int b_event, b_ch;
  int threclock;
  float integral;
  int clockwidth;
  int riseclock, peakclock;
  float peak_hight;
  float photonwave[CLOCK_WIDTH]={};
  
  tree_threentry -> Branch("event",&b_event);
  tree_threentry -> Branch("ch",&b_ch);
  tree_threentry -> Branch("threclock",&threclock);
  tree_threentry -> Branch("integral",&integral);
  tree_threentry -> Branch("clockwidth",&clockwidth);
  tree_threentry -> Branch("riseclock",&riseclock);
  tree_threentry -> Branch("peakclock",&peakclock);
  tree_threentry -> Branch("peak_hight",&peak_hight);
  tree_threentry -> Branch("photonwave",photonwave,Form("photonwave[%d]/F",CLOCK_WIDTH));
    
  int win_smpwidth=(int)(WIN_WIDTH*100);
  int b_windownum;
  int win_abssmp;
  float window_integral;  
  int window_photonflag;
  float window_wave[win_smpwidth];
  tree_winentry -> Branch("event",&b_event);
  tree_winentry -> Branch("window_num",&b_windownum);
  tree_winentry -> Branch("ch",&b_ch);
  tree_winentry -> Branch("integral",&window_integral);
  tree_winentry -> Branch("photonflag",&window_photonflag);
  tree_winentry -> Branch("windowwave",window_wave,Form("window_wave[%d]/F",win_smpwidth));

  int intervalclock;
  int prethreclock;
  tree_interval -> Branch("event",&b_event);
  tree_interval -> Branch("ch",&b_ch);
  tree_interval -> Branch("intervalclock",&intervalclock);
  
  
  int photon_winnum[T_NREADCH]={};
  int threlow_darknum[T_NREADCH]={}, window_darknum[T_NREADCH]={};
  float avephoton[T_NREADCH][40]={}, avephoton_time[T_NREADCH][40], avephoton_num[T_NREADCH]={};
 
  /*  for(int ch=0; ch<T_NREADCH; ch++){
    //    tree_hybrid -> Branch(Form("ch%d", ch),&window_integral[ch]);
    //    tree_photonnum -> Branch(Form("ch%d", ch),&photon_avenum[ch]);
    //tree_hybrid -> Branch(Form("ch%d", ch), &eff_gain[ch]);
    }*/
  
  for(int event=0 ; event<T_NEVENT ; event++){ //Loop of ALL events (for analisys of dark current)
    if(event==0){cout<<"#### start! ####"<<endl;}
    if(event%1000==0 && event!=0){cout<<"#####"<<event<<"/"<<T_NEVENT<<"#####"<<endl;}
    t_wave->GetEntry(event);
    for(int ch=0 ; ch<T_NREADCH ; ch++){
      for(int i=0 ; i<T_NSAMPLE ; i++){
	t_fadc[ch][i] = t_wave->GetLeaf(Form("ch%d",ch))->GetValue(i);
      }
    }

    //Baseline determination method
    int hist_peak;
    for(int ch=0 ; ch<T_NREADCH ; ch++){
      h_fadc[ch] = new TH1F(Form("h_fadc_%d",ch),Form("h_fadc_%d",ch),2000,-2000,0); 
      for(int i=0 ; i<T_NSAMPLE; i++){
	h_fadc[ch]->Fill(t_fadc[ch][i]);
      }
      hist_peak = 0;
      for(int i=0 ; i<h_fadc[ch]->GetXaxis()->GetNbins() ; i++){
	if(h_fadc[ch]->GetBinContent(hist_peak) < h_fadc[ch]->GetBinContent(i)){
	  hist_peak = i;
	}
      }

      hist_peak = (int)(h_fadc[ch]->GetXaxis()->GetXmin() + hist_peak*(h_fadc[ch]->GetXaxis()->GetXmax() - h_fadc[ch]->GetXaxis()->GetXmin())/(h_fadc[ch]->GetXaxis()->GetNbins()));
      f_baseline[ch] = new TF1(Form("f_baseline_%d",ch),"[0]/[2]/sqrt(2.*TMath::Pi())*exp(-0.5*pow((x-[1])/[2],2))", -2000, 0);
      f_baseline[ch]->SetParameters(h_fadc[ch]->GetBinContent(hist_peak),hist_peak,4);
      h_fadc[ch]->Fit(Form("f_baseline_%d",ch),"Q","",hist_peak-4,hist_peak+4);
      
      baseline[ch]       = f_baseline[ch]->GetParameter(1);
      baseline_sigma[ch] = f_baseline[ch]->GetParameter(2);
      h_baselinesigma[ch]->Fill(baseline_sigma[ch]);
      delete h_fadc[ch];
    }


    int led_flag=0;
    int noise_flag=0;    
    for(int ch=1; ch<T_NREADCH; ch++){
      for(int i=0; i<T_NSAMPLE; i++){
	if(t_fadc[ch][i]<baseline[ch]-1000)  led_flag++;
      }
      if(baseline_sigma[ch]>2.5) noise_flag++;
    }
    if(led_flag>0){
      cout<<Form("Error LED event  ### skip EVENT %d ###",event)<<endl;
      skip_event++;
      continue;
    }else if(noise_flag>0 && SKIP_EVENT==1){
      //cout<<Form("Error noisy event  ### skip EVENT %d ###",event)<<endl;
      skip_event++;
      continue;
    }

    //threshold section
    for(int ch=0 ; ch<T_NREADCH ; ch++){ //Loop of all channel
      for(int i=0 ; i<T_NSAMPLE ; i++){ //Loop in AN event
	if(i==0) prethreclock=0;
	if(t_fadc[ch][i]<baseline[ch]-THRESHOLD){
	  b_event=event;
	  b_ch=ch;
	  integral=0;
	  clockwidth=0;
	  riseclock=0;
	  threclock=i;
	  peakclock=i;
	  peak_hight=baseline[ch]-t_fadc[ch][i];
	  for(int j=-2; j<0; j++){
	  integral += baseline[ch] - t_fadc[ch][i+j];
	  clockwidth++;
	  }
	  for(int j=0; j<CLOCK_WIDTH; j++){
	    photonwave[j]=baseline[ch] - t_fadc[ch][i+j-2];
	  }

	  while(t_fadc[ch][i]<baseline[ch]/*-baseline_sigma[ch]*/){
	    integral += baseline[ch] - t_fadc[ch][i];
	    if(t_fadc[ch][i-1]>baseline[ch]-THRESHOLD && t_fadc[ch][i]<baseline[ch]-THRESHOLD){
	      threlow_darknum[ch]++;
	      intervalclock=i-prethreclock;
	      tree_interval -> Fill();
	      prethreclock=i;
	    }
	    if(baseline[ch]-t_fadc[ch][i]>peak_hight){
	      peak_hight=baseline[ch]-t_fadc[ch][i];
	      peakclock=i;
	    }	  
	    clockwidth++;
	    i++;
	    if(i==T_NSAMPLE) break;
	  }
	  riseclock=peakclock-threclock;
	  tree_threentry -> Fill();

	  if(threclock>3 && threclock<2023 && integral>100 && integral<160){
	    for(int j=0; j<clockwidth; j++){
	      avephoton[ch][j]=(avephoton[ch][j]*avephoton_num[ch]+baseline[ch]-t_fadc[ch][threclock-5+j])/(avephoton_num[ch]+1);
	    }
	    avephoton_num[ch]++;
	  }

	  if(threclock>3 && threclock<2023){
	    h_thre[ch] -> Fill(integral);
	    h_threlow[ch] -> Fill(integral);
	    h_riseinteg[ch] -> Fill(integral, (float)(riseclock*10));
	    }
	}
      } //Loop in AN event 
    }//Loop in all channel


    //hybrid section
      for(int win=0; win<WIN_NUM; win++){//Loop of all windows in AN event
	for(int ch=0; ch<T_NREADCH; ch++){//Loop of all channel
	  b_event=event;
	  b_windownum=win;
	  b_ch=ch;
	  window_integral=0;
	  window_photonflag=0;
	  for(int win_smp=0; win_smp<win_smpwidth; win_smp++){
	    win_abssmp=20+T_NSAMPLE/WIN_NUM*win+win_smp;
	    window_integral+=baseline[ch]-t_fadc[ch][win_abssmp];
	    window_wave[win_smp] = baseline[ch]-t_fadc[ch][win_abssmp];
	    if(t_fadc[ch][win_abssmp]<baseline[ch]-THRESHOLD){
	      window_photonflag++;
	    }
 	    if(t_fadc[ch][win_abssmp-1]>baseline[ch]-THRESHOLD && t_fadc[ch][win_abssmp]<baseline[ch]-THRESHOLD){
	      window_darknum[ch]++;
	    }
	  }
	  if(window_photonflag==0){
	  h_window[ch] -> Fill(window_integral);
	  }
	  if(window_photonflag>0){
	    h_hybrid[ch] -> Fill(window_integral);
	    photon_winnum[ch]++;
	  }
	  tree_winentry -> Fill();
	}//Loop of all channel
      } //Loop of all windows in AN event 
  } //Loop of ALL events

  for(int ch=0; ch<T_NREADCH; ch++){
    for(int i=0; i<40; i++){
      avephoton_time[ch][i]=10*i;
    }
    g_avephoton[ch] = new TGraph(40, avephoton_time[ch], avephoton[ch]);
    g_avephoton[ch] -> SetName(Form("g_avephoton_ch%d",ch));
    g_avephoton[ch] -> SetTitle(Form("g_avephoton_ch%d",ch));
  }

  int channel[T_NREADCH]={};
  for(int ch=0; ch<T_NREADCH; ch++) channel[ch]=ch; 
  TF1 *f_threlow1[T_NREADCH], *f_threlow2[T_NREADCH];
  TF1 *f_window[T_NREADCH];
  TF1 *f_hybrid[T_NREADCH];
  int threlow_peak[T_NREADCH]={}, threlow_peakbin[T_NREADCH]={};
  int window_peak[T_NREADCH]={}, window_peakbin[T_NREADCH]={};
  int hybrid_peak[T_NREADCH]={}, hybrid_peakbin[T_NREADCH]={}; 
  float threlow_peakfit1[T_NREADCH]={}, threlow_peakfit2[T_NREADCH]={}, window_peakfit[T_NREADCH]={}, hybrid_peakfit[T_NREADCH];
  float threlow_opadc[T_NREADCH]={}, threlow_histadcmean[T_NREADCH]={}, threlow_opgain[T_NREADCH]={}, threlow_effgain[T_NREADCH]={}, threlow_carate[T_NREADCH]={};
  float hybrid_opadc[T_NREADCH]={}, hybrid_histadcmean[T_NREADCH]={}, hybrid_opgain[T_NREADCH]={}, hybrid_effgain[T_NREADCH]={}, hybrid_carate[T_NREADCH]={};
  int allwindow_num=(T_NEVENT-skip_event)*WIN_NUM;
  int zerophoton_winnum[T_NREADCH]={};
  float average_photonnum[T_NREADCH]={};
  float threlow_darkrate[T_NREADCH]={}, window_darkrate[T_NREADCH]={};

  tree_threlow -> Branch("ch",channel,Form("channel[%d]/I",T_NREADCH));
  tree_threlow -> Branch("onephoton_gain",threlow_opgain,Form("threlow_opgain[%d]/F",T_NREADCH));
  tree_threlow -> Branch("eff_gain",threlow_effgain,Form("threlow_effgain[%d]/F",T_NREADCH));
  tree_threlow -> Branch("ca_rate",threlow_carate,Form("threlow_carate[%d]/F",T_NREADCH));
  tree_threlow -> Branch("dark_rate",threlow_darkrate,Form("threlow_darkrate[%d]/F",T_NREADCH));

  tree_hybrid -> Branch("ch",channel,Form("channel[%d]/I",T_NREADCH));
  tree_hybrid -> Branch("onephoton_gain",hybrid_opgain,Form("hybrid_opgain[%d]/F",T_NREADCH));
  tree_hybrid -> Branch("eff_gain",hybrid_effgain,Form("hybrid_effgain[%d]/F",T_NREADCH));
  tree_hybrid -> Branch("ca_rate",hybrid_carate,Form("hybrid_carate[%d]/F",T_NREADCH));
  tree_hybrid -> Branch("dark_rate",window_darkrate,Form("window_darkrate[%d]/F",T_NREADCH));
  tree_hybrid -> Branch("abverage_photonnum",average_photonnum,Form("average_photonnum[%d]/F",T_NREADCH));
    



  for(int ch=1; ch<T_NREADCH; ch++){
    for(int i=0; i<h_threlow[ch]->GetXaxis()->GetNbins(); i++){
      if(h_threlow[ch]->GetBinContent(threlow_peak[ch]) < h_threlow[ch]->GetBinContent(i)){
	threlow_peak[ch] = i;
	threlow_peakbin[ch]=(int)h_threlow[ch]->GetBinContent(threlow_peak[ch]);
      }
    }
    threlow_peak[ch] = (int)(h_threlow[ch]->GetXaxis()->GetXmin() + threlow_peak[ch]*(h_threlow[ch]->GetXaxis()->GetXmax() - h_threlow[ch]->GetXaxis()->GetXmin())/h_threlow[ch]->GetXaxis()->GetNbins());
    f_threlow1[ch] = new TF1(Form("f_threlow1_%d",ch),"[0]/[2]/sqrt(2.*TMath::Pi())*exp(-0.5*pow((x-[1])/[2],2))", threlow_peak[ch]-40,threlow_peak[ch]+40);
    f_threlow2[ch] = new TF1(Form("f_threlow2_%d",ch),"[0]/[2]/sqrt(2.*TMath::Pi())*exp(-0.5*pow((x-[1])/[2],2))", 2*threlow_peak[ch]-40,2*threlow_peak[ch]+40);
    f_threlow1[ch]->SetParameters(threlow_peakbin[ch],threlow_peak[ch],25);
    f_threlow2[ch]->SetParameters(threlow_peakbin[ch]/10,2*threlow_peak[ch],25);
    h_threlow[ch] -> Fit(Form("f_threlow1_%d",ch),"+Q","",threlow_peak[ch]-range_threlow_fit,threlow_peak[ch]+range_threlow_fit);
    h_threlow[ch] -> Fit(Form("f_threlow2_%d",ch),"+Q","",2*threlow_peak[ch]-range_threlow_fit,2*threlow_peak[ch]+range_threlow_fit);
    threlow_peakfit1[ch] = f_threlow1[ch]->GetParameter(1);
    threlow_peakfit2[ch] = f_threlow2[ch]->GetParameter(1);
    threlow_opadc[ch] = threlow_peakfit2[ch] - threlow_peakfit1[ch];
    threlow_opgain[ch] = threlow_opadc[ch]/(50*AMP*1.602)*pow(10.,8.);
    threlow_histadcmean[ch] = h_threlow[ch]->GetMean();
    threlow_effgain[ch] = (threlow_histadcmean[ch]-threlow_peakfit1[ch]+threlow_opadc[ch])/(50*AMP*1.602)*pow(10.,8.);
    threlow_carate[ch] = (1-threlow_opgain[ch]/threlow_effgain[ch])*100;

    threlow_darkrate[ch] = (float)threlow_darknum[ch]/(T_NSAMPLE*(T_NEVENT-skip_event))*pow(10.,8.);

    cout<<"ch"<<ch<<" threlow_opgain="<<threlow_opgain[ch]<<" threlow_effgain="<<threlow_effgain[ch]
	<<" carate="<<threlow_carate[ch]<<"%"<<" threlow_darkrate="<<threlow_darkrate[ch]<<"Hz"<<endl;


    for(int i=0; i<h_window[ch]->GetXaxis()->GetNbins(); i++){
      if(h_window[ch]->GetBinContent(window_peak[ch]) < h_window[ch]->GetBinContent(i)){
	window_peak[ch] = i;
	window_peakbin[ch]=(int)h_window[ch]->GetBinContent(window_peak[ch]);
      }
    }
    window_peak[ch] = (int)(h_window[ch]->GetXaxis()->GetXmin() + window_peak[ch]*(h_window[ch]->GetXaxis()->GetXmax() - h_window[ch]->GetXaxis()->GetXmin())/h_window[ch]->GetXaxis()->GetNbins());
    f_window[ch] = new TF1(Form("f_window_%d",ch),"[0]/[2]/sqrt(2.*TMath::Pi())*exp(-0.5*pow((x-[1])/[2],2))", window_peak[ch]-40,window_peak[ch]+40);
    f_window[ch]->SetParameters(window_peakbin[ch],window_peak[ch],25);
    h_window[ch] -> Fit(Form("f_window_%d",ch),"Q","",window_peak[ch]-range_window_fit,window_peak[ch]+range_window_fit);
    window_peakfit[ch] = f_window[ch]->GetParameter(1);
    
    for(int i=0; i<h_hybrid[ch]->GetXaxis()->GetNbins(); i++){
      if(h_hybrid[ch]->GetBinContent(hybrid_peak[ch]) < h_hybrid[ch]->GetBinContent(i)){
	hybrid_peak[ch] = i;
	hybrid_peakbin[ch]=(int)h_hybrid[ch]->GetBinContent(hybrid_peak[ch]);
      }
    }
    hybrid_peak[ch] = (int)(h_hybrid[ch]->GetXaxis()->GetXmin() + hybrid_peak[ch]*(h_hybrid[ch]->GetXaxis()->GetXmax() - h_hybrid[ch]->GetXaxis()->GetXmin())/h_hybrid[ch]->GetXaxis()->GetNbins());
    f_hybrid[ch] = new TF1(Form("f_hybrid_%d",ch),"[0]/[2]/sqrt(2.*TMath::Pi())*exp(-0.5*pow((x-[1])/[2],2))", hybrid_peak[ch]-40,hybrid_peak[ch]+40);
    f_hybrid[ch]->SetParameters(hybrid_peakbin[ch],hybrid_peak[ch],25);
    h_hybrid[ch] -> Fit(Form("f_hybrid_%d",ch),"Q","",hybrid_peak[ch]-range_hybrid_fit,hybrid_peak[ch]+range_hybrid_fit);
    hybrid_peakfit[ch] = f_hybrid[ch]->GetParameter(1);
    hybrid_opadc[ch] = hybrid_peakfit[ch] - window_peakfit[ch];
    hybrid_opgain[ch] = hybrid_opadc[ch]/(50*AMP*1.602)*pow(10.,8.);

    zerophoton_winnum[ch] = allwindow_num-photon_winnum[ch];
    average_photonnum[ch] = -log(zerophoton_winnum[ch]/(float)allwindow_num);
    hybrid_histadcmean[ch] = h_hybrid[ch]->GetMean();
    hybrid_effgain[ch] = (hybrid_histadcmean[ch]-window_peakfit[ch])*(h_hybrid[ch]->GetEntries())/(h_hybrid[ch]->GetEntries()+h_window[ch]->GetEntries())/(50*AMP*average_photonnum[ch]*1.602)*pow(10.,8.);
    hybrid_carate[ch] = (1-hybrid_opgain[ch]/hybrid_effgain[ch])*100;

    
    window_darkrate[ch]=(float)window_darknum[ch]/(WIN_WIDTH*allwindow_num)*pow(10.,6.);

    //    h_hybrid[ch]->Write();
    //h_window[ch]->Write();
    cout<<"ch"<<ch<<" hybrid_opgain="<<hybrid_opgain[ch]<<" hybrid_effgain="<<hybrid_effgain[ch]<<" carate="<<hybrid_carate[ch]<<"%"
	<<" average_photonnum="<<average_photonnum[ch]/*<<" hybrid_darknum="<<window_darknum[ch]*/<<" hybrid_darkrate="<<window_darkrate[ch]<<"Hz"<<endl;
    //    cout<<"ch"<<ch<<" threlow one peak="<<threlow_peakfit1[ch]<<" threlow two peak="<<threlow_peakfit2[ch]<<" window peak="<<window_peakfit[ch]<<" hybrid peak="<<hybrid_peakfit[ch]<<endl;
  }

    tree_threlow -> Fill();
    tree_hybrid -> Fill();
  
    tree_threlow -> Write();
    tree_hybrid -> Write();
    tree_threentry -> Write();
    tree_winentry -> Write();
    tree_interval -> Write();
      
  //Canvas
  TCanvas *c[T_NREADCH];
  TCanvas *c_threlow[T_NREADCH], *c_window[T_NREADCH], *c_hybrid[T_NREADCH];
  for(int ch=1 ; ch<T_NREADCH ; ch++){
    c[ch] = new TCanvas(Form("c_%d",ch),Form("ch%d",ch),100*(ch-1),0,1800,800);
    c_threlow[ch] = new TCanvas(Form("c_threlow%d",ch),Form("ch%d",ch),100*(ch-1),100,600,450);
    c_window[ch] = new TCanvas(Form("c_window%d",ch),Form("ch%d",ch),100*(ch-1),200,600,450);
    c_hybrid[ch] = new TCanvas(Form("c_hybrid%d",ch),Form("ch%d",ch),100*(ch-1),300,600,450);
  }

  for(int ch=1 ; ch<T_NREADCH ; ch++){
    c[ch]->Divide(3,2,0.01,0.01);
  }
   for(int ch=1 ; ch<T_NREADCH ; ch++){
    c[ch]->cd(1)->SetLogy(0);
    h_threlow[ch]->GetXaxis()->SetTitle("Integral");
    h_threlow[ch]->GetYaxis()->SetTitle("count");
    h_threlow[ch]->Draw();
  }
   for(int ch=1 ; ch<T_NREADCH ; ch++){
    c[ch]->cd(2)->SetLogy(0);
    h_window[ch]->GetXaxis()->SetTitle("Integral");
    h_window[ch]->GetYaxis()->SetTitle("count");
    h_window[ch]->Draw();
  }
  for(int ch=1 ; ch<T_NREADCH ; ch++){
    c[ch]->cd(3)->SetLogy(0);
    h_hybrid[ch]->GetXaxis()->SetTitle("Integral");
    h_hybrid[ch]->GetYaxis()->SetTitle("count");
    h_hybrid[ch]->Draw();
  }
  for(int ch=1 ; ch<T_NREADCH ; ch++){
    c[ch]->cd(4)->SetLogy(1);
    h_baselinesigma[ch]->GetXaxis()->SetTitle("Baseline sigma");
    h_baselinesigma[ch]->GetYaxis()->SetTitle("count");
    h_baselinesigma[ch]->Draw();
  }
  for(int ch=1; ch<T_NREADCH; ch++){
    c[ch]->cd(5);
    h_riseinteg[ch]->GetXaxis()->SetTitle("Integral");
    h_riseinteg[ch]->GetYaxis()->SetTitle("risetime");
    h_riseinteg[ch]->SetMaximum(100);
    h_riseinteg[ch]->Draw("colz");
  }
  for(int ch=1; ch<T_NREADCH; ch++){
    c[ch]->cd(6);
    g_avephoton[ch]->GetXaxis()->SetTitle("time (ns)");
    g_avephoton[ch]->GetYaxis()->SetTitle("amplitude (mV)");
    g_avephoton[ch]->SetMarkerStyle(20);
    g_avephoton[ch]->Draw("ap");
  }

  for(int ch=1; ch<T_NREADCH; ch++){
    c_threlow[ch]->cd()->SetLogy(0);
    h_threlow[ch]->GetXaxis()->SetTitle("Integral");
    h_threlow[ch]->GetYaxis()->SetTitle("count");
    h_threlow[ch]->Draw();
  }
  for(int ch=1; ch<T_NREADCH; ch++){
    c_window[ch]->cd()->SetLogy(0);
    h_window[ch]->GetXaxis()->SetTitle("Integral");
    h_window[ch]->GetYaxis()->SetTitle("count");
    h_window[ch]->Draw();
  }
  for(int ch=1; ch<T_NREADCH; ch++){
    c_hybrid[ch]->cd()->SetLogy(0);
    h_hybrid[ch]->GetXaxis()->SetTitle("Integral");
    h_hybrid[ch]->GetYaxis()->SetTitle("count");
    h_hybrid[ch]->Draw();
  }

  
  for(int ch=1; ch<T_NREADCH; ch++){
    c[ch]->Write();
  }
  new_file -> Close();
  
  system("mkdir -p canvas/dark/each");
  
    for(int ch=1; ch<T_NREADCH; ch++){
      c[ch]->Print(Form("canvas/dark/c_%d.pdf",ch));
      c_threlow[ch]->Print(Form("canvas/dark/each/c_threlow%d.pdf",ch));
      c_window[ch]->Print(Form("canvas/dark/each/c_window%d.pdf",ch));
      c_hybrid[ch]->Print(Form("canvas/dark/each/c_hybrid%d.pdf",ch));
     }

    //    app.Run();
    return 0;
}
