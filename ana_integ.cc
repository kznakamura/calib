#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <TString.h>
#include <iomanip>
#include "Root_h.h"

using namespace std;

const int T_NSAMPLE=2048;
const int T_NREADCH=2;
const int T_NEVENT=1000;


int main(int argc, char* argv[]){
  if(argc!=3){
    cout<<"usage : ./ana_dark (analisysfile name) (integral cycle)"<<endl;
    exit(1);
  }
  char* inputfile = argv[1];
  cout<<"input file : "<<Form("%s",inputfile)<<endl;
  int CYCLE = atoi(argv[2]);
  
  //TApplication app( "app", &argc, argv );
  
  TTree *t_wave[CYCLE];//, *t_time[CYCLE];
  TFile *readfile_ana[CYCLE];
  
  //  TFile *readfile_ana = new TFile(inputfile);
  for(int cy=0; cy<CYCLE; cy++){
    //    if(cy<10){
      readfile_ana[cy] = new TFile(Form("%s_%02d.root",inputfile,cy+1));
      cout<<Form("%s_%02d.root",inputfile,cy+1)<<endl;
      t_wave[cy] = (TTree*)readfile_ana[cy]->Get("rawwave");
      //t_time[cy] = (TTree*)readfile_ana[cy]->Get("timearray");
  }
  int range_baseline_min, range_baseline_max;
  int range_integral_min, range_integral_max;
  int range_integralhist_bin, range_integralhist_min, range_integralhist_max;
  range_baseline_min     = 0; range_baseline_max     = 900;
  range_integral_min     = 900; range_integral_max     = 1800;
  range_integralhist_bin = 100; range_integralhist_min = 0; range_integralhist_max = 10000;
  
  // read parameter.dat file
  string tmpc;
  fstream file_param("parameter.dat",std::ios::in);
  if(!file_param){
    cout<<" parameter.dat not found !"<<endl;
  }else{
    cout<<" Found parameter.dat!"<<endl;
    while(file_param >> tmpc){
      if(tmpc=="range_baseline") file_param >> range_baseline_min     >> range_baseline_max;
      if(tmpc=="range_integral") file_param >> range_integral_min     >> range_integral_max;
      if(tmpc=="range_integralhist"){ file_param >> range_integralhist_bin     >> range_integralhist_min     >> range_integralhist_max;
	cout<<"ana_integ changed"<<endl;
      }
    }
  }

  
  
  TFile *rootfile = new TFile(Form("sum_%s.root",inputfile),"recreate");
  TTree *tree_integ = new TTree("integral", "integral");
  // TTree *tree_integdev = new TTree("integraldev", "integral stddev");
  TTree *tree_integerr = new TTree("integralerr", "integral err");
  TTree *tree_integentry = new TTree("integentry","integral entry");
  
  //Histogram
  TH1F* h_integral[T_NREADCH];
  
  float t_fadc[T_NREADCH][T_NSAMPLE];
  int b_cy, b_event, b_ch;
  float baseline;
  float integral;
  float integral_wave[range_integral_max-range_integral_min];
  float integral_mean[T_NREADCH]={}, integral_err[T_NREADCH]={}; //integral_dev[T_NREADCH]={};
  for(int ch=0; ch<T_NREADCH; ch++){
    tree_integ->Branch(Form("ch%d",ch),&integral_mean[ch]);//,"integral_mean[T_NREADCH]/F");
    //   tree_integdev->Branch(Form("ch%d",ch),&integral_dev[ch]);
    tree_integerr->Branch(Form("ch%d",ch),&integral_err[ch]);
  }
  tree_integentry -> Branch("cycle",&b_cy);
  tree_integentry -> Branch("event",&b_event);
  tree_integentry -> Branch("ch",&b_ch);
  tree_integentry -> Branch("integ",&integral);
  tree_integentry -> Branch("integwave",integral_wave,Form("integral_wave[%d]/F",range_integral_max-range_integral_min));
  
  
  for(int cy=0; cy<CYCLE; cy++){
    b_cy=cy;
    cout<<Form("analysys start! %s_%02d.root",inputfile,cy+1)<<endl;
    for( int ch=0 ; ch<T_NREADCH ; ch++ ) {
      h_integral[ch] = new TH1F(Form("h_integral_ch%d",ch),Form("h_integral_ch%d",ch),range_integralhist_bin,range_integralhist_min,range_integralhist_max);
    }
    
    for(int event=10 ; event<T_NEVENT ; event++){ //Loop of ALL events (for analisys of dark current)
      t_wave[cy]->GetEntry(event);
      b_event=event;
      for(int ch=0 ; ch<T_NREADCH ; ch++){
	for(int i=0 ; i<T_NSAMPLE ; i++){
	  t_fadc[ch][i] = t_wave[cy]->GetLeaf(Form("ch%d",ch))->GetValue(i);
	}
      }
      
      //Baseline determination method
      for(int ch=0 ; ch<T_NREADCH ; ch++){ //loop of all ch
	b_ch=ch;
	baseline=t_fadc[ch][range_baseline_min];
	integral=0;
	//TH1D *h_baseline = new TH1D("h_baseline","h_baseline",2000,-2000,0);
	for(int i=range_baseline_min ; i<range_baseline_max; i++){
	  baseline=(baseline*(i-range_baseline_min)+t_fadc[ch][i])/(i-range_baseline_min+1);
	  //h_baseline -> Fill(t_fadc[ch][i]);
	}
	/*	int max_bin = h_baseline -> GetMaximumBin();
	int hist_peak = -2000 + 2000*max_bin/2000;
	TF1 *f_baseline = new TF1("f_baseline","gaus",-2000,0);
	h_baseline -> Fit("f_baseline","Q","",hist_peak-4,hist_peak+4);
	baseline = f_baseline -> GetParameter(1);
	*/
	for(int i=range_integral_min; i<range_integral_max; i++){
	  integral+=baseline-t_fadc[ch][i];
	  integral_wave[i-range_integral_min]=baseline-t_fadc[ch][i];
	}
	h_integral[ch] -> Fill(integral);
	tree_integentry -> Fill();
	//delete h_baseline;
	//delete f_baseline;
      } //loop of all ch
    } //loop of all event
   
    for(int ch=0 ; ch<T_NREADCH ; ch++){
      integral_mean[ch]=(float)h_integral[ch]->GetMean();
      //      integral_dev[ch]=(float)h_integral[ch]->GetStdDev();
      integral_err[ch]=(float)(h_integral[ch]->GetStdDev()/(float)sqrt(T_NEVENT-10));
      //      cout<<"ch="<<ch<<" integral_mean="<<integral_mean[ch]<<" integral_dev="<<integral_dev[ch]<<endl;
      cout<<"ch="<<ch<<" integral_mean="<<integral_mean[ch]<<" integral_err="<<integral_err[ch]<<endl;
      
      //h_integral[ch]->Write();
      delete h_integral[ch];
    }
    tree_integ->Fill();
    //tree_integdev->Fill();
    tree_integerr->Fill();
  }

  tree_integ->Write();
  //tree_integdev->Write();
  tree_integerr->Write();
  tree_integentry -> Write();
  rootfile->Close();
    
  //  app.Run();
  return 0;
}
