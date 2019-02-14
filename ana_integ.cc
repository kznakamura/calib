#include <iostream>
#include <string>
//#include <fstream>
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <TString.h>
//#include <iomanip>
//#include "Root_h.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TLeaf.h>

using namespace std;
/*
const int T_NSAMPLE=2048;
const int T_NREADCH=2;
const int T_NEVENT=1000;
*/
const int MAXCYCLE=50;
const int MAXSAMPLING=2048;

class Analizer{
public:
  Analizer(string filename, int input_cycle);
  //~Analizer();

private:
  int m_MAXCYCLE = MAXCYCLE;
  int m_MAXSAMPLING = MAXSAMPLING;
  string m_ifilename;
  int m_input_cycle;
  int m_readch;

  int m_integral_min = 1000, m_integral_max = 1200;
  int m_dark_min = 400, m_dark_max = 600;
  
  int m_ch;
  double m_integral[MAXCYCLE] = {};
  double m_integral_error[MAXCYCLE] = {};
  double m_dark[MAXCYCLE] = {}, m_dark_error[MAXCYCLE] = {};
  
  int m_cycle_event;
  int m_ch_event;
  double m_integral_event;
  double m_dark_event;
  double m_baseline, m_baseline_sigma;
};

Analizer::Analizer(string filename, int input_cycle){
  
  m_ifilename = filename;
  m_input_cycle = input_cycle;

  TFile *ifile[m_input_cycle];
  TTree *rawwave_tree[m_input_cycle];
  
  for(int cy=0; cy<m_input_cycle; cy++){
    ifile[cy] = new TFile(Form("%s_%02d.root",m_ifilename.c_str(),cy+1));
    cout<<Form("%s_%02d.root",m_ifilename.c_str(),cy+1)<<endl;
    rawwave_tree[cy] = (TTree*)ifile[cy]->Get("rawwave");
    m_readch = rawwave_tree[cy] -> GetNbranches();
  }
  
  TFile *ofile = new TFile(Form("sum_%s.root", m_ifilename.c_str()), "recreate");
  TTree *integral_header = new TTree("integral_header","integral_header");
  TTree *integral_tree = new TTree("integral_tree", "integral_tree");
  TTree *integral_event_tree = new TTree("integral_event_tree", "integral_event_tree");
  
  integral_header -> Branch("integral_min", &m_integral_min);
  integral_header -> Branch("integral_max", &m_integral_max);
  integral_header -> Branch("dark_min", &m_dark_min);
  integral_header -> Branch("dark_max", &m_dark_max);
  integral_header -> Fill();
  integral_header -> Write();
    
  integral_tree -> Branch("cycle", &m_input_cycle, "cycle/I");
  integral_tree -> Branch("ch", &m_ch);
  integral_tree -> Branch("integral", m_integral, "integral[cycle]/D");
  integral_tree -> Branch("integral_error", m_integral_error, "integral_error[cycle]/D");
  integral_tree -> Branch("dark", m_dark, "dark[cycle]/D");
  integral_tree -> Branch("dark_error", m_dark_error, "dark_error[cycle]/D");

  integral_event_tree -> Branch("cycle_event", &m_cycle_event);
  integral_event_tree -> Branch("ch", &m_ch);
  integral_event_tree -> Branch("integral_event", &m_integral_event);
  integral_event_tree -> Branch("dark_event", &m_dark_event);
  integral_event_tree -> Branch("baseline", &m_baseline);
  integral_event_tree -> Branch("baseline_sigma", &m_baseline_sigma);
    
  TApplication app("app",0,0,0,0);
  TCanvas *c = new TCanvas(1);
  // c -> DrawFrame(-400,0,-200,500);
  
  for(int ch=0; ch<m_readch; ch++){
    m_ch = ch;
    for(int cy=0; cy<m_input_cycle; cy++){
      m_cycle_event = cy;
      //cout << cy << endl;
      int max_event = rawwave_tree[cy] -> GetEntries();
      float wave[MAXSAMPLING] = {};
      rawwave_tree[cy] -> SetBranchAddress(Form("ch%d",ch), wave);
      int wave_length = rawwave_tree[cy] -> GetLeaf(Form("ch%d",ch)) -> GetLen();
      
      TH1D *h_integral = new TH1D("h_integral", "h_integral", 20020, -200, 200000);
      TH1D *h_dark = new TH1D("h_dark", "h_dark", 20020, -200, 200000);
      
      for(int ev=0; ev<max_event; ev++){
	rawwave_tree[cy] -> GetEntry(ev);
	TH1D *h_baseline = new TH1D("h_baseline", "h_baseline", pow(2,14), -2250, 0);	
	for(int smp=0; smp<wave_length; smp++){
	  h_baseline -> Fill(wave[smp]);
	}
	TF1 *f_baseline = new TF1("f_baseline", "gaus", -2250, 0);
        double baseline_peak_bin = -2250.0 + h_baseline->GetBinWidth(0) * h_baseline -> GetMaximumBin();
	h_baseline -> Fit("f_baseline", "Q", "", baseline_peak_bin-5, baseline_peak_bin+5);
	m_baseline = f_baseline -> GetParameter(1);
	m_baseline_sigma = f_baseline -> GetParameter(2);
      
	m_integral_event = 0;
	for(int smp=m_integral_min; smp<m_integral_max; smp++){
	  m_integral_event += m_baseline - (double)wave[smp];
	}
	m_dark_event = 0;
	for(int smp=m_dark_min; smp<m_dark_max; smp++){
	  m_dark_event += m_baseline - (double)wave[smp];
	}
	/*	if(ev>500&&cy==5){
	  cout << "event=" << ev <<  endl;
	  rawwave_tree[cy] -> Draw("ch1:Iteration$",Form("Entry$==%d",ev),"l");
	  c-> Update();
	  gSystem -> ProcessEvents();
	  }*/
      
	integral_event_tree -> Fill();
	h_integral -> Fill(m_integral_event);
	h_dark -> Fill(m_dark_event);
	delete h_baseline;
	delete f_baseline;
	}
      
      double integral_peak_bin = -200 + h_integral->GetBinWidth(0) * h_integral->GetMaximumBin();
      double dark_peak_bin = -200 + h_dark->GetBinWidth(0) * h_dark->GetMaximumBin();
      double hist_integral_sigma = h_integral->GetStdDev();
      double hist_dark_sigma = h_dark->GetStdDev();
      TF1 *f_integral = new TF1("f_integral", "gaus", -200, 40000);
      TF1 *f_dark = new TF1("f_dark", "gaus", -200, 40000);
      h_integral -> Fit("f_integral", "Q", "", integral_peak_bin-hist_integral_sigma*3, integral_peak_bin+hist_integral_sigma*3);
      h_dark -> Fit("f_dark", "Q", "", dark_peak_bin-hist_dark_sigma*3, dark_peak_bin+hist_dark_sigma*3);
      double integral_mean = f_integral->GetParameter(1);
      double integral_sigma = f_integral->GetParameter(2);
      double dark_mean = f_dark->GetParameter(1);
      double dark_sigma = f_dark->GetParameter(2);
      if(cy<13){
	c->cd()->DrawFrame(-200,0,integral_peak_bin*2,200);
	h_integral -> Draw("same");
	c -> Update();
	gSystem -> ProcessEvents();
      }
      m_integral[cy] = integral_mean;
      m_integral_error[cy] = integral_sigma/sqrt(max_event);
      m_dark[cy] = dark_mean;
      m_dark_error[cy] = dark_sigma/sqrt(max_event);
      delete h_integral;
      delete h_dark;
      delete f_integral;
      delete f_dark;
    }
    integral_tree -> Fill();
  }
  integral_tree -> Write();
  integral_event_tree -> Write();
  
  for(int cy=0; cy<m_input_cycle; cy++){
    ifile[cy] -> Close();
  }
  ofile -> Close();

}


int main(int argc, char* argv[]){
  if(argc!=3){
    cout<<"usage : ./ana_dark (analisysfile name) (integral cycle)"<<endl;
    return -1;
  }
 
  string ifilename = argv[1];
  cout<<"input file : "<< ifilename <<endl;
  int cycle = atoi(argv[2]);

  Analizer *anal = new Analizer(ifilename, cycle);

  delete anal;
}

/* 
  
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
  
  
    
  TFile *rootfile = new TFile(Form("sum_%s.root",ifile),"recreate");
  TTree *tree_integ = new TTree("integral", "integral");
  TTree *tree_integerr = new TTree("integralerr", "integral err");
  TTree *tree_integentry = new TTree("integentry","integral entry");
  
  //Histogram
  TH1F* h_integral[T_NREADCH];
  
  float t_fadc[T_NREADCH][T_NSAMPLE];
  int b_cy, b_event, b_ch;
  float baseline;
  float integral;
  float integral_wave[range_integral_max-range_integral_min];
  float integral_mean[T_NREADCH]={}, integral_err[T_NREADCH]={};
  for(int ch=0; ch<T_NREADCH; ch++){
    tree_integ->Branch(Form("ch%d",ch),&integral_mean[ch]);

    tree_integerr->Branch(Form("ch%d",ch),&integral_err[ch]);
  }
  tree_integentry -> Branch("cycle",&b_cy);
  tree_integentry -> Branch("event",&b_event);
  tree_integentry -> Branch("ch",&b_ch);
  tree_integentry -> Branch("integ",&integral);
  tree_integentry -> Branch("integwave",integral_wave,Form("integral_wave[%d]/F",range_integral_max-range_integral_min));
  
  
  for(int cy=0; cy<CYCLE; cy++){
    b_cy=cy;
    cout<<Form("analysys start! %s_%02d.root",ifile,cy+1)<<endl;
    for( int ch=0 ; ch<T_NREADCH ; ch++ ) {
      h_integral[ch] = new TH1F(Form("h_integral_ch%d",ch),Form("h_integral_ch%d",ch),range_integralhist_bin,range_integralhist_min,range_integralhist_max);
    }
    
    for(int event=10 ; event<T_NEVENT ; event++){ //Loop of ALL events (for analisys of dark current)
      rawwave[cy]->GetEntry(event);
      b_event=event;
      for(int ch=0 ; ch<T_NREADCH ; ch++){
	for(int i=0 ; i<T_NSAMPLE ; i++){
	  t_fadc[ch][i] = rawwave[cy]->GetLeaf(Form("ch%d",ch))->GetValue(i);
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
*/
