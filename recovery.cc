#include <iostream>
#include <string>
#include <regex>
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <utility>
#include <TPaveText.h>

using namespace std;

const int MAXCH = 6;
const int MAXCYCLE = 25;

class Analizer{
public:
  Analizer(string target, string setup_num, string anadark_name, double fit_max, string anadate);
  ~Analizer();
  void saveTree();
  void drawCanvas();
  void saveCanvas();
private:
  int mem_readch = 2;
  double mem_fit_max;
  int mem_linearfit_num = 5;
  double mem_mppc_amp = 10.;
  double mem_sys_impedance = 50.;
  double mem_e_charge = 1.602e-19;
  double mem_mppc_pixel = 3600;
  double mem_led_width = 1000; //ns unit
  int mem_ch[MAXCH] = {};
  int mem_serial[MAXCH] = {};
  int mem_cycle = 0;
  double mem_effgain[MAXCH] = {};
  double mem_onepar_tau[MAXCH] = {};
  //  double mem_onepar_slope[MAXCH] = {};
  double mem_twopar_tau1[MAXCH] = {}, mem_twopar_tau2[MAXCH] = {};
  double mem_twopar_alpha[MAXCH] = {}, mem_twopar_beta[MAXCH] = {};
  double mem_integral_mean[MAXCH][MAXCYCLE] = {};
  double mem_integral_err[MAXCH][MAXCYCLE] = {};
  double mem_linear_residual[MAXCH][MAXCYCLE] = {};
  double mem_linear_residual_err[MAXCH][MAXCYCLE] = {};
  double mem_nobs_mean[MAXCH][MAXCYCLE] = {};
  double mem_nref_mean[MAXCH][MAXCYCLE] = {};
  double mem_nobs_err[MAXCH][MAXCYCLE] = {};
  double mem_nref_err[MAXCH][MAXCYCLE] = {};
  double mem_onepar_residual[MAXCH][MAXCYCLE] = {};
  double mem_twopar_residual[MAXCH][MAXCYCLE] = {};
  double mem_onepar_residual_err[MAXCH][MAXCYCLE] = {};
  double mem_twopar_residual_err[MAXCH][MAXCYCLE] = {};

  TGraphErrors *g_adc[MAXCH];
  TF1 *f_linear[MAXCH];
  TGraphErrors *g_linear_residual[MAXCH];
  TGraphErrors *g_photon[MAXCH];
  TF1 *f_onepar[MAXCH], *f_twopar[MAXCH];
  TGraphErrors *g_onepar_resi[MAXCH], *g_twopar_resi[MAXCH];
  TCanvas *canvas[MAXCH];
    
  TString mem_target;
  TString mem_setup_num;
  string mem_anadark_name;
  string mem_onepar_func = "x/(1+[0]*x)";
  string mem_twopar_func = "[2]*x/(1+[0]*x)+[3]*x/(1+[1]*x)";
  string mem_anadate;
  
  TString mem_input_anadark_name;
  TString mem_input_integ_name;
  TString mem_output_file_name;
  string mem_setupdir_path = "/home/kznakamura/DAQ/setup/";
  string mem_analyzeddir_path = "/home/kznakamura/DAQ/analyzed/";
  TString mem_analyzed_file_path;
  void setFileName(); 
  void inputSetupFile();
  void setValue();
  void setAnaFilePath(){mem_analyzed_file_path = mem_analyzeddir_path + mem_target + "/ana" + mem_anadate;} 
  void anaAdc2Photon();
  void anaRecovery();
  void makeCanvas();
};


Analizer::Analizer(string target, string setup_num, string anadark_name, double fit_max, string anadate){
  mem_target = target;
  mem_setup_num = setup_num;
  mem_anadark_name = anadark_name;
  mem_fit_max = fit_max;
  mem_anadate = anadate;
  cout << "target: " << target << ", setup" << setup_num << ", anadark: " << anadark_name << ", fit_max: " << fit_max <<", anadate:" << anadate << endl;
  setFileName();
  inputSetupFile();
  setValue();
  setAnaFilePath();
  anaAdc2Photon();
  anaRecovery();  
  makeCanvas();
}

Analizer::~Analizer(){
  cout << "Analizer is deleted" << endl;
}

void Analizer::setFileName(){
  mem_input_anadark_name =  mem_anadark_name;
  string dark_header1, dark_header2, dark_date, dark_condition, dark_num;
  regex rx(R"(_)");
  sregex_token_iterator it(mem_anadark_name.begin(), mem_anadark_name.end(), rx, -1);
  dark_header1 = *it++;
  dark_header2 = *it++;
  dark_date = *it++;
  dark_condition = *it++;
  dark_num = *it++;
  mem_input_integ_name = "sum_integ_" + dark_date + "_" + dark_condition + ".root";
  mem_output_file_name = mem_target + "_" + mem_setup_num + ".root";
  //cout << mem_input_anadark_name << endl;
  //cout << mem_input_integ_name << endl;
  //  cout << mem_output_file_name << endl;
  
}

void Analizer::inputSetupFile(){

  string setup_file =  (string)mem_target + "_setup.dat";
  string setup_file_path = mem_setupdir_path + setup_file;
  
  ifstream finsetup(setup_file_path);
  if(!finsetup){
    cout << "### can not find " << setup_file << " ###" << endl;
    exit(1);
  }else{
    cout << setup_file << " is found!!" << endl;
    string setup_scan;
    while(finsetup >> setup_scan) if(setup_scan == "setup" + mem_setup_num) break;
    for(int ch=1; ch<mem_readch; ch++){
      finsetup >> mem_ch[ch] >> mem_serial[ch];
      //cout << mem_ch[ch] << " " << mem_serial[ch] << endl;
    }
  }
  finsetup.close();  
}

void Analizer::setValue(){

  TFile *dark_file = new TFile(mem_input_anadark_name);
  TFile *integ_file = new TFile(mem_input_integ_name);
  TTree *dark_tree = (TTree*)dark_file -> Get("hybrid");
  TTree *integ_tree = (TTree*)integ_file -> Get("integral_tree");
  //  TTree *integ_err_tree = (TTree*)integ_file -> Get("integralerr");

  dark_tree -> GetEntry(0);
  for(int ch=1; ch<mem_readch; ch++){
    mem_effgain[ch] = (double)(dark_tree -> GetLeaf("hybrid_effgain") -> GetValue(ch));
    cout << mem_ch[ch] << " " << mem_serial[ch] << " " << mem_effgain[ch] << endl;
  }
  
  double integral_mean[MAXCYCLE] = {};
  double integral_error[MAXCYCLE] = {};
  double dark_mean[MAXCYCLE] = {};
  integ_tree -> SetBranchAddress("integral_mean",integral_mean);
  integ_tree -> SetBranchAddress("integral_error",integral_error);
  integ_tree -> SetBranchAddress("dark_mean", dark_mean);

  for(int ch=0; ch<mem_readch; ch++){
    integ_tree -> GetEntry(ch);
    mem_cycle = integ_tree -> GetLeaf("max_cycle") -> GetValue(0);
    for(int cy=0; cy<mem_cycle; cy++){
      mem_integral_mean[ch][cy] = integral_mean[cy] - dark_mean[cy];
      mem_integral_err[ch][cy] = integral_error[cy];
    }
  }
  dark_file -> Close();
  integ_file -> Close();
}

void Analizer::anaAdc2Photon(){
  for(int ch=1; ch<mem_readch; ch++){
    g_adc[ch] = new TGraphErrors(mem_cycle, mem_integral_mean[0], mem_integral_mean[ch], mem_integral_err[0], mem_integral_err[ch]);
    f_linear[ch] = new TF1(Form("f_linear_ch%d",ch),"[0]*x+[1]", -1000, 1.0e5);
    g_adc[ch] -> Fit(Form("f_linear_ch%d",ch),"Q0","",5,mem_integral_mean[0][mem_linearfit_num-1]);
    double slope  = f_linear[ch] -> GetParameter(0);
    double seppen = f_linear[ch] -> GetParameter(1);
    for(int cy=0; cy<mem_cycle; cy++){
      mem_linear_residual[ch][cy] = mem_integral_mean[ch][cy]/f_linear[ch]->Eval(mem_integral_mean[0][cy]);
      mem_linear_residual_err[ch][cy] = mem_integral_err[ch][cy]/f_linear[ch]->Eval(mem_integral_mean[0][cy]);
    }
    g_linear_residual[ch] = new TGraphErrors(mem_cycle, mem_integral_mean[0], mem_linear_residual[ch], mem_integral_err[ch], mem_linear_residual_err[ch]);
    for(int cy=0; cy<mem_cycle; cy++){
      mem_nobs_mean[ch][cy] 
	= mem_integral_mean[ch][cy]/(mem_e_charge*mem_effgain[ch]*mem_sys_impedance*mem_mppc_amp)*pow(10,-11);
      mem_nref_mean[ch][cy] 
	= (mem_integral_mean[0][cy]*slope + seppen)/(mem_e_charge*mem_effgain[ch]*mem_sys_impedance*mem_mppc_amp)*pow(10,-11);
      mem_nobs_err[ch][cy] 
	= mem_integral_err[ch][cy]/(mem_e_charge*mem_effgain[ch]*mem_sys_impedance*mem_mppc_amp)*pow(10,-11);
      mem_nref_err[ch][cy] 
	= (mem_integral_err[0][cy]*slope + seppen)/(mem_e_charge*mem_effgain[ch]*mem_sys_impedance*mem_mppc_amp)*pow(10,-11);
    }
    g_photon[ch] = new TGraphErrors(mem_cycle, mem_nref_mean[ch], mem_nobs_mean[ch], mem_nref_err[ch], mem_nobs_err[ch]);
  }
}

void Analizer::anaRecovery(){
  if(mem_fit_max<0){
    cout << "mem_fit_max is out of range" << endl;
    exit(1);
  }
  for(int ch=1; ch<mem_readch; ch++){
    if(mem_fit_max > mem_nref_mean[ch][mem_cycle-1]){
      cout << "mem_fit_max is changed to mem_nref_mean[ch][mem_cycle-1]" << endl;
      mem_fit_max = mem_nref_mean[ch][mem_cycle-1];
    }
    f_onepar[ch] = new TF1(Form("f_onepar_ch%d",ch), mem_onepar_func.c_str(), 0, 1.0e5);
    f_onepar[ch] -> SetParameters(120/mem_mppc_pixel/mem_led_width, 1.0);
    g_photon[ch] -> Fit(Form("f_onepar_ch%d",ch), "Q0", "", 5, mem_fit_max);
    f_twopar[ch] = new TF1(Form("f_twopar_ch%d",ch), mem_twopar_func.c_str(), 0, 1.0e5);
    f_twopar[ch] -> SetParameters(50/mem_mppc_pixel/mem_led_width, 100/mem_mppc_pixel/mem_led_width, 0.5, 0.5);
    f_twopar[ch] -> SetParLimits(2,0,1);
    f_twopar[ch] -> SetParLimits(3,0,1);
    g_photon[ch] -> Fit(Form("f_twopar_ch%d",ch), "Q0", "", 5, mem_fit_max);
    for(int cy=0; cy<mem_cycle; cy++){
      mem_onepar_residual[ch][cy] = mem_nobs_mean[ch][cy]/f_onepar[ch]->Eval(mem_nref_mean[ch][cy]);
      mem_twopar_residual[ch][cy] = mem_nobs_mean[ch][cy]/f_twopar[ch]->Eval(mem_nref_mean[ch][cy]);
      mem_onepar_residual_err[ch][cy] = mem_nobs_err[ch][cy]/f_onepar[ch]->Eval(mem_nref_mean[ch][cy]);
      mem_twopar_residual_err[ch][cy] = mem_nobs_err[ch][cy]/f_twopar[ch]->Eval(mem_nref_mean[ch][cy]);
    }
    g_onepar_resi[ch] = new TGraphErrors(mem_cycle, mem_nref_mean[ch], mem_onepar_residual[ch], mem_nref_err[ch], mem_onepar_residual_err[ch]);
    g_twopar_resi[ch] = new TGraphErrors(mem_cycle, mem_nref_mean[ch], mem_twopar_residual[ch], mem_nref_err[ch], mem_twopar_residual_err[ch]);
    mem_onepar_tau[ch] = (f_onepar[ch]->GetParameter(0))*mem_mppc_pixel*mem_led_width;
    //    mem_onepar_slope[ch] = f_onepar[ch]->GetParameter(1);
    mem_twopar_tau1[ch] = (f_twopar[ch]->GetParameter(0))*mem_mppc_pixel*mem_led_width;
    mem_twopar_tau2[ch] = (f_twopar[ch]->GetParameter(1))*mem_mppc_pixel*mem_led_width;
    mem_twopar_alpha[ch] = f_twopar[ch]->GetParameter(2);
    mem_twopar_beta[ch] = f_twopar[ch]->GetParameter(3);
    if(mem_twopar_tau1[ch] > mem_twopar_tau2[ch]){
      swap(mem_twopar_tau1[ch], mem_twopar_tau2[ch]);
      swap(mem_twopar_alpha[ch], mem_twopar_beta[ch]);
    }
  }  
  cout << "### analyzed result ###" << endl;
  for(int ch=1; ch<mem_readch; ch++){
    cout << "[ch: " << mem_ch[ch] << ", serial: " << mem_serial[ch] <<"]\n";
    cout << "one parameter: " << "tau=" << mem_onepar_tau[ch] << "ns\n";
    cout << "two parameter: " << "tau1=" << mem_twopar_tau1[ch] << "ns, tau2=" << mem_twopar_tau2[ch] << "ns (alpha=" << mem_twopar_alpha[ch] << ", beta=" << mem_twopar_beta[ch] << ")";
    cout << endl;
  }
}

void Analizer::saveTree(){
  string make_ofile_dir = "mkdir -p " + (string)mem_analyzed_file_path;
  system(make_ofile_dir.c_str());
  TString output_file_path = mem_analyzed_file_path + "/" + mem_output_file_name;
  TFile *ofile = new TFile(output_file_path, "recreate");
  TTree *header = new TTree("header", "header info");
  header -> Branch("fit_max", &mem_fit_max);
  header -> Branch("linearfit_num", &mem_linearfit_num);
  header -> Branch("cycle_num", &mem_cycle);
  header -> Branch("mppc_amp", &mem_mppc_amp);
  header -> Branch("mppc_pixel", &mem_mppc_pixel);
  header -> Branch("led_width", &mem_led_width);
  header -> Branch("e_charge", &mem_e_charge);
  header -> Fill();
  header -> Write();
  
  TTree *onepar_tree = new TTree("onepar_tree", "onepar_tree");
  TTree *twopar_tree = new TTree("twopar_tree", "twopar_tree");
  
  int b_ch=0, b_serial=0;
  double b_effgain=0, b_tau=0,/* b_slope=0,*/ b_tau1=0, b_tau2=0, b_alpha=0, b_beta=0;
  double b_nobs_mean[MAXCYCLE]={}, b_nref_mean[MAXCYCLE]={}, b_nobs_err[MAXCYCLE]={}, b_nref_err[MAXCYCLE]={}, b_onepar_residual[MAXCYCLE]={}, b_twopar_residual[MAXCYCLE]={};

  onepar_tree -> Branch("ch", &b_ch);
  onepar_tree -> Branch("serial", &b_serial);
  onepar_tree -> Branch("effgain", &b_effgain);
  onepar_tree -> Branch("tau", &b_tau);
  // onepar_tree -> Branch("slope", &b_slope);
  onepar_tree -> Branch("cycle", &mem_cycle);
  onepar_tree -> Branch("nobs_mean", b_nobs_mean, "nobs_mean[cycle]/D");
  onepar_tree -> Branch("nref_mean", b_nref_mean, "nref_mean[cycle]/D");
  onepar_tree -> Branch("nobs_err", b_nobs_err, "nobs_err[cycle]/D");
  onepar_tree -> Branch("nref_err", b_nref_err, "nref_err[cycle]/D");
  onepar_tree -> Branch("residual", b_onepar_residual, "residual[cycle]/D");
  
  twopar_tree -> Branch("ch", &b_ch);
  twopar_tree -> Branch("serial", &b_serial);
  twopar_tree -> Branch("effgain", &b_effgain);
  twopar_tree -> Branch("tau1", &b_tau1);
  twopar_tree -> Branch("tau2", &b_tau2);
  twopar_tree -> Branch("alpha", &b_alpha);
  twopar_tree -> Branch("beta", &b_beta);
  twopar_tree -> Branch("cycle", &mem_cycle);
  twopar_tree -> Branch("nobs_mean", b_nobs_mean, "nobs_mean[cycle]/D");
  twopar_tree -> Branch("nref_mean", b_nref_mean, "nref_mean[cycle]/D");
  twopar_tree -> Branch("nobs_err", b_nobs_err, "nobs_err[cycle]/D");
  twopar_tree -> Branch("nref_err", b_nref_err, "nref_err[cycle]/D");
  twopar_tree -> Branch("residual", b_twopar_residual, "residual[cycle]/D");
  
  for(int ch=1; ch<mem_readch; ch++){
    b_ch = mem_ch[ch];
    b_serial = mem_serial[ch];
    b_effgain = mem_effgain[ch];
    b_tau = mem_onepar_tau[ch];
    //    b_slope = mem_onepar_slope[ch];
    b_tau1 = mem_twopar_tau1[ch];
    b_tau2 = mem_twopar_tau2[ch];
    b_alpha = mem_twopar_alpha[ch];
    b_beta = mem_twopar_beta[ch];
    for(int cy=0; cy<mem_cycle; cy++){
      b_nobs_mean[cy] = mem_nobs_mean[ch][cy];
      b_nref_mean[cy] = mem_nref_mean[ch][cy];
      b_nobs_err[cy] = mem_nobs_err[ch][cy];
      b_nref_err[cy] = mem_nref_err[ch][cy];
      b_onepar_residual[cy] = mem_onepar_residual[ch][cy];
      b_twopar_residual[cy] = mem_twopar_residual[ch][cy];
    }
    onepar_tree -> Fill();
    twopar_tree -> Fill();
  }  
  onepar_tree -> Write();
  twopar_tree -> Write();

  TTree *tailer = new TTree("tailer", "analysis info");
  tailer -> Branch("target", &mem_target);
  tailer -> Branch("setup_num", &mem_setup_num);
  tailer -> Branch("input_dark", &mem_input_anadark_name);
  tailer -> Branch("input_integ", &mem_input_integ_name);
  tailer -> Branch("onepar_func", &mem_onepar_func);
  tailer -> Branch("twopar_func", &mem_twopar_func);
  tailer -> Branch("anadate", &mem_anadate);
  tailer -> Fill();
  tailer -> Write();
  
  ofile -> Close();
  cout <<"### " <<  mem_output_file_name << " is saved in (" << mem_analyzed_file_path <<") ###" << endl; 
}


void Analizer::makeCanvas(){
  TF1 *f_residual_base = new TF1("f_residual_base","1", 0, 1.0e5);
  f_residual_base -> SetLineColor(1);
  f_residual_base -> SetLineStyle(7);
  TPaveText *tp_anainfo[MAXCH], *tp_mppcinfo[MAXCH], *tp_oneparinfo[MAXCH], *tp_twoparinfo[MAXCH];  
  
  for(int ch=1; ch<mem_readch; ch++){
    canvas[ch] = new TCanvas("test","test",50,50,1500,900);
    canvas[ch] -> Divide(5,3,0.001,0.001);
    canvas[ch] -> cd(1);
    g_adc[ch] -> Draw("ap");
    g_adc[ch] -> SetTitle("ADC count;PMT ADC count;MPPC ADC count");
    f_linear[ch] -> Draw("same");
    canvas[ch] -> cd(2) -> DrawFrame(0,0,1.2*mem_integral_mean[0][mem_linearfit_num-1],1.2*mem_integral_mean[ch][mem_linearfit_num-1],"linear range ADC count;PMT ADC count;MPPC ADC count");
    g_adc[ch] -> Draw("p");
    f_linear[ch] -> Draw("same");
    canvas[ch] -> cd(3) -> DrawFrame(0,0.9,mem_integral_mean[0][mem_linearfit_num+1], 1.1, "linear fit residual(small region);PMT ADC count;(MPPC ADC count)/(linear fit)");
    g_linear_residual[ch] -> Draw("p");
    f_residual_base -> Draw("same");
    canvas[ch] -> cd(4);
    tp_mppcinfo[ch] = new TPaveText(0.1,0.1,0.9,0.9);
    tp_mppcinfo[ch] -> AddText("MPPC information");
    tp_mppcinfo[ch] -> AddText(Form("ch: %d, serial: %d",mem_ch[ch], mem_serial[ch]));
    tp_mppcinfo[ch] -> AddText(Form("effgain: %d",(int)mem_effgain[ch]));
    tp_mppcinfo[ch] -> AddText(Form("pixel: %dpix",(int)mem_mppc_pixel));
    tp_mppcinfo[ch] -> Draw();
    canvas[ch] -> cd(5);
    tp_anainfo[ch] = new TPaveText(0.1,0.1,0.9,0.9);
    tp_anainfo[ch] -> AddText("analysis information");
    tp_anainfo[ch] -> AddText(Form("target: %s (setup%s)",mem_target.Data(), mem_setup_num.Data()));
    tp_anainfo[ch] -> AddText(Form("anadate: %s",((TString)mem_anadate).Data()));
    tp_anainfo[ch] -> AddText(Form("max fit range: %d (#mu s^{-1})",(int)mem_fit_max));
    tp_anainfo[ch] -> AddText("input files:");
    tp_anainfo[ch] -> AddText(mem_input_anadark_name.Data());
    tp_anainfo[ch] -> AddText(mem_input_integ_name.Data());
    tp_anainfo[ch] -> Draw();
    canvas[ch] -> cd(6) -> DrawFrame(0,0,1.1*mem_nref_mean[ch][mem_cycle-1],1.1*mem_nobs_mean[ch][mem_cycle-1],"one parameter;Nref (#mu s^{-1});Nobs (#mu s^{-1})");
    g_photon[ch] -> Draw("p");
    f_onepar[ch] -> SetLineStyle(7);
    f_onepar[ch] -> SetLineColor(4);
    f_onepar[ch] -> DrawCopy("same");
    f_onepar[ch] -> SetLineStyle(1);
    f_onepar[ch] -> SetLineColor(2);
    f_onepar[ch] -> DrawF1(0,mem_fit_max,"same");
    canvas[ch] -> cd(7) -> DrawFrame(0,0.9,1.1*mem_nref_mean[ch][mem_cycle-1],1.1,"one parameter residual;Nref (#mu s^{-1});Data/Fit");
    g_onepar_resi[ch] -> Draw("p");
    f_residual_base -> Draw("same");
    canvas[ch] -> cd(8) -> DrawFrame(0,0.9,1.1*mem_nref_mean[ch][mem_linearfit_num+1],1.1,"one parameter residual (small region);Nref (#mu s^{-1});Data/Fit");
    g_onepar_resi[ch] -> Draw("p");
    f_residual_base -> Draw("same");
    canvas[ch] -> cd(9) -> DrawFrame(0,0.9,mem_fit_max,1.1,"one parameter residual (fit region);Nref (#mu s^{-1});Data/Fit");
    g_onepar_resi[ch] -> Draw("p");
    f_residual_base -> Draw("same");
    canvas[ch] -> cd(10);
    tp_oneparinfo[ch] = new TPaveText(0.1,0.1,0.9,0.9);
    tp_oneparinfo[ch] -> AddText("one parameter");
    tp_oneparinfo[ch] -> AddText(Form("tau=%.2fns",mem_onepar_tau[ch]));
    //    tp_oneparinfo[ch] -> AddText(Form("slope=%.3f",mem_onepar_slope[ch]));
    tp_oneparinfo[ch] -> Draw();
    canvas[ch] -> cd(11) -> DrawFrame(0,0,1.1*mem_nref_mean[ch][mem_cycle-1],1.1*mem_nobs_mean[ch][mem_cycle-1],"two parameter;Nref (#mu s^{-1});Nobs (#mu s^{-1})");
    g_photon[ch] -> Draw("p");
    f_twopar[ch] -> SetLineStyle(7);
    f_twopar[ch] -> SetLineColor(4);
    f_twopar[ch] -> DrawCopy("same");
    f_twopar[ch] -> SetLineStyle(1);
    f_twopar[ch] -> SetLineColor(2);
    f_twopar[ch] -> DrawF1(0,mem_fit_max,"same");
    canvas[ch] -> cd(12) -> DrawFrame(0,0.9,1.1*mem_nref_mean[ch][mem_cycle-1],1.1,"two parameter residual;Nref (#mu s^{-1});Data/Fit");
    g_twopar_resi[ch] -> Draw("p");
    f_residual_base -> Draw("same");
    canvas[ch] -> cd(13) -> DrawFrame(0,0.9,1.1*mem_nref_mean[ch][mem_linearfit_num+1],1.1,"two parameter residual (small region);Nref (#mu s^{-1});Data/Fit");
    g_twopar_resi[ch] -> Draw("p");
    f_residual_base -> Draw("same");
    canvas[ch] -> cd(14) -> DrawFrame(0,0.9,mem_fit_max,1.1,"two parameter residual (fit region);Nref (#mu s^{-1});Data/Fit");
    g_twopar_resi[ch] -> Draw("p");
    f_residual_base -> Draw("same");
    canvas[ch] -> cd(15);
    tp_twoparinfo[ch] = new TPaveText(0.1,0.1,0.9,0.9);
    tp_twoparinfo[ch] -> AddText("two parameter");
    tp_twoparinfo[ch] -> AddText(Form("tau1=%.1fns, tau2=%.0fns",mem_twopar_tau1[ch],mem_twopar_tau2[ch]));
    tp_twoparinfo[ch] -> AddText(Form("(alpha=%.3f, beta=%.3f)",mem_twopar_alpha[ch],mem_twopar_beta[ch]));
    tp_twoparinfo[ch] -> Draw();
  }
}

  void Analizer::drawCanvas(){
    TApplication app("app",0,0,0,0);

    for(int ch=1; ch<mem_readch; ch++){
      canvas[ch] -> DrawClone();
    }
    app.Run();
    /* string flag;
    while(1){
      cout << "input q to quit" << endl;
      cin>>flag;
      if(flag=="q") break;
    }
    */
  }
  
void Analizer::saveCanvas(){
  string make_canvas_dir = "mkdir -p " + (string)mem_analyzed_file_path + "/canvas/" + (string)mem_setup_num;
  system(make_canvas_dir.c_str());
  TString canvas_path = mem_analyzed_file_path + "/canvas/" + mem_setup_num + "/";
  for(int ch=1; ch<mem_readch; ch++){
    canvas[ch] -> Print(Form("%smaxfit%d_ch%d.pdf",canvas_path.Data(),(int)mem_fit_max,ch));
  }
}

int main(int argc, char* argv[]){
  if( argc<4 || argc>6 ){
    cout << "usage: ./recovery target (anadarkfile).root setupnum [maxfitrange] [anadate(1*****)]\n";
    cout << "target files are shown below\n" << endl;
    system("ls ~/DAQ/setup");
    return 0;
  }
    string target = string(argv[1]);
    string anadark_name = string(argv[2]);
    string setup_num = string(argv[3]);
    double fit_max;
    string anadate;
    if(argc==6){
      fit_max = atof(argv[4]); 
      anadate = string(argv[5]);
      cout << "maxfit: " << (int)fit_max << ", anadate: " << anadate << endl;
    }else if(argc==5){
      fit_max = atof(argv[4]); 
      cout << "maxfit: " << (int)fit_max << ", anadate: ********" << endl;
      anadate = "********";
    }else if(argc==4){
      cout << "maxfit: [no maxfit. max-cycle is used.], anadate: ********" << endl;
      fit_max = 1.0e5;
      anadate = "********";
    }
    gStyle -> SetMarkerStyle(20);
    gStyle -> SetTitleYOffset(1.3);
    Analizer *anal = new Analizer(target, setup_num, anadark_name, fit_max, anadate);
    if(argc>=6){
      anal -> saveTree();
      anal -> saveCanvas();
    }
    anal -> drawCanvas(); 
    delete anal; 
}
