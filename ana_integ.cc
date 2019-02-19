#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TLeaf.h>
#include <TStyle.h>

using namespace std;

const int MAXCYCLE=50;
const int MAXSAMPLING=2048;

class Analizer{
public:
  Analizer(string filename, int max_cycle, const bool is_debug=false);

private:
  int m_MAXCYCLE = MAXCYCLE;
  int m_MAXSAMPLING = MAXSAMPLING;
  string m_ifilename;
  int m_max_cycle;
  bool m_is_debug;

  int m_max_ch;

  int m_integral_range_min = 1000, m_integral_range_max = 1200;
  int m_dark_range_min = 400, m_dark_range_max = 600;
  
  int m_integral_hist_min = -200, m_integral_hist_max = 200000;
  int m_dark_hist_min = -200, m_dark_hist_max = 200000;
  
  int m_ch;
  double m_integral_mean[MAXCYCLE] = {};
  double m_integral_error[MAXCYCLE] = {};
  double m_dark_mean[MAXCYCLE] = {}, m_dark_error[MAXCYCLE] = {};
  
  int m_cycle;
  int m_event;
  double m_integral_event;
  double m_dark_event;
  double m_baseline, m_baseline_sigma;
};

Analizer::Analizer(string filename, int max_cycle, const bool is_debug){
  
  m_ifilename = filename;
  m_max_cycle = max_cycle;
  m_is_debug = is_debug;
  if(m_is_debug){
    cout << "# Message: debug mode" << endl;
  }

  gStyle -> SetOptFit();

  TFile* ifile[m_max_cycle];
  TTree* rawwave_tree[m_max_cycle];
  
  cout << "//---- input files ----//" << endl;
  for(int cy=0; cy<m_max_cycle; cy++){
    ifile[cy] = new TFile(Form("%s_%02d.root",m_ifilename.c_str(),cy+1));
    cout<<Form("%s_%02d.root",m_ifilename.c_str(),cy+1)<<endl;
    rawwave_tree[cy] = (TTree*)ifile[cy]->Get("rawwave");
    m_max_ch = rawwave_tree[cy] -> GetNbranches();
  }
  cout << "// ---- end of files ---- //" << endl;
  
  TFile *ofile = new TFile(Form("sum_%s.root", m_ifilename.c_str()), "recreate");
  TTree *integral_header = new TTree("integral_header","integral_header");
  TTree *integral_tree = new TTree("integral_tree", "integral_tree");
  TTree *integral_event_tree = new TTree("integral_event_tree", "integral_event_tree");
  
  integral_header -> Branch("integral_range_min", &m_integral_range_min);
  integral_header -> Branch("integral_range_max", &m_integral_range_max);
  integral_header -> Branch("dark_range_min", &m_dark_range_min);
  integral_header -> Branch("dark_range_max", &m_dark_range_max);
  integral_header -> Branch("MAXCYCLE", &m_MAXCYCLE);
  integral_header -> Branch("MAXSAMPLING", &m_MAXSAMPLING);
  integral_header -> Fill();
  integral_header -> Write();
    
  integral_tree -> Branch("max_cycle", &m_max_cycle, "max_cycle/I");
  integral_tree -> Branch("ch", &m_ch);
  integral_tree -> Branch("integral_mean", m_integral_mean, "integral_mean[max_cycle]/D");
  integral_tree -> Branch("integral_error", m_integral_error, "integral_error[max_cycle]/D");
  integral_tree -> Branch("dark_mean", m_dark_mean, "dark_mean[max_cycle]/D");
  integral_tree -> Branch("dark_error", m_dark_error, "dark_error[max_cycle]/D");

  integral_event_tree -> Branch("cycle", &m_cycle);
  integral_event_tree -> Branch("event", &m_event);
  integral_event_tree -> Branch("ch", &m_ch);
  integral_event_tree -> Branch("integral_event", &m_integral_event);
  integral_event_tree -> Branch("dark_event", &m_dark_event);
  integral_event_tree -> Branch("baseline", &m_baseline);
  integral_event_tree -> Branch("baseline_sigma", &m_baseline_sigma);
    
  TCanvas *c_baseline = new TCanvas("c_baseline", "c_baseline",700,500);
  TCanvas *c_integral = new TCanvas("c_integral", "c_integral",700,500);
  TCanvas *c_dark = new TCanvas("c_dark", "c_dark",700,500);
  gSystem -> Unlink("debug_baseline.gif");
  gSystem -> Unlink("debug_integral.gif");
  gSystem -> Unlink("debug_dark.gif");
    
  
  for(int ch=0; ch<m_max_ch; ch++){
    m_ch = ch;
    cout << "### ch" << ch << " / " << m_max_ch-1 << " ###" << endl; 
    cout << "cycle: ";
    for(int cy=0; cy<m_max_cycle; cy++){
      m_cycle = cy;
      cout << cy << ", ";
      cout << flush;
      int max_event = rawwave_tree[cy] -> GetEntries();
      float input_wave[MAXSAMPLING] = {};
      rawwave_tree[cy] -> SetBranchAddress(Form("ch%d",ch), input_wave);
      int input_wave_length = rawwave_tree[cy] -> GetLeaf(Form("ch%d",ch)) -> GetLen();
      
      TH1D *h_integral = new TH1D("h_integral", 
				  "h_integral", 
				  (m_integral_hist_max-m_integral_hist_min)/10, 
				  m_integral_hist_min, 
				  m_integral_hist_max);
      TH1D *h_dark = new TH1D("h_dark", 
			      "h_dark", 
			      (m_dark_hist_max-m_dark_hist_min)/10, 
			      m_dark_hist_min, 
			      m_dark_hist_max);
      
      for(int ev=0; ev<max_event; ev++){
	rawwave_tree[cy] -> GetEntry(ev);
	m_event = ev;
	TH1D *h_baseline = new TH1D("h_baseline", "h_baseline", pow(2,14), -2250, 0);	
	for(int smp=0; smp<input_wave_length; smp++){
	  h_baseline -> Fill(input_wave[smp]);
	}
	TF1 *f_baseline = new TF1("f_baseline", "gaus", -2250, 0);
        double baseline_peak_bin = -2250.0 + h_baseline->GetBinWidth(0) * h_baseline -> GetMaximumBin();
	f_baseline -> SetParLimits(1,baseline_peak_bin-5, baseline_peak_bin+5);
	f_baseline -> SetParameter(2,2);
	h_baseline -> Fit("f_baseline", "Q", "", baseline_peak_bin-5, baseline_peak_bin+5);
	m_baseline = f_baseline -> GetParameter(1);
	m_baseline_sigma = f_baseline -> GetParameter(2);
	
	m_integral_event = 0;
	for(int smp=m_integral_range_min; smp<m_integral_range_max; smp++){
	  m_integral_event += m_baseline - (double)input_wave[smp];
	}
	m_dark_event = 0;
	for(int smp=m_dark_range_min; smp<m_dark_range_max; smp++){
	  m_dark_event += m_baseline - (double)input_wave[smp];
	}

	if(m_is_debug&&ev%100==0){
	  c_baseline -> cd() -> DrawFrame(m_baseline-20,0,m_baseline+20,600,Form("baseline: ch%d, cycle%d, event%d;ADC count; # of entries",ch,cy,ev));
	  h_baseline -> Draw("sames");
	  c_baseline -> Modified();
	  c_baseline -> Update();
	  c_baseline -> Print("debug_baseline.gif+");
	  }
	
	integral_event_tree -> Fill();
	if(h_integral->Fill(m_integral_event) == -1){
	  cout << "#Error: h_integral is overflowed or underflowed" << endl;
	  cout << "cycle=" << cy+1
	       << ", ch=" << ch
	       << ", event=" << ev
	       << ", baseline_peak_bin=" << baseline_peak_bin
	       << ", baseline=" << m_baseline 
	       << ", value=" << m_integral_event << endl;
	  exit(-1);
	}
	if(h_dark->Fill(m_dark_event) == -1){
	  cout << "#Error: h_dark_event is oveflowed or underflowed" << endl;
	  cout << "ch=" << ch
	       << ", event=" << ev
	       << ", value=" << m_dark_event << endl;
	  exit(-1);
	}
	delete h_baseline;
	delete f_baseline;
      }
      
      double integral_peak_bin = m_integral_hist_min + h_integral->GetBinWidth(0) * h_integral->GetMaximumBin();
      double dark_peak_bin = m_dark_hist_min + h_dark->GetBinWidth(0) * h_dark->GetMaximumBin();
      double hist_integral_sigma = h_integral->GetStdDev();
      double hist_dark_sigma = h_dark->GetStdDev();
      TF1 *f_integral = new TF1("f_integral", "gaus", m_integral_hist_min, m_integral_hist_max);
      TF1 *f_dark = new TF1("f_dark", "gaus", m_dark_hist_min, m_dark_hist_max);
      h_integral -> Fit("f_integral", "Q", "", integral_peak_bin-hist_integral_sigma*3, integral_peak_bin+hist_integral_sigma*3);
      h_dark -> Fit("f_dark", "Q", "", dark_peak_bin-hist_dark_sigma*3, dark_peak_bin+hist_dark_sigma*3);
      double integral_fit_mean= f_integral->GetParameter(1);
      double integral_fit_sigma = f_integral->GetParameter(2);
      double dark_fit_mean = f_dark->GetParameter(1);
      double dark_fit_sigma = f_dark->GetParameter(2);
      
      if(m_is_debug){
	c_integral -> cd() -> DrawFrame(-50,0,integral_peak_bin*2.0,h_integral->GetMaximum()*1.3,Form("Integral: ch%d, cycle%d; integral, # of entries",ch,cy));
	h_integral -> Draw("sames");
	c_integral -> Modified();
	c_integral -> Update();
	c_integral -> Print("debug_integral.gif+");
	c_dark -> cd() -> DrawFrame(-50,0,100,600,Form("Dark: ch%d, cycle%d; integral, # of entries",ch,cy));
	h_dark -> Draw("sames");
	c_dark -> Modified();
	c_dark -> Update();
	c_dark -> Print("debug_dark.gif+");
	}

      m_integral_mean[cy] = integral_fit_mean;
      m_integral_error[cy] = integral_fit_sigma/sqrt(max_event);
      m_dark_mean[cy] = dark_fit_mean;
      m_dark_error[cy] = dark_fit_sigma/sqrt(max_event);
      if(cy!=0 && m_integral_mean[cy]<m_integral_mean[cy-1]){
	cout << "\n# Warning: current integral value is smaller than the previous integral value." << endl;
	m_integral_mean[cy] = m_integral_mean[cy-1];
	m_integral_error[cy] = m_integral_error[cy-1];
	m_dark_mean[cy] = m_dark_mean[cy-1];
	m_dark_error[cy] = m_dark_error[cy-1];
      }
      delete h_integral;
      delete h_dark;
      delete f_integral;
      delete f_dark;
    }
    integral_tree -> Fill();
    cout << "\n";
  }
  integral_tree -> Write();
  integral_event_tree -> Write();
  
  for(int cy=0; cy<m_max_cycle; cy++){
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
  bool is_debug = false;

  Analizer *anal = new Analizer(ifilename, cycle, is_debug);

  delete anal;
}
