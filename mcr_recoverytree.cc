void mcr_recoverytree(string inputfile){
  const TString RECANA = "modiv"; //org, two, modiv, twomodiv
  const int SETUP_NUM = 5;
  const int READCH = 1;

  string filename, extension;
  
  istringstream filestream(inputfile);
  getline(filestream,filename,'.');
  getline(filestream,extension,'.');

  TString INPUTFILE = inputfile;
  TString OUTPUTFILE = "tree_"+filename+".root";
  
  TFile *f = new TFile(OUTPUTFILE,"recreate");
  TTree *tree_recovery = new TTree("recovery","recovery time");
  

  ifstream finrectime(INPUTFILE);
  if(!finrectime){
    cout<<"#### can not find "<<INPUTFILE<<" ####"<<endl;
    return 1;
  }
  cout<<INPUTFILE<<" was found!!"<<endl;

  int ch, serial;
  float effgain, chisq;
  int setup_index=1;
  
  tree_recovery -> Branch("ch", &ch);
  tree_recovery -> Branch("serial", &serial);
  tree_recovery -> Branch("effgain", &effgain);

  
  string dummy;
  if(RECANA == "org" || RECANA == "modiv"){
    float tau, tauerr;
    tree_recovery -> Branch("tau", &tau);
    tree_recovery -> Branch("tauerr", &tauerr);
    tree_recovery -> Branch("chisq", &chisq);
    
    while(setup_index<=SETUP_NUM){
      while(finrectime >> dummy) if(dummy == "chisq") break;
      int ch_index=1;
      while(finrectime >> ch >> serial >> effgain >> tau >> tauerr >> chisq){
	tree_recovery -> Fill();
	ch_index++;
	if(ch_index > READCH) break;
      }
      setup_index++;
    }
    
  }else if(RECANA== "two" || RECANA == "twomodiv"){
    float alpha, alphaerr, beta, betaerr, tau1, tauerr1, tau2, tauerr2;
    tree_recovery -> Branch("alpha", &alpha);
    tree_recovery -> Branch("alphaerr", &alphaerr);
    tree_recovery -> Branch("beta", &beta);
    tree_recovery -> Branch("betaerr", &betaerr);
    tree_recovery -> Branch("tau1", &tau1);
    tree_recovery -> Branch("tauerr1", &tauerr1);
    tree_recovery -> Branch("tau2", &tau2);
    tree_recovery -> Branch("tauerr2", &tauerr2);
    tree_recovery -> Branch("chisq", &chisq);
    
    while(setup_index<=SETUP_NUM){
      while(finrectime >> dummy) if(dummy == "chisq") break;
      int ch_index=1;
      while(finrectime >> ch >> serial >> effgain >> alpha >> alphaerr >> beta >> betaerr >> tau1 >> tauerr1 >> tau2 >> tauerr2 >> chisq){
	//	cout<<ch<<" "<<serial<<" "<<effgain<<" "<<tau<<" "<<tauerr<<endl;
	tree_recovery -> Fill();
	ch_index++;
	if(ch_index > READCH) break;
      }
      setup_index++;
    }
  }
 
  tree_recovery -> Write();
  f -> Close();
  
}
