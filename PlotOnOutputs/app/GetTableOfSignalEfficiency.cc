#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>
//ROOT INCLUDES
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>

//#include <tree.hh>

const bool _debug = true;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

using namespace std;

int main( int argc, char* argv[] )
{

  ifstream file;// read input list directly
  std::string line;
  int n = 0;
  std::cout << "[Usage]: ./GetTableOfSignalEfficiency inputList analysisTag<Razor2015_76X,Razor2016_80X,Razor2017_92X> \n" << std::endl;
  srand(time(NULL));
  gROOT->Reset();
  gStyle->SetOptStat(0);
  //if ( argc == 3 )
  //{ 
        //file.open (argv[1], ios::in | ios::binary);
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return -1;
        } else {
                while(getline(file,line)) ++n;
                std::cout << "n = " << n << "\n" << std::endl;
        }
  
        std::string analysisTag = argv[2];
        if ( analysisTag == "" )
        {
                std::cerr << "[ERROR]: please provide the analysisTag<Razor2015_76X,Razor2016_80X,Razor2017_92X>" << std::endl;
                return -1;
        } 
  //}

  std::ifstream ifs( argv[1], std::ifstream::in );
  assert(ifs);

  TFile* fin; 
  TTree* tree; 
  vector<TFile*> fins; 
  vector<TTree*> trees; 
  double denominator;
  vector<double> den;
  //memset(denominator,0,n); 

  int i = 0;
  std::string process, rootFileName;
  //std::string labelProcess[n];
  vector<string> labelProcess;
  while ( ifs.good() ){
          ifs >> process >> rootFileName;
          if ( ifs.eof() ) continue;
          if ( process.find("#") != std::string::npos ) continue;
          if ( _debug ) std::cout << process << " " << rootFileName << std::endl;
          //labelProcess[i] = process;
          labelProcess.push_back(process);
          fin = new TFile( rootFileName.c_str(), "READ");
          fins.push_back(fin);
          if ( _debug ) std::cout << "[INFO]: file: " << rootFileName << " passed check\n\n"<< std::endl;

          //------------------------
          //Getting TTree and Histos
          //------------------------
          tree = (TTree*)fin->Get("HggRazorLeptons");
          trees.push_back(tree);
          denominator = ((TH1F*)fin->Get("NEvents"))->Integral();
          den.push_back(denominator);

          i++;
          std::cout << "i = " << i << "\n" << std::endl;
          //std::cout << "denominator = " << denominator << "\n" << std::endl;
  }


          /*TCanvas* c = new TCanvas( "c", "c", 800, 700 );
          TTree* tree0 = trees.at(0);
          TH1F* h0 = new TH1F("h0","",100,0,100000000);
          tree0->Draw("pTGammaGamma>>h0","mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40) && pho1Pt> 25. && pho2Pt>25.");
          double num = h0->Integral(); 
          std::cout << "num = " << num << "\n" << std::endl;
          TH1F* H0 = (TH1F*)gDirectory->Get("h0");
          double num = H0->Integral(); 
          std::cout << "i = " << i << "\n" << std::endl;
          std::cout << "num = " << num << "\n" << std::endl;*/
          
  //Put cuts
  vector<string> cut;
  vector<string> cut_name;
  vector<string> cut_name1;
  vector<string> cut_name2;
  vector<string> cut_name3;
  //baseline cut
  string baseline;
  /*
   * if(analysisTag!="Razor2017_92X")
          baseline = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40) && pho1Pt> 25. && pho2Pt>25.";
  else
          baseline = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5";
*/
          //baseline = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt>40. || pho2Pt>40.) && pho1Pt>25. && pho2Pt>25.";
  baseline = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1SC_Eta) <1.4442 && abs(pho2SC_Eta)<1.4442 && (pho1Pt/mGammaGamma>1./3. || pho2Pt/mGammaGamma>1./3.) && pho1Pt/mGammaGamma>1./4. && pho2Pt/mGammaGamma>1./4. && pho1R9>0.5 && pho2R9>0.5";

  cut.push_back(baseline);
  cut_name.push_back("Baseline");
  cut_name1.push_back("Baseline");
  cut_name2.push_back("Baseline");
  cut_name3.push_back("Baseline");
  
  //MR>150 cut
  string MR150 = baseline + " && MR>150"; 
  cut.push_back(MR150);
  cut_name.push_back("MR150");
  cut_name1.push_back("MR150");
  cut_name2.push_back("MR150");
  cut_name3.push_back("MR150");
  //box=0 Zmm
  string Zmm = MR150+ " && box == 0  ";
  //string Zmm = MR150+ " && box == 0  && lep1PassSelection > 1 && lep2PassSelection > 1";
  cut.push_back(Zmm);
  cut_name.push_back("Zmm");
  //box=1 Zee
  string Zee = MR150+ " && box == 1  ";
  //string Zee = MR150+ " && box == 1  && lep1PassSelection > 1 && lep2PassSelection > 1";
  cut.push_back(Zee);
  cut_name.push_back("Zee");
  cut_name1.push_back("TwoLeptons");
  cut_name3.push_back("TwoLeptons");
  //box=2 Emu
  string Emu = MR150+ " && box == 2  ";
  //string Emu = MR150+ " && box == 2  && lep1PassSelection > 1 && lep2PassSelection > 1";
  cut.push_back(Emu);
  cut_name.push_back("Emu");
  cut_name1.push_back("Emu");
  cut_name3.push_back("Emu");
  //box=3 OneMu
  string OneMu = MR150+ " && box == 3  ";
  //string OneMu = MR150+ " && box == 3  && lep1PassSelection > 1";
  cut.push_back(OneMu);
  cut_name.push_back("OneMu");
  //box=4 OneEle
  string OneEle = MR150+ " && box == 4  ";
  //string OneEle = MR150+ " && box == 4  && lep1PassSelection > 1";
  cut.push_back(OneEle);
  cut_name.push_back("OneEle");
  cut_name1.push_back("OneLepton");
  cut_name2.push_back("Leptons");
  cut_name3.push_back("OneLeptonHighPt");
  cut_name3.push_back("OneLeptonLowPt");
  //box=5 HighPt
  string HighPt = MR150+ " && box == 5";
  cut.push_back(HighPt);
  cut_name.push_back("HighPt");
  cut_name1.push_back("HighPt");
  cut_name3.push_back("HighPt");
  //box=6 Hbb
  string Hbb = MR150+ " && box == 6";
  cut.push_back(Hbb);
  cut_name.push_back("Hbb");
  //box=7 Zbb
  string Zbb = MR150+ " && box == 7";
  cut.push_back(Zbb);
  cut_name.push_back("Zbb");
  cut_name1.push_back("HZbb");
  cut_name3.push_back("HZbb");
  //box=8 HighRes
  string HighRes = MR150+ " && box == 8";
  cut.push_back(HighRes);
  cut_name.push_back("HighRes");
  //box=9 LowRes
  string LowRes = MR150+ " && box == 9";
  cut.push_back(LowRes);
  cut_name.push_back("LowRes");
  cut_name1.push_back("LowPt");
  cut_name2.push_back("Hadrons");
  cut_name3.push_back("LowPt");

/*  
  //Spilt OneLepton by pTGammaGamma cut at 110GeV
  //Remember to use cut_name3 in line 301
  cut.clear();
  cut.push_back(baseline);
  cut.push_back(MR150);
  string TwoLeptons = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40) && pho1Pt> 25. && pho2Pt>25. && MR > 150 && (box == 0 || box == 1) ";
  cut.push_back(TwoLeptons);
  cut.push_back(Emu);
  string OneLeptonHighPt = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40) && pho1Pt> 25. && pho2Pt>25. && MR > 150 && (box == 3 || box ==4) && pTGammaGamma>110.";
  string OneLeptonLowPt = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40) && pho1Pt> 25. && pho2Pt>25. && MR > 150 && (box == 3 || box == 4) && pTGammaGamma<110.";
  cut.push_back(OneLeptonHighPt);
  cut.push_back(OneLeptonLowPt);
  cut.push_back(HighPt);
  string HZbb = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40) && pho1Pt> 25. && pho2Pt>25. && MR > 150 && (box == 6 || box == 7) ";
  cut.push_back(HZbb);
  string LowPt = "mGammaGamma > 103. && mGammaGamma < 160. && pho1passIso == 1 && pho2passIso == 1 && pho1passEleVeto == 1 && pho2passEleVeto == 1 && abs(pho1Eta) <1.48 && abs(pho2Eta)<1.48 && (pho1Pt>40||pho2Pt>40) && pho1Pt> 25. && pho2Pt>25. && MR > 150 && (box == 8 || box == 9) ";
  cut.push_back(LowPt);
*/

  int m = cut.size();
  TH1F* h[n][m];
  double num[n][m];
  double ratio[n][m];
  for(int i = 0; i < n; i++)
  {
          TTree* thistree = trees.at(i);
          TH1F* hist = new TH1F("hist","",100,0,100000000);
          double thisdenominator = den.at(i);
          for(int j = 0; j < m; j++)
          {
                  string thisstring = cut.at(j);
                  thistree->Draw("pTGammaGamma>>hist",thisstring.c_str());
                  double numerator = hist->Integral();
                  num[i][j] = numerator; 
                  //std::cout << "numerator = " << numerator << "\n" << std::endl;
                  //std::cout << "num = " << num[i][j] << "\n" << std::endl;
                  ratio[i][j] = num[i][j]/thisdenominator; 
                  //std::cout << "ratio = " << ratio[i][j] << "\n" << std::endl;
          }

  }

  //leptons; inclusive
  //twoLeptons, Emu, OneLeptons; HighPt, Hbb/Zbb, LowPt(High/LowRes)
  //box 0/1,2,3/4;5,6/7,8/9
  //j = box+2 : 2/3,4....
  double ratio_1[n][8]; // Baseline, MR150, twoLeptons...
  double ratio_2[n][4]; // Baseline, MR150, leptons, inclusive
  double ratio_3[n][9]; // Spilit One Lepton category by pTGammaGamma cut at 110
  for(int i = 0; i < n; i++)
  {
          ratio_1[i][0] = ratio[i][0];
          ratio_1[i][1] = ratio[i][1];
          ratio_1[i][2] = ratio[i][2] + ratio[i][3];
          ratio_1[i][3] = ratio[i][4];
          ratio_1[i][4] = ratio[i][5] + ratio[i][6];
          ratio_1[i][5] = ratio[i][7];
          ratio_1[i][6] = ratio[i][8] + ratio[i][9];
          ratio_1[i][7] = ratio[i][10] + ratio[i][11];

          ratio_2[i][0] = ratio[i][0];
          ratio_2[i][1] = ratio[i][1];
          ratio_2[i][2] = ratio[i][2] + ratio[i][3] + ratio[i][4] + ratio[i][5] + ratio[i][6];
          ratio_2[i][3] = ratio[i][7] + ratio[i][8] + ratio[i][9] + ratio[i][10] + ratio[i][11];
  }
/*
  for(int i = 0; i < n; i++)
          for(int j = 0; j < 12; j++)
                  std::cout << "ratio 0 = " << ratio[i][j] << "\n" << std::endl;

  for(int i = 0; i < n; i++)
          for(int j = 0; j < 8; j++)
                  std::cout << "ratio 1 = " << ratio_1[i][j] << "\n" << std::endl;

  for(int i = 0; i < n; i++)
          for(int j = 0; j < 4; j++)
                  std::cout << "ratio 2 = " << ratio_2[i][j] << "\n" << std::endl;

*/

  //TABLE 1

  std::cout << "\n\n";
  std::cout << "\\begin{table*}[htb]\n\\begin{center}\n\\caption{Signal efficiency of all categories\\label{tab:categories}}\n\\def\\arraystretch{1.5}";
  std::cout << "\n\\begin{tabular}{|c|c|c|c|c|c|c|c|}\n\\hline\n \\multirow{}{}{Cut} & \\multicolumn{" << n << "}{c|}{SignalEfficiency \(\\%\)} \\\\";
  std::cout << "\n\\cline{2-" << n+1 << "}\n";
  //std::cout << "\nCut & HZ_127 & HZ_300 & HZ_600 & HZ_900 \\\\\n\\hline\n";
  
  string form = ""; 
  string form0 ;
  string form1 = " & ";
  string form3;
  string underscore = "_";
  string underscore1 = "\\_";
  for(int i = 0; i < n; i++)
  {
          form += form1;
          form3 = labelProcess.at(i); 
          form3.replace(2,1,underscore1);
          form3.replace(7,1,underscore1);
          form += form3; 
  }
  string form2 = " \\\\";
  form += form2;
  std::cout << form << std::endl;
  std::cout << "\\hline\n";

  string tmp0, tmp1;
  std::stringstream ss;
  double tmp = 0.;
  for(int j = 0; j < m; j++)
  {
          //tmp0 = cut_name3.at(j);
          tmp0 = cut_name.at(j);
          for(int i = 0; i < n; i++)
          {
                  tmp0 += form1;
                  tmp = ratio[i][j]*100;
                  ss<<tmp;
                  tmp1 = ss.str();
                  tmp0 += tmp1;
                  ss.str("");
          }
          tmp0 += form2;
          std::cout << tmp0 << std::endl;
  }

  std::cout << "\\hline\n\\end{tabular}\n\\end{center}\n\\end{table*}" << std::endl;
  std::cout << "\n\n";

  //TABLE 2

  std::cout << "\\begin{table*}[htb]\n\\begin{center}\n\\caption{Signal efficiency of subgroups\\label{tab:subgroups}}\n\\def\\arraystretch{1.5}";
  std::cout << "\n\\begin{tabular}{|c|c|c|c|c|c|c|c|}\n\\hline\n \\multirow{}{}{Cut} & \\multicolumn{" << n << "}{c|}{SignalEfficiency \(\\%\)} \\\\";
  std::cout << "\n\\cline{2-" << n+1 << "}\n";
  
  form = "";
  for(int i = 0; i < n; i++)
  {
          form += form1;
          form3 = labelProcess.at(i); 
          form3.replace(2,1,underscore1);
          form3.replace(7,1,underscore1);
          form += form3; 
  }
  form += form2;
  std::cout << form << std::endl;
  std::cout << "\\hline\n";

  for(int j = 0; j < 8; j++)
  {
          tmp0 = cut_name1.at(j);
          for(int i = 0; i < n; i++)
          {
                  tmp0 += form1;
                  tmp = ratio_1[i][j]*100;
                  ss<<tmp;
                  tmp1 = ss.str();
                  tmp0 += tmp1;
                  ss.str("");
          }
          tmp0 += form2;
          std::cout << tmp0 << std::endl;
  }

  std::cout << "\\hline\n\\end{tabular}\n\\end{center}\n\\end{table*}" << std::endl;
  std::cout << "\n\n";

  //TABLE 3

  std::cout << "\\begin{table*}[htb]\n\\begin{center}\n\\caption{Signal Efficiency of leptons and inclusive\\label{tab:2groups}}\n\\def\\arraystretch{1.5}";
  std::cout << "\n\\begin{tabular}{|c|c|c|c|c|c|c|c|}\n\\hline\n \\multirow{}{}{Cut} & \\multicolumn{" << n << "}{c|}{SignalEfficiency \(\\%\)} \\\\";
  std::cout << "\n\\cline{2-" << n+1 << "}\n";
  
  form = "";
  for(int i = 0; i < n; i++)
  {
          form += form1;
          form3 = labelProcess.at(i); 
          form3.replace(2,1,underscore1);
          form3.replace(7,1,underscore1);
          form += form3; 
          //form += labelProcess.at(i); 
  }
  form += form2;
  std::cout << form << std::endl;
  std::cout << "\\hline\n";

  for(int j = 0; j < 4; j++)
  {
          tmp0 = cut_name2.at(j);
          for(int i = 0; i < n; i++)
          {
                  tmp0 += form1;
                  double tmp = ratio_2[i][j]*100;
                  ss<<tmp;
                  tmp1 = ss.str();
                  tmp0 += tmp1;
                  ss.str("");
          }
          tmp0 += form2;
          std::cout << tmp0 << std::endl;
  }

  std::cout << "\\hline\n\\end{tabular}\n\\end{center}\n\\end{table*}" << std::endl;
  std::cout << "\n\n";




  
  return 0;
}
