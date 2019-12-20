#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <utility>
#include <TGraphAsymmErrors.h>
#include <TEfficiency.h>
#include <iostream> 
#include <iomanip> 


  TH1F* rebin(TH1F* hist, vector<double> xlowedges);
  TH1F* rebin(TH1F* hist, Int_t nbins);
  TH2F* rebin(TH2F* hist, vector<double> xlowedges, vector<double> ylowedges);
  TH3F* rebin(TH3F* hist, vector<double> xlowedges, vector<double> ylowedges, vector<double> zlowedges);
  // TGraphAsymmErrors* createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
  //                                                           string histname, 
  //                                                           vector<Double_t> bins, 
  //                                                           Double_t xlow = -99, Double_t xhigh = -99, 
  //                                                           Double_t ylow = -99, Double_t yhigh = -99); 
  // TH2F* createEfficiencyHist2D(TH2F* numerator, TH2F* denominator,
  //                                               std::string histname, 
  //                                               vector<Double_t> xbins, vector<Double_t> ybins );


//--------------------------------------------------------------------------------------------------
// Rebin for 1D hists
//--------------------------------------------------------------------------------------------------
TH1F* rebin(TH1F* hist, vector<double> xlowedges) {

  Double_t *xbins;
  Double_t *BinContent;
  Double_t *BinError;
 
  xbins = new Double_t[xlowedges.size()];
  for (UInt_t i=0;i<xlowedges.size();i++) {
    xbins[i] = xlowedges[i];
//     cout << "bins: " << xbins[i] << endl;
  }

//   cout << "Nbins: " << xlowedges.size() << endl;

  BinContent = new Double_t[xlowedges.size()+2];
  BinError = new Double_t[xlowedges.size()+2];
  for (UInt_t binx=0; binx < xlowedges.size()+2 ;++binx) {
      BinContent[binx] = 0;
      BinError[binx] = 0;
  }

  TH1F *rebinHist = new TH1F(hist->GetName(), hist->GetTitle(), xlowedges.size()-1, xbins);

//   cout << "new nbins : " << rebinHist->GetXaxis()->GetNbins() << endl;

  //refill
  for (UInt_t i=0;int(i)<hist->GetXaxis()->GetNbins()+2;++i) {
    Double_t x = hist->GetXaxis()->GetBinCenter(i);
    UInt_t xbin = 0;

//     cout << "old bin " << i << " : " << x << endl;
    //Find which x rebinned bin we are in
    for (UInt_t binx=0; binx < xlowedges.size()+1 ; ++binx) {
      if (binx == 0) { //underflow 
        if (x <= xlowedges[0]) {
          xbin = 0;
          break;
        }
      } else if (binx < xlowedges.size()) {
        if (x > xlowedges[binx-1] && x <= xlowedges[binx]) {
          xbin = binx;
          break;
        }
      } else { //overflow
        if (x > xlowedges[binx]) {
          xbin = binx;
          break;
        }
      }
    }
    BinContent[xbin] += hist->GetBinContent(i);
    BinError[xbin] += hist->GetBinError(i)*hist->GetBinError(i);
  }
  
  for (UInt_t binx=0; binx < xlowedges.size()+2 ;++binx) {
    rebinHist->SetBinContent(binx,BinContent[binx]);
    rebinHist->SetBinError(binx,TMath::Sqrt(BinError[binx]));     
  } 
  
  delete [] xbins;
  delete [] BinContent;
  delete [] BinError;
  return rebinHist;
}

//--------------------------------------------------------------------------------------------------
// Rebin histogram into NBins
//--------------------------------------------------------------------------------------------------
TH1F* rebin(TH1F* hist, Int_t nbins) {
  TH1F *result = hist;
  if (nbins > 0) {

//     cout << "nbins = " << nbins << endl;

    vector<Double_t> bins;
    for (int b=0; b<nbins+1; ++b) {
      Double_t binsize = (hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins()) - 
                          hist->GetXaxis()->GetBinLowEdge(1)) / nbins;
      bins.push_back(hist->GetXaxis()->GetBinLowEdge(1) + binsize * b);
//       cout << "binsize = " << binsize << " : " << hist->GetXaxis()->GetBinLowEdge(1) + binsize * b << endl;
    }
    result = rebin(hist, bins);
  }
  return result;
}

//--------------------------------------------------------------------------------------------------
// Rebin for 2D hists
//--------------------------------------------------------------------------------------------------
TH2F* rebin(TH2F* hist, vector<double> xlowedges, vector<double> ylowedges) {

  Double_t *xLow;
  Double_t *yLow;
  Double_t **BinContent;
  Double_t **BinError;

  xLow = new Double_t[xlowedges.size()];
  yLow = new Double_t[ylowedges.size()];
  BinContent = new Double_t*[xlowedges.size()+1];
  BinError = new Double_t*[xlowedges.size()+1];
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {    
    BinContent[binx] = new Double_t[ylowedges.size()+1];
    BinError[binx] = new Double_t[ylowedges.size()+1];
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
      BinContent[binx][biny] = 0;
      BinError[binx][biny] = 0;
    }
  }

  for (UInt_t i=0; i<xlowedges.size();++i)
    xLow[i] = xlowedges[i];
  for (UInt_t i=0; i<ylowedges.size();++i) {
    yLow[i] = ylowedges[i];
  }

  TH2F *rebinHist = new TH2F(hist->GetName(), hist->GetTitle(), xlowedges.size() - 1, xLow, ylowedges.size() - 1, yLow);

  //refill the histogram
  for (UInt_t i=0;Int_t(i)<=hist->GetXaxis()->GetNbins()+1;i++) {
    for (UInt_t j=0;Int_t(j)<=hist->GetYaxis()->GetNbins()+1;j++) {
      

      Double_t x = hist->GetXaxis()->GetBinCenter(i);
      Double_t y = hist->GetYaxis()->GetBinCenter(j);
      UInt_t xbin = 0;
      UInt_t ybin = 0;

      for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {
        if (binx == 0) { //underflow 
          if (x <= xlowedges[0]) {
            xbin = 0;
            break;
          }
        } else if (binx < xlowedges.size()) {
          if (x > xlowedges[binx-1] && x <= xlowedges[binx]) {
            xbin = binx;
            break;
          }
        } else { //overflow
          if (x > xlowedges[binx]) {
            xbin = binx;
            break;
          }
        }
      }

      for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
        if (biny == 0) { //underflow 
          if (y <= ylowedges[0]) {
            ybin = 0;
            break;
          }
        } else if (biny < ylowedges.size()) {
          if (y > ylowedges[biny-1] && y <= ylowedges[biny]) {
            ybin = biny;
            break;
          }
        } else { //overflow
          if (y > ylowedges[biny]) {
            ybin = biny;
            break;
          }
        }
      }
  
      BinContent[xbin][ybin] += hist->GetBinContent(i,j);
      BinError[xbin][ybin] += hist->GetBinError(i,j)*hist->GetBinError(i,j);     
    }
  }

  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
      rebinHist->SetBinContent(binx,biny,BinContent[binx][biny]);
      rebinHist->SetBinError(binx,biny,TMath::Sqrt(BinError[binx][biny])); 
    }
  }

  delete [] xLow;
  delete [] yLow;
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {    
    delete [] BinContent[binx];
    delete [] BinError[binx];
  }
  delete [] BinContent;
  delete [] BinError;

  return rebinHist;
}


//--------------------------------------------------------------------------------------------------
// Rebin for 3D hists
//--------------------------------------------------------------------------------------------------
TH3F* rebin(TH3F* hist, vector<double> xlowedges, vector<double> ylowedges, vector<double> zlowedges) {

  Double_t *xLow;
  Double_t *yLow;
  Double_t *zLow;
  Double_t ***BinContent;
  Double_t ***BinError;

  xLow = new Double_t[xlowedges.size()];
  yLow = new Double_t[ylowedges.size()];
  zLow = new Double_t[zlowedges.size()];
  BinContent = new Double_t**[xlowedges.size()+1];
  BinError = new Double_t**[xlowedges.size()+1];
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {    
    BinContent[binx] = new Double_t*[ylowedges.size()+1];
    BinError[binx] = new Double_t*[ylowedges.size()+1];
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
      BinContent[binx][biny] = new Double_t[zlowedges.size()+1];
      BinError[binx][biny] = new Double_t[zlowedges.size()+1];
      for (UInt_t binz=0; binz < zlowedges.size()+1 ;++binz) {
        BinContent[binx][biny][binz] = 0;
        BinError[binx][biny][binz] = 0;
      }
    }
  }

  for (UInt_t i=0; i<xlowedges.size();++i)
    xLow[i] = xlowedges[i];
  for (UInt_t i=0; i<ylowedges.size();++i) {
    yLow[i] = ylowedges[i];
  }
  for (UInt_t i=0; i<zlowedges.size();++i) {
    zLow[i] = zlowedges[i];
  }

  TH3F *rebinHist = new TH3F(hist->GetName(), hist->GetTitle(), xlowedges.size() - 1, xLow, ylowedges.size() - 1, yLow, zlowedges.size() - 1, zLow);

  //refill the histogram
  for (UInt_t i=0;Int_t(i)<=hist->GetXaxis()->GetNbins()+1;++i) {
    for (UInt_t j=0;Int_t(j)<=hist->GetYaxis()->GetNbins()+1;++j) {
      for (UInt_t k=0;Int_t(k)<=hist->GetZaxis()->GetNbins()+1;++k) {
      
        Double_t x = hist->GetXaxis()->GetBinCenter(i);
        Double_t y = hist->GetYaxis()->GetBinCenter(j);
        Double_t z = hist->GetZaxis()->GetBinCenter(k);
        UInt_t xbin = 0;
        UInt_t ybin = 0;
        UInt_t zbin = 0;
        
        for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {
          if (binx == 0) { //underflow 
            if (x <= xlowedges[0]) {
              xbin = 0;
              break;
            }
          } else if (binx < xlowedges.size()) {
            if (x > xlowedges[binx-1] && x <= xlowedges[binx]) {
              xbin = binx;
              break;
            }
          } else { //overflow
            if (x > xlowedges[binx]) {
              xbin = binx;
              break;
            }
          }
        }
        
        for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
          if (biny == 0) { //underflow 
            if (y <= ylowedges[0]) {
              ybin = 0;
              break;
            }
          } else if (biny < ylowedges.size()) {
            if (y > ylowedges[biny-1] && y <= ylowedges[biny]) {
              ybin = biny;
              break;
            }
          } else { //overflow
            if (y > ylowedges[biny]) {
              ybin = biny;
              break;
            }
          }
        }
        
        for (UInt_t binz=0; binz < zlowedges.size()+1 ;++binz) {
          if (binz == 0) { //underflow 
            if (z <= zlowedges[0]) {
              zbin = 0;
              break;
            }
          } else if (binz < zlowedges.size()) {
            if (z > zlowedges[binz-1] && z <= zlowedges[binz]) {
              zbin = binz;
              break;
            }
          } else { //overflow
            if (z > zlowedges[binz]) {
              zbin = binz;
              break;
            }
          }
        }
        
        BinContent[xbin][ybin][zbin] += hist->GetBinContent(i,j,k);
        BinError[xbin][ybin][zbin] += hist->GetBinError(i,j,k)*hist->GetBinError(i,j,k);
      }
    }
  }
  
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
      for (UInt_t binz=0; binz < zlowedges.size()+1 ;++binz) {
        rebinHist->SetBinContent(binx,biny,binz,BinContent[binx][biny][binz]);
        rebinHist->SetBinError(binx,biny,binz,TMath::Sqrt(BinError[binx][biny][binz])); 
      }
    }
  }
  
  delete [] xLow;
  delete [] yLow;
  delete [] zLow;
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {    
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {    
      delete [] BinContent[binx][biny];
      delete [] BinError[binx][biny];
    }
    delete [] BinContent[binx];
    delete [] BinError[binx];
  }
  delete [] BinContent;
  delete [] BinError;
  
  return rebinHist;
}


//--------------------------------------------------------------------------------------------------
// Create Efficiency Histogram. 
//--------------------------------------------------------------------------------------------------

/*
TGraphAsymmErrors* createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
                                         string histname, 
                                         vector<Double_t> bins, 
                                         Double_t xlow, Double_t xhigh, 
                                         Double_t ylow, Double_t yhigh,
                                         bool smallerr = false
                                         ) {
  TH1F *n = numerator;
  TH1F *d = denominator;  
  if (bins.size() > 0) {   
    n = rebin(numerator,bins);
    d = rebin(denominator,bins);
  }


  Int_t nbins = n->GetNbinsX();

  assert(nbins <= 200);
  Double_t x[200];
  Double_t y[200];
  Double_t xErr[200];
  Double_t yErrLow[200];
  Double_t yErrHigh[200];

  for (int i=0;i < 200; i++) {
    x[i] = 0;
    y[i] = 0;
    xErr[i] = 0;
    yErrLow[i] = 0;
    yErrHigh[i] = 0;
  }

  //don't take the overflow bins
  for (int b=0; b<nbins ; ++b) {

    x[b] = n->GetXaxis()->GetBinCenter(b+1);   
    xErr[b] = 0.0;

    Double_t ratio = 0;
    Double_t errLow = 0;
    Double_t errHigh = 0;     

    Double_t n1 = TMath::Nint(n->GetBinContent(b+1));
    Double_t n2 = TMath::Nint(d->GetBinContent(b+1));
    cout << "numerator: " << n1 << " and denominator: " << n2 << endl;
    if (n1 > n2) n1 = n2;

    if (n2>0) {
      ratio = n1/n2;
      if (ratio > 1) ratio = 1;
      errLow = ratio - TEfficiency::ClopperPearson((UInt_t)n2, (UInt_t)n1, 0.68269, kFALSE);
      errHigh = TEfficiency::ClopperPearson((UInt_t)n2, (UInt_t)n1, 0.68269, kTRUE) - ratio;
    }

//     cerr << " done bin " << b << " " << x[b] << " : " << n1 << "(" << n->GetBinContent(b+1) << ")" << " / " << n2 << "(" << d->GetBinContent(b+1) << ")" << " = " << ratio << " " << errLow << " " << errHigh << endl;
    y[b] = ratio;
    yErrLow[b] = errLow;
    yErrHigh[b] = errHigh;

    cout << "x value: " << x[b] << " and y value: " << y[b] << endl; 

  }

  //get rid of points with large relative errorbars
  if(smallerr) { 
    for (int b=0; b<nbins ; ++b) { 
      if((yErrHigh[b] >= .8*y[b]) | (yErrLow[b] >= .8*y[b])) { 
        yErrHigh[b] = 0.0; 
        yErrLow[b] = 0.0; 
        y[b] = 0.0; 
        xErr[b] = 0.0; 
        x[b] = 0.0; 
      }
      //if(x[b] > 160) { 
      //  yErrHigh[b] = 0.0; 
      //  yErrLow[b] = 0.0; 
      //  y[b] = 0.0; 
      //  xErr[b] = 0.0; 
      //  x[b] = 0.0; 
      //}
    }
  } 

  TGraphAsymmErrors *efficiency = new TGraphAsymmErrors(nbins, x, y, xErr, xErr, yErrLow,yErrHigh );
//  efficiency->SetName(histname.c_str());
//  efficiency->SetTitle(histname.c_str());
//  efficiency->GetXaxis()->SetTitle(n->GetXaxis()->GetTitle());
//  efficiency->GetYaxis()->SetTitle("Efficiency");

  if (yhigh != -99)
    efficiency->SetMaximum(yhigh);
  if (ylow != -99)
    efficiency->SetMinimum(ylow);
  if (xlow != -99 && xhigh != -99) 
    efficiency->GetXaxis()->SetRangeUser(xlow,xhigh);

//  efficiency->SetMarkerSize(1);
//  efficiency->SetLineWidth(2);

  return efficiency;
}
*/

TGraphAsymmErrors* createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
                                         string histname, 
                                         vector<Double_t> bins, 
                                         Double_t xlow, Double_t xhigh, 
                                         Double_t ylow, Double_t yhigh,
                                         bool smallerr = false
                                         ) {
  TH1F *n = numerator;
  TH1F *d = denominator;  
  if (bins.size() > 0) {   
    n = rebin(numerator,bins);
    d = rebin(denominator,bins);
  }

  Int_t nbins = n->GetNbinsX();

  assert(nbins <= 200);

  //define vectors to temporarily hold graph values
  vector<double> X; 
  vector<double> Y; 
  vector<double> Xerrlow; 
  vector<double> Xerrhigh; 
  vector<double> Yerrhigh;
  vector<double> Yerrlow; 
  
  //don't take the overflow bins
  for (int b=0; b<nbins ; ++b) {

    Double_t xtemp = n->GetXaxis()->GetBinCenter(b+1); 
    Double_t xerrlow =  n->GetXaxis()->GetBinCenter(b+1) - n->GetXaxis()->GetBinLowEdge(b+1);  
    Double_t xerrhigh = n->GetXaxis()->GetBinUpEdge(b+1) - n->GetXaxis()->GetBinCenter(b+1);  

    Double_t ratio = 0;
    Double_t errLow = 0;
    Double_t errHigh = 0;

    Double_t n1 = TMath::Nint(n->GetBinContent(b+1));
    Double_t n2 = TMath::Nint(d->GetBinContent(b+1));
    cout << "numerator: " << n1 << " and denominator: " << n2 << endl;
    if (n1 > n2) n1 = n2;

    if (n2>0) {
      ratio = n1/n2;
      if (ratio > 1) ratio = 1;
      errLow = ratio - TEfficiency::ClopperPearson((UInt_t)n2, (UInt_t)n1, 0.68269, kFALSE);
      errHigh = TEfficiency::ClopperPearson((UInt_t)n2, (UInt_t)n1, 0.68269, kTRUE) - ratio;
    }

    cout << " done bin " << b << " " << xtemp << " : " << n1 << "(" << n->GetBinContent(b+1) << ")" << " / " << n2 << "(" << d->GetBinContent(b+1) << ")" << " = " << ratio << " " << errLow << " " << errHigh << endl;
    Double_t ytemp = ratio;
    Double_t yerrlowtemp = errLow;
    Double_t yerrhightemp = errHigh;
    
    if(smallerr) { 
      if((yerrhightemp >= .5*ytemp) | (yerrlowtemp >= .5*ytemp)) continue; 
      //if(xtemp > 160) continue; 
    }

    cout << "x: " << xtemp << " and y: " << ytemp << endl;

    X.push_back(xtemp);  
    Y.push_back(ytemp); 
    Xerrlow.push_back(xerrlow);
    Xerrhigh.push_back(xerrhigh);
    Yerrlow.push_back(yerrlowtemp);
    Yerrhigh.push_back(yerrhightemp);
  
  }

  //count entries in vector
  int nindices = X.size();

  Double_t x[nindices];
  Double_t y[nindices];
  Double_t xErrLow[nindices];
  Double_t xErrHigh[nindices];
  Double_t yErrLow[nindices];
  Double_t yErrHigh[nindices];

  //fill array from vector
  for (int i=0;i < nindices; i++) {
    x[i] = X.at(i);
    y[i] = Y.at(i);
    xErrLow[i] = Xerrlow.at(i);
    xErrHigh[i] = Xerrhigh.at(i);
    yErrLow[i] = Yerrlow.at(i);
    yErrHigh[i] = Yerrhigh.at(i);
  }

  TGraphAsymmErrors *efficiency = new TGraphAsymmErrors(nbins, x, y, xErrLow, xErrHigh, yErrLow,yErrHigh );
//  efficiency->SetName(histname.c_str());
//  efficiency->SetTitle(histname.c_str());
//  efficiency->GetXaxis()->SetTitle(n->GetXaxis()->GetTitle());
//  efficiency->GetYaxis()->SetTitle("Efficiency");

  if (yhigh != -99)
    efficiency->SetMaximum(yhigh);
  if (ylow != -99)
    efficiency->SetMinimum(ylow);
  if (xlow != -99 && xhigh != -99) 
    efficiency->GetXaxis()->SetRangeUser(xlow,xhigh);

//  efficiency->SetMarkerSize(1);
//  efficiency->SetLineWidth(2);

  return efficiency;
}

//--------------------------------------------------------------------------------------------------
// Create Efficiency Graph. Use Clopper Pearson uncertainty intervals.
//--------------------------------------------------------------------------------------------------
  TH2F* createEfficiencyHist2D(TH2F* numerator, TH2F* denominator,
                               std::string histname, 
                               vector<Double_t> xbins, vector<Double_t> ybins ) {
    
  TH2F *n = numerator;
  TH2F *d = denominator;
  if (xbins.size() > 0 && ybins.size() > 0) {
    n = rebin(numerator,xbins,ybins);
    d = rebin(denominator,xbins,ybins);
  }
  
  TH2F *eff = (TH2F*)n->Clone(histname.c_str());
  assert(eff);

  for (int b=1; b<eff->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<eff->GetYaxis()->GetNbins()+1 ; ++c) {

      Double_t num = TMath::Nint(n->GetBinContent(b,c));
      Double_t den = TMath::Nint(d->GetBinContent(b,c));
      Double_t ratio = 0;
      Double_t errLow = 0;
      Double_t errHigh = 0;
      if (den > 0) {
        ratio = num / den;
        if (ratio > 1) ratio = 1;
        errLow = ratio - TEfficiency::ClopperPearson((UInt_t)den, (UInt_t)num, 0.68269, kFALSE);
        errHigh = TEfficiency::ClopperPearson((UInt_t)den, (UInt_t)num, 0.68269, kTRUE) - ratio;
      }
      
      eff->SetBinContent(b,c,ratio);
      eff->SetBinError(b,c,(errLow+errHigh)/2);

      cout << b << "," << c << " : " << num << " / " << den << " = " << ratio << " | " << errLow << " " << errHigh << " " << (errLow+errHigh)/2 << " : " << eff->GetBinError(b,c) << endl;

    }
  }

  return eff;
}
