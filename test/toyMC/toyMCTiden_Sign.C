#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TLine.h"
#include "TVectorF.h"
#include "TClonesArray.h"
#include "TTreeStream.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using std::cout;
using std::setw;


void PlotdEdx(Int_t sign);

/*

cd /u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/
aliroot -l
.L toyMCTiden_Sign.C+
toyMCTiden_Sign()
PlotdEdx(0)

*/
Int_t fUsedSign = 0;
ULong64_t nEvents=10000;
const Int_t nParticles=5;
const Int_t nSignBins =3;
const Int_t nMoments = 19;
const Int_t elMean  = 8, piMean =14, kaMean =8, prMean =10, deMean =6;
const Int_t elBMean = 6, piBMean=10, kaBMean=6, prBMean=8, deBMean=5;

Int_t nTracksPerEventArr[nSignBins][nParticles]={0};
Double_t elParams[]={8.,1.5};
Double_t piParams[]={4.,1.};
Double_t kaParams[]={11.,1.5};
Double_t prParams[]={17.,1.6};
Double_t deParams[]={25.,1.8};
const Int_t nMaxTracksPerEvent = 10000;
Float_t trackdEdx[nMaxTracksPerEvent]={0.};
Int_t trackSign[nMaxTracksPerEvent]={0};

TH1D *hParticles[nSignBins][nParticles];
TH1D *hFirstMoms[nSignBins][nParticles];
TF1 *fParticles[nSignBins][nParticles];
Int_t sign=0;
Int_t signBin=-100;
UInt_t cutBit=0;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kDe=4,
  kElEl=5,kPiPi=6,kKaKa=7,kPrPr=8,kDeDe=9,
  kElPi=10,kElKa=11,kElPr=12,kElDe=13,
  kPiKa=14,kPiPr=15,kPiDe=16,
  kKaPr=17,kKaDe=18,
};
//
//
// ================================================================================================
//
//
Double_t myfunction(Double_t *x, Double_t *par)
{
   Double_t f = par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
   return f;
}
//
void toyMCTiden_Sign()
{

  TTreeSRedirector *outputData = new TTreeSRedirector("DataTree_Sign.root", "recreate");
  TTreeSRedirector *outputFits = new TTreeSRedirector("LineShapes_Sign.root", "recreate");
  funcLineShapesCArr.SetOwner(kTRUE);
  for (Int_t j=0;j<nSignBins;j++) {
    for (Int_t i=0;i<nParticles;i++) {
      hParticles[j][i] = new TH1D(Form("hist_%d_bin_%d",i,j),Form("hist_%d_bin_%d",i,j),1500,0.,50.);
      hFirstMoms[j][i] = new TH1D(Form("firstMom_%d_bin_%d",i,j),Form("firstMom_%d_bin_%d",i,j),1500,0.,50.);
      fParticles[j][i] = new TF1(Form("particle_%d_bin_%d",i,j),"gaus",0,50);
      // fParticles[j][i] = new TF1(Form("particle_%d_bin_%d",i,j),"[0]*exp(-0.5*TMath::Power((x-[1])/[2],2))",0,50);
      // fParticles[j][i] = new TF1(Form("particle_%d_bin_%d",i,j),myfunction,0,50,3);
      // fParticles[j][i] = new TF1(Form("particle_%d_bin_%d",i,j),[&](double*x, double *p){ return p[0]*exp(-0.5*TMath::Power((x[0]-p[1])/p[2],2)); }, 0, 50, 3);
    }
  }

  if (fUsedSign==-1) signBin = 0;
  if (fUsedSign== 1) signBin = 1;
  if (fUsedSign== 0) signBin = 2;

  Double_t pidVar = 0;
  TRandom randomGen;

  for (ULong64_t ievent=0;ievent<nEvents;ievent++){  // event loop

    // rest particle counters
    for (Int_t i=0;i<nSignBins;i++){
      for (Int_t j=0;j<nParticles;j++){
        nTracksPerEventArr[i][j]=0;
      }
    }
    // reset track counter
    for (Int_t i=0;i<nMaxTracksPerEvent;i++) { trackdEdx[i]=0.; trackSign[i]=0; }
    //
    Float_t cent=randomGen.Uniform(0,10);
    Int_t trCount=0;
    nTracksPerEventArr[1][kEl] = randomGen.Poisson(elMean);
    nTracksPerEventArr[1][kPi] = randomGen.Poisson(piMean);
    nTracksPerEventArr[1][kKa] = randomGen.Poisson(kaMean);
    nTracksPerEventArr[1][kPr] = randomGen.Poisson(prMean);
    nTracksPerEventArr[1][kDe] = randomGen.Poisson(deMean);

    nTracksPerEventArr[0][kEl] = randomGen.Poisson(elBMean);
    nTracksPerEventArr[0][kPi] = randomGen.Poisson(piBMean);
    nTracksPerEventArr[0][kKa] = randomGen.Poisson(kaBMean);
    nTracksPerEventArr[0][kPr] = randomGen.Poisson(prBMean);
    nTracksPerEventArr[0][kDe] = randomGen.Poisson(deBMean);

    nTracksPerEventArr[2][kEl] = randomGen.Poisson(nTracksPerEventArr[0][kEl]+nTracksPerEventArr[1][kEl]);
    nTracksPerEventArr[2][kPi] = randomGen.Poisson(nTracksPerEventArr[0][kPi]+nTracksPerEventArr[1][kPi]);
    nTracksPerEventArr[2][kKa] = randomGen.Poisson(nTracksPerEventArr[0][kKa]+nTracksPerEventArr[1][kKa]);
    nTracksPerEventArr[2][kPr] = randomGen.Poisson(nTracksPerEventArr[0][kPr]+nTracksPerEventArr[1][kPr]);
    nTracksPerEventArr[2][kDe] = randomGen.Poisson(nTracksPerEventArr[0][kDe]+nTracksPerEventArr[1][kDe]);
    //
    // put extra correlation
    nTracksPerEventArr[0][kEl]=nTracksPerEventArr[0][kKa];
    nTracksPerEventArr[1][kEl]=nTracksPerEventArr[1][kKa];

    if(ievent%10000==0) {
      cout << " event = " << ievent << "  ----  " << trCount;
      cout <<  "  " << nTracksPerEventArr[1][kEl] <<  "  " << nTracksPerEventArr[0][kEl] ;
      cout <<  "  " << nTracksPerEventArr[1][kPi] <<  "  " << nTracksPerEventArr[0][kPi] ;
      cout <<  "  " << nTracksPerEventArr[1][kKa] <<  "  " << nTracksPerEventArr[0][kKa] ;
      cout <<  "  " << nTracksPerEventArr[1][kPr] <<  "  " << nTracksPerEventArr[0][kPr] ;
      cout <<  "  " << nTracksPerEventArr[1][kDe] <<  "  " << nTracksPerEventArr[0][kDe] ;
      cout << endl;
    }
    //
    // Generate particles
    for (Int_t ipart=0; ipart<nParticles; ipart++){

      for (Int_t i=0; i<nTracksPerEventArr[1][ipart];i++){

        if (ipart == kEl) {trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]); hParticles[1][ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPi) {trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]); hParticles[1][ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kKa) {trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]); hParticles[1][ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPr) {trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]); hParticles[1][ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kDe) {trackdEdx[trCount] = randomGen.Gaus(deParams[0],deParams[1]); hParticles[1][ipart]->Fill(trackdEdx[trCount]);}
        hParticles[2][ipart]->Fill(trackdEdx[trCount]);
        trackSign[trCount] = 1;
        trCount++;
      }

    }
    //
    // Generate anti-particles
    for (Int_t ipart=0; ipart<nParticles; ipart++){

      for (Int_t i=0; i<nTracksPerEventArr[0][ipart];i++){

        if (ipart == kEl) {trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]); hParticles[0][ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPi) {trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]); hParticles[0][ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kKa) {trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]); hParticles[0][ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kPr) {trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]); hParticles[0][ipart]->Fill(trackdEdx[trCount]);}
        if (ipart == kDe) {trackdEdx[trCount] = randomGen.Gaus(deParams[0],deParams[1]); hParticles[0][ipart]->Fill(trackdEdx[trCount]);}
        hParticles[2][ipart]->Fill(trackdEdx[trCount]);
        trackSign[trCount] = -1;
        trCount++;
      }

    }

    //
    // Fill histograms first moms
    for (Int_t i=0;i<nParticles;i++){
        hFirstMoms[0][i]->Fill(nTracksPerEventArr[0][i]);
        hFirstMoms[1][i]->Fill(nTracksPerEventArr[1][i]);
        hFirstMoms[2][i]->Fill(nTracksPerEventArr[2][i]);
    }

    //
    // Fill data tree
    for (Int_t itr=0;itr<10000;itr++){
      outputData->GetFile()->cd();
      if(trackdEdx[itr]!=0){
        *outputData << "tracks"
        << "gid="         << ievent <<
        "dEdx="           << trackdEdx[itr] <<
        "cutBit="         << cutBit <<
        "sign="           << trackSign[itr] <<
        "cent="           << cent <<
        "\n";
      }
    }
    //
  } // event loop ends
  //
  // dump line shape
  Int_t objcounter=0;
  for (Int_t j=0;j<nSignBins;j++) {
    for (Int_t i=0;i<nParticles;i++) {
      hParticles[j][i] -> SetLineColor(colors[i+1]);
      hFirstMoms[j][i] -> SetLineColor(colors[i+1]);
      fParticles[j][i] -> SetLineColor(colors[i+1]);
      fParticles[j][i] -> SetLineWidth(2);
      fParticles[j][i]->SetNpx(1000);
      fParticles[j][i]->SetParameters(hParticles[j][i]->GetBinContent(hParticles[j][i]->GetMaximumBin()),hParticles[j][i]->GetMean(),hParticles[j][i]->GetRMS());
      if (hParticles[j][i]) hParticles[j][i] -> Fit(fParticles[j][i],"QN");
      new (funcLineShapesCArr[objcounter++]) TF1(*fParticles[j][i]);
      // funcLineShapesCArr[objcounter++] = (TF1*)fParticles[j][i];
    }
  }
  outputFits->GetFile()->cd();
  funcLineShapesCArr . Write("funcLineShapesCArr",TObject::kSingleKey);
  funcLineShapesCArr . Clear("C");
  for (Int_t j=0;j<nSignBins;j++){
    for (Int_t i=0;i<nParticles;i++){
      hParticles[j][i] -> Write();
    }
  }
  for (Int_t j=0;j<nSignBins;j++){
    for (Int_t i=0;i<nParticles;i++){
      hFirstMoms[j][i] -> Write();
    }
  }
  delete outputFits;
  delete outputData;

}

void PlotdEdx(Int_t sign)
{


  TFile *f = new TFile("DataTree_Sign.root");
  TTree *tree = (TTree*)f->Get("tracks");
  TFile *g = new TFile("LineShapes_Sign.root");
  TClonesArray *cloneArrFunc = (TClonesArray*)g->Get("funcLineShapesCArr");
  TF1 *fShape[nParticles];

  Int_t binIndex=0;
  if (sign==-1) binIndex = 0;
  if (sign== 1) binIndex = 1;
  if (sign== 0) binIndex = 2;

  for (Int_t ipart = 0; ipart<nParticles; ipart++) {
    TString objName = Form("particle_%d_bin_%d",ipart,binIndex);
    fShape[ipart] = (TF1*)cloneArrFunc->FindObject(objName);
  }
  TString cutSign = (TMath::Abs(sign)==1) ? Form("sign==%d",sign) : "";
  tree->Draw("dEdx>>h(1500,0,50)",cutSign);
  for (Int_t ipart = 0; ipart<nParticles; ipart++){
    fShape[ipart]->Draw("same");
  }

}
