#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TVectorF.h"
#include "TClonesArray.h"
#include "TTreeStream.h"
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using std::cout;
using std::setw;


void PrintMoments();

/*

cd /u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/
aliroot -l
.L toyMCEffCheck.C+
toyMCEffCheck()

*/
Int_t switchOffParticle=-100;   // default = -100 for 4 particles
Int_t fUsedSign = 0;
Int_t fSkipEvent=10;   // every fSkipEvent th event is defined distorted
Double_t fSkipPhi=18.;    // phi< 2pi/fSkipPhi is switched of

ULong64_t nEvents=10000000;
const Int_t nParticles=5;
const Int_t nSignBins =3;
const Int_t nMoments = 27;
Int_t elMean = 20, piMean=100, kaMean=20, prMean=30, deMean=8;
Int_t nTracksPerEventArr[nSignBins][nParticles]={0};
Double_t elParams[]={8.,1.5};
Double_t piParams[]={4.,1.};
Double_t kaParams[]={11.,1.5};
Double_t prParams[]={17.,1.6};
Double_t deParams[]={25.,1.8};

Float_t trackdEdx[1000]={0.};
Float_t trackPhi[1000]={0};
Int_t trackSign[1000]={0};


TH1D *hParticles[nSignBins][nParticles];
TH1D *hFirstMoms[nSignBins][nParticles];
TF1 *fParticles[nSignBins][nParticles];
UInt_t cutBit=0;
Int_t sign=0;
Int_t signBin=-100;
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
TClonesArray funcLineShapesCArr("TF1",50000);
enum momentType{kEl=0,kPi=1,kKa=2,kPr=3,kDe=4,
  kElEl=5,kPiPi=6,kKaKa=7,kPrPr=8,kDeDe=9,
  kElPi=10,kElKa=11,kElPr=12,kElDe=13,
  kPiKa=14,kPiPr=15,kPiDe=16,
  kKaPr=17,kKaDe=18,
  kPr3=19,kPr4=20,
  kEl3=21,kEl4=22,
  kDe3=23,kDe4=24,
  kPi3=25,kPi4=26,};
TString momNames[nMoments] = {"El","Pi","Ka","Pr","De",
  "El2","Pi2","Ka2","Pr2","De2",
  "ElPi","ElKa","ElPr","ElDe",
  "PiKa","PiPr","PiDe",
  "KaPr","KaDe",
  "Pr3","Pr4",
  "El3","El4",
  "De3","De4",
  "Pi3","Pi4"};
  enum momentTypeUnlikeSign {
    kPiPosPiNeg=0,
    kPiPosKaNeg=1,
    kPiPosPrNeg=2,
    kKaPosPiNeg=3,
    kKaPosKaNeg=4,
    kKaPosPrNeg=5,
    kPrPosPiNeg=6,
    kPrPosKaNeg=7,
    kPrPosPrNeg=8,
    kLaPosLaNeg=9,
    kChPosChNeg=10,
    kBaPosBaNeg=11,
  };
//
//
// ================================================================================================
//
//
void toyMCEffCheck(){

  TTreeSRedirector *outputData = new TTreeSRedirector(Form("TIden_DataTree_%d.root",fSkipEvent), "recreate");
  TTreeSRedirector *outputFits = new TTreeSRedirector(Form("TIden_LineShapes_%d.root",fSkipEvent), "recreate");
  funcLineShapesCArr.SetOwner(kTRUE);
  for (Int_t j=0;j<nSignBins;j++) {
    for (Int_t i=0;i<nParticles;i++) {
      hParticles[j][i] = new TH1D(Form("hist_%d_bin_%d",i,j),Form("hist_%d_bin_%d",i,j),1500,0.,50.);
      hFirstMoms[j][i] = new TH1D(Form("firstMom_%d_bin_%d",i,j),Form("firstMom_%d_bin_%d",i,j),1500,0.,50.);
      fParticles[j][i] = new TF1(Form("particle_%d_bin_%d",i,j),"gaus",0,50);
    }
  }
  if (fUsedSign==-1) signBin = 0;
  if (fUsedSign== 1) signBin = 1;
  if (fUsedSign== 0) signBin = 2;

  Double_t twoPi = 2*TMath::Pi();
  Double_t pidVar = 0;
  TRandom randomGen;
  TVectorF recMoments(nMoments);
  for(Int_t i=0;i<nMoments; i++) recMoments[i]=0.;

  for (ULong64_t ievent=0;ievent<nEvents;ievent++){  // event loop

    // rest particle counters
    for (Int_t i=0;i<nSignBins;i++){
      for (Int_t j=0;j<nParticles;j++){
        nTracksPerEventArr[i][j]=0;
      }
    }
    // reset track counter
    for (Int_t i=0;i<1000;i++) { trackdEdx[i]=0.; trackSign[i]=0; trackPhi[i]=0.; }
    //
    Float_t cent=randomGen.Uniform(0,10);
    Float_t phi = -100.;
    Int_t trCount=0;
    //
    Int_t nTracksPerEventsPerParticle[nParticles]={0};
    nTracksPerEventsPerParticle[kEl] = randomGen.Poisson(elMean);
    nTracksPerEventsPerParticle[kPi] = randomGen.Poisson(piMean);
    nTracksPerEventsPerParticle[kKa] = randomGen.Poisson(kaMean);
    nTracksPerEventsPerParticle[kPr] = randomGen.Poisson(prMean);
    nTracksPerEventsPerParticle[kDe] = randomGen.Poisson(deMean);
    //
    // Switch of some particles
    for (Int_t i=0;i<nParticles;i++) {
      if ( switchOffParticle>-1 && i>=switchOffParticle ) nTracksPerEventsPerParticle[i]=0;
    }
    //
    // Generate electron dEdx
    for (Int_t i=0; i<nTracksPerEventsPerParticle[kEl];i++){
      sign = (i%2==0) ? sign = 1 : sign = -1;
      phi=randomGen.Uniform(0,twoPi);
      if (fSkipEvent>0) {
        if ( i%fSkipEvent==0 && (phi<(twoPi/fSkipPhi)) ) continue;
      }
      trackdEdx[trCount] = randomGen.Gaus(elParams[0],elParams[1]);
      hParticles[2][kEl]->Fill(trackdEdx[trCount]);
      trackPhi[trCount]=phi;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kEl]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kEl]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kEl]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kEl]++; trackSign[trCount]=1; }
      nTracksPerEventArr[2][kEl]++;
      trCount++;
    }
    //
    // Generate pion dEdx
    for (Int_t i=0; i<nTracksPerEventsPerParticle[kPi];i++){
      sign = (i%2==0) ? sign = 1 : sign = -1;
      phi=randomGen.Uniform(0,twoPi);
      if (fSkipEvent>0) {
        if (i%fSkipEvent==0 && phi<(twoPi/fSkipPhi)) continue;
      }
      trackdEdx[trCount] = randomGen.Gaus(piParams[0],piParams[1]);
      hParticles[2][kPi]->Fill(trackdEdx[trCount]);
      trackPhi[trCount]=phi;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kPi]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kPi]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kPi]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kPi]++; trackSign[trCount]=1; }
      nTracksPerEventArr[2][kPi]++;
      trCount++;
    }
    //
    // Generate kaon dEdx
    for (Int_t i=0; i<nTracksPerEventsPerParticle[kKa];i++){
      sign = (i%2==0) ? sign = 1 : sign = -1;
      phi=randomGen.Uniform(0,twoPi);
      if (fSkipEvent>0) {
        if (i%fSkipEvent==0 && phi<(twoPi/fSkipPhi)) continue;
      }
      trackdEdx[trCount] = randomGen.Gaus(kaParams[0],kaParams[1]);
      hParticles[2][kKa]->Fill(trackdEdx[trCount]);
      trackPhi[trCount]=phi;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kKa]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kKa]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kKa]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kKa]++; trackSign[trCount]=1; }
      nTracksPerEventArr[2][kKa]++;
      trCount++;
    }
    //
    // Generate proton dEdx
    for (Int_t i=0; i<nTracksPerEventsPerParticle[kPr];i++){
      sign = (i%2==0) ? sign = 1 : sign = -1;
      phi=randomGen.Uniform(0,twoPi);
      if (fSkipEvent>0) {
        if (i%fSkipEvent==0 && phi<(twoPi/fSkipPhi)) continue;
      }
      trackdEdx[trCount] = randomGen.Gaus(prParams[0],prParams[1]);
      hParticles[2][kPr]->Fill(trackdEdx[trCount]);
      trackPhi[trCount]=phi;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kPr]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kPr]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kPr]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kPr]++; trackSign[trCount]=1; }
      nTracksPerEventArr[2][kPr]++;
      trCount++;
    }
    // Generate proton dEdx
    for (Int_t i=0; i<nTracksPerEventsPerParticle[kDe];i++){
      sign = (i%2==0) ? sign = 1 : sign = -1;
      phi=randomGen.Uniform(0,twoPi);
      if (fSkipEvent>0) {
        if (i%fSkipEvent==0 && phi<(twoPi/fSkipPhi)) continue;
      }
      trackdEdx[trCount] = randomGen.Gaus(deParams[0],deParams[1]);
      hParticles[2][kDe]->Fill(trackdEdx[trCount]);
      trackPhi[trCount]=phi;
      // neg --> 0, neutrol --> 2, pos --> 1
      if (sign==-1) {hParticles[0][kDe]->Fill(trackdEdx[trCount]); nTracksPerEventArr[0][kDe]++; trackSign[trCount]=-1;}
      if (sign== 1) {hParticles[1][kDe]->Fill(trackdEdx[trCount]); nTracksPerEventArr[1][kDe]++; trackSign[trCount]=1; }
      nTracksPerEventArr[2][kDe]++;
      trCount++;
    }
    //
    // Fill histograms first moms
    for (Int_t i=0;i<nParticles;i++){
      if ( !(switchOffParticle>-1 && i>=switchOffParticle) ) {
        hFirstMoms[0][i]->Fill(nTracksPerEventArr[0][i]);
        hFirstMoms[1][i]->Fill(nTracksPerEventArr[1][i]);
        hFirstMoms[2][i]->Fill(nTracksPerEventArr[2][i]);
      }
    }

    if(ievent%(nEvents/10)==0) {
      cout << ievent << "  " << trCount << " " <<  nTracksPerEventArr[signBin][kEl];
      cout <<  "  " << nTracksPerEventArr[signBin][kPi] ;
      cout <<  "  " << nTracksPerEventArr[signBin][kKa] ;
      cout <<  "  " << nTracksPerEventArr[signBin][kPr] ;
      cout <<  "  " << nTracksPerEventArr[signBin][kDe] << endl;
    }
    //
    // Fill data tree
    for (Int_t itr=0;itr<1000;itr++){
      outputData->GetFile()->cd();
      if(trackdEdx[itr]>0){
        *outputData << "tracks"
        << "gid="   << ievent <<
        "dEdx="     << trackdEdx[itr] <<
        "cutBit="   << cutBit <<
        "sign="     << trackSign[itr] <<
        "cent="     << cent <<
        "phi="      << trackPhi[itr] <<
        "\n";
      }
    }
    //
    // Calculate Moments
    recMoments[kEl]=Float_t(nTracksPerEventArr[signBin][kEl]);
    recMoments[kPi]=Float_t(nTracksPerEventArr[signBin][kPi]);
    recMoments[kKa]=Float_t(nTracksPerEventArr[signBin][kKa]);
    recMoments[kPr]=Float_t(nTracksPerEventArr[signBin][kPr]);
    recMoments[kDe]=Float_t(nTracksPerEventArr[signBin][kDe]);

    recMoments[kElEl]=recMoments[kEl]*recMoments[kEl];
    recMoments[kPiPi]=recMoments[kPi]*recMoments[kPi];
    recMoments[kKaKa]=recMoments[kKa]*recMoments[kKa];
    recMoments[kPrPr]=recMoments[kPr]*recMoments[kPr];
    recMoments[kDeDe]=recMoments[kDe]*recMoments[kDe];

    recMoments[kElPi]=recMoments[kEl]*recMoments[kPi];
    recMoments[kElKa]=recMoments[kEl]*recMoments[kKa];
    recMoments[kElPr]=recMoments[kEl]*recMoments[kPr];
    recMoments[kElDe]=recMoments[kEl]*recMoments[kDe];

    recMoments[kPiKa]=recMoments[kPi]*recMoments[kKa];
    recMoments[kPiPr]=recMoments[kPi]*recMoments[kPr];
    recMoments[kPiDe]=recMoments[kPi]*recMoments[kDe];

    recMoments[kKaPr]=recMoments[kKa]*recMoments[kPr];
    recMoments[kKaDe]=recMoments[kKa]*recMoments[kDe];

    recMoments[kPr3]=recMoments[kPr]*recMoments[kPr]*recMoments[kPr];
    recMoments[kEl3]=recMoments[kEl]*recMoments[kEl]*recMoments[kEl];
    recMoments[kDe3]=recMoments[kDe]*recMoments[kDe]*recMoments[kDe];
    recMoments[kPi3]=recMoments[kPi]*recMoments[kPi]*recMoments[kPi];

    recMoments[kPr4]=recMoments[kPr]*recMoments[kPr]*recMoments[kPr]*recMoments[kPr];
    recMoments[kEl4]=recMoments[kEl]*recMoments[kEl]*recMoments[kEl]*recMoments[kEl];
    recMoments[kDe4]=recMoments[kDe]*recMoments[kDe]*recMoments[kDe]*recMoments[kDe];
    recMoments[kPi4]=recMoments[kPi]*recMoments[kPi]*recMoments[kPi]*recMoments[kPi];
    //
    // Dump moments to tree
    outputData->GetFile()->cd();
    *outputData << "events"
    << "gid=" << ievent <<
    "moment.=" << &recMoments <<             // second moments for particle+antiparticle
    "\n";

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
      fParticles[j][i]->SetParameters(hParticles[j][i]->GetMean(),hParticles[j][i]->GetRMS());
      if (hParticles[j][i]) hParticles[j][i] -> Fit(fParticles[j][i],"QN");
      funcLineShapesCArr[objcounter] = (TF1*)fParticles[j][i];
      objcounter++;
    }
  }
  outputFits->GetFile()->cd();
  funcLineShapesCArr . Write("funcLineShapesCArr",TObject::kSingleKey);
  funcLineShapesCArr.Clear("C");
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
  // for (Int_t i=0;i<4;i++) fParticles[i] -> Write();
  delete outputFits;
  delete outputData;

  PrintMoments();

}

void PrintMoments(){

  cout << " ================================================================================== "<<endl;
  cout << " ================================================================================== "<<endl;
  TFile *f = new TFile(Form("TIden_DataTree_%d.root",fSkipEvent));
  TTree *tree = (TTree*)f->Get("events");
  TH1D *h = new TH1D("hGen","hGen",nMoments,0.,nMoments);
  TH1D *htemp[nMoments];
  for (Int_t i=1;i<nMoments+1;i++) {
    tree->Draw(Form("moment.fElements[%d]",i-1),"","goff");
    htemp[i-1] = (TH1D*)tree->GetHistogram()->Clone();
    htemp[i-1]->SetName(momNames[i-1]);
    cout << momNames[i-1] << "  " << htemp[i-1]->GetMean() << endl;
    h->SetBinContent(i,htemp[i-1]->GetMean());
  }
  cout << " ================================================================================== "<<endl;
  cout << " ================================================================================== "<<endl;
  for (Int_t i=1;i<nMoments+1;i++) h->GetXaxis()->SetBinLabel(i,momNames[i-1]);
  TFile *outFile = new TFile(Form("TIden_ToyMC_Gen_%d.root",fSkipEvent),"recreate");
  h->Write();
  for (Int_t i=0;i<nMoments;i++) htemp[i]->Write();
  outFile -> Close();
  delete outFile; //yeni eklave etdim.

}

void CalculateRatios(){

  /*

  cd /u/marsland/PHD/macros/marsland_EbyeRatios/schleching/TIdentity/TIdentity/test/
  aliroot -l
  .L toyMCEffCheck.C+
  CalculateRatios()

  */

  TFile *fIn[8];
  TH1D  *hIn[8];
  TH1D  *hRatio[8];

  TFile *fRef = new TFile("/lustre/nyx/alice/users/marsland/TIdenTOYMC/TIden_ToyMC_Gen_0.root");
  TH1D  *hRef = (TH1D*)fRef->Get("hGen");

  for (Int_t i=0;i<8;i++){
    fIn[i] = new TFile(Form("/lustre/nyx/alice/users/marsland/TIdenTOYMC/TIden_ToyMC_Gen_%d.root",i+3));
    hIn[i] = (TH1D*)fIn[i]->Get("hGen");
    hIn[i]->SetName(Form("hIn_%d",i));
  }

  for (Int_t i=0;i<8;i++) hRatio[i] = (TH1D*)hRef->Clone();
  for (Int_t i=0;i<8;i++) {
    hIn[i]->SetName(Form("ratio_%d",i));
    hIn[i]->SetLineColor(colors[i]);
    hIn[i]->SetLineWidth(3);
    hIn[i]->SetMarkerColor(colors[i]);
    hIn[i]->Divide(hRatio[i]);
    hIn[i]->GetYaxis()->SetRangeUser(0.90,1.01);
    hIn[i]->GetYaxis()->SetTitle("ratio to ref");
    for (Int_t j=1;j<nMoments+1;j++) hIn[i]->GetXaxis()->SetBinLabel(j,momNames[j-1]);
  }

  hIn[0]->Draw();
  for (Int_t i=1;i<8;i++) hIn[i]->Draw("same");


  TFile *outFile = new TFile("TIden_ToyMC_Ratios.root","recreate");
  for (Int_t i=0;i<8;i++) hIn[i]->Write();
  outFile -> Close();
  delete outFile; //yeni eklave etdim.



}
