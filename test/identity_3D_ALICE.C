
//Author: Anar Rustamov
//identity method for third moments in case of 4 particles.

#include "TIdentity2D.h"
#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include "TStopwatch.h"
#include <fstream>
#include <sstream>
#include <cstdio>
#include <TF1.h>
#include <TH2F.h>
#include <math.h>
#include "TRandom3.h"
#include "TRandom1.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include <algorithm>
#include <math.h>
#include <numeric>
#include "TStyle.h"
#include "TVectorF.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TMatrixD.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TLorentzVector.h"

using namespace std;

Int_t numAllEvents;
Int_t numAllCutEvents;

TF1 *funProton, *funBackground, *funAProton;

Double_t evtNumber;
Double_t WK_sum, WP_sum, WAP_sum, WMult_sum;
Int_t multEv;

vector<int> sumMult;

vector<double> W2P_sum_vec;

vector<double> W3P_sum_vec;
vector<double> W3AP_sum_vec;
vector<double> W3K_sum_vec;

vector<double> W2K_sum_vec;
vector<double> W2AP_sum_vec;
vector<double> WPK_sum_vec;
vector<double> WPAP_sum_vec;
vector<double> WAPK_sum_vec;
vector<double> WAPPrK_sum_vec;

vector<double> WPr2AP_sum_vec;
vector<double> WPr2K_sum_vec;
vector<double> WAP2K_sum_vec;

vector<double> WAP2Pr_sum_vec;
vector<double> WK2Pr_sum_vec;
vector<double> WK2AP_sum_vec;

Double_t recP2_av, recAP2_av, recK2_av, recPAP_av, recPK_av, recAPK_av;
Double_t wP_P, wK_P, wAP_P, wP_P2, wK_P2, wAP_P2, wP_P3, wK_P3, wAP_P3, wP_K, wK_K, wAP_K, wP_K2, wK_K2, wAP_K2;
Double_t wP_K3, wK_K3, wAP_K3, wP_AP, wK_AP, wAP_AP, wP_AP2, wK_AP2, wAP_AP2, wP_AP3, wK_AP3, wAP_AP3;

///mixed
Double_t wPK_P, wPAP_P, wAPK_P, wPK_K, wPAP_K, wAPK_K, wPK_AP, wPAP_AP, wAPK_AP, wP2AP_P, wP2AP_AP, wP2AP_K, wAP2P_P, wAP2P_AP;
Double_t wAP2P_K, wP2K_P, wP2K_AP, wP2K_K, wK2P_P, wK2P_AP, wK2P_K, wAP2K_P, wAP2K_AP, wAP2K_K, wK2AP_P, wK2AP_AP, wK2AP_K, wAPPrK_Pr, wAPPrK_K, wAPPrK_AP;
//
TF1 *funP_P = NULL,    *funK_P = NULL,     *funAP_P = NULL,    *funP_P2 = NULL,   *funK_P2 = NULL,    *funAP_P2 = NULL,  *funP_P3 = NULL,   *funK_P3 = NULL,   *funAP_P3 = NULL;
TF1 *funP_K = NULL,    *funK_K = NULL,     *funAP_K = NULL,    *funP_K2 = NULL,   *funK_K2 = NULL,    *funAP_K2 = NULL,  *funP_K3 = NULL,   *funK_K3 = NULL,   *funAP_K3 = NULL;
TF1 *funP_AP = NULL,   *funK_AP = NULL,    *funAP_AP = NULL,   *funP_AP2 = NULL,  *funK_AP2 = NULL,   *funAP_AP2 = NULL, *funP_AP3 = NULL,  *funK_AP3 = NULL,  *funAP_AP3 = NULL;
TF1 *funPK_P = NULL,   *funPAP_P = NULL,   *funAPK_P = NULL,   *funPK_K = NULL,   *funPAP_K = NULL,   *funAPK_K = NULL,  *funPK_AP = NULL,  *funPAP_AP = NULL, *funAPK_AP = NULL;
TF1 *funP2AP_P = NULL, *funP2AP_AP = NULL, *funP2AP_K = NULL,  *funAP2P_P = NULL, *funAP2P_AP = NULL, *funAP2P_K = NULL, *funP2K_P = NULL,  *funP2K_AP = NULL, *funP2K_K = NULL;
TF1 *funK2P_P = NULL,  *funK2P_AP = NULL,  *funK2P_K = NULL,   *funAP2K_P = NULL, *funAP2K_AP = NULL, *funAP2K_K = NULL, *funK2AP_P = NULL, *funK2AP_AP = NULL;
TF1 *funK2AP_K = NULL, *funAPPrK_P = NULL, *funAPPrK_K = NULL, *funAPPrK_AP = NULL;
////////////////////////////

Int_t size_size = 6000;
Int_t fMyBin[3];

Float_t wmean[3][3];
Float_t wmean2[3][3];
Float_t wmean3[3][3];
Float_t wmix[3][3];

vector<double> WPM_sum_vec;
vector<double> WKM_sum_vec;
vector<double> WAPM_sum_vec;

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================
// ======= Modification Part =============================================================================

const Int_t fnEtaBins      = 16;  // MC: 8, Data:16
const Int_t fnCentBins     = 9;
const Int_t fnParticleBins = 3;
const Int_t fnMomBins      = 150;
const Int_t fnParticleBinsInTree = 8;
Int_t fMindEdx = -1020;
Int_t fMaxdEdx =  1020;
Double_t nSubSample = 25.;
Double_t fCentBinCenter = 0.;
//
const Int_t nBinsLineShape = 8160; // 6120
Int_t fTestNtracks = -1;
Bool_t lookUpTableForLine = kFALSE;
Int_t lookUpTableLineMode = 0; // 0 for hist and 1 for func
Int_t   fNthFitIteration = 6;
TString treeLineShapes   = "treeId";
TString treeIdentity     = "fIdenTree";   // for data "fIdenTree",    for MC "fIdenTreeMC"  , new data "tracks"

const Float_t fEtaRangeDown = -0.8;
const Float_t fEtaRangeUp   = 0.8;
const Float_t fMomRangeDown = 0.2;
const Float_t fMomRangeUp   = 3.2;
Float_t xCentBins[] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
const Int_t colors[]   = {kBlack, kRed+1 , kBlue+1, kGreen+3, kMagenta+1, kOrange-1,kCyan+2,kYellow+2, kRed, kGreen};
//
//
// Some Global declarations
TFile *fLineShapesLookUpTable = NULL;
TH1D *fChi2    = NULL;
TH1D *fcRows   = NULL;
TH1D *fhPtot   = NULL;
TH1D *fhEta    = NULL;
TH1D *fhCent   = NULL;
static TH1D *****hParticles;
static TF1 *****fParticles;
TTree *momTree = NULL;
TH1D **hDedxDebug;
TH1D **fHistWs;
TH1D **fHistOmegas;
//
Char_t  inputfileNameDataTree[255];     //file name of tree
Char_t  inputfileNameLineShapes[255];   // file name for fit functions
Int_t   fSubsample=-100;
Float_t fCentInput=-100;
Float_t fEtaDownInput=-100;
Float_t fEtaUpInput=-100;
Float_t fpDownInput=-100;
Float_t fpUpInput=-100;
//
//
// TString fileNameDataTree = "/home/marsland/Desktop/TIdentityCMake_Dev/test/inputALICE/SubSample_cent0.00_ss2.root";
// TString fileNameLineShapes = "/home/marsland/Desktop/TIdentityCMake_Dev/test/inputALICE/FitResults_piS0.7_kaS0.5_prS0.5_pikaprKauto.root";
TString fileNameDataTree ="";
TString fileNameLineShapes = "";
Int_t fCentInputBin = 0;
Int_t fpDownBin     = 0;
Int_t fpUpBin       = 0;
Int_t fEtaDownBin   = 0;
Int_t fEtaUpBin     = 0;
//
//
Double_t nEvents = 0;
Double_t nnorm   = 1.;
Int_t fEtaBin, fCentBin, fMomBin;
UInt_t fCutBit;
TVectorF *fIntegrals;
TVectorF *fMoments1st,  *fMoments2nd, *fMoments3rd;
//
//
TString fitFunctionGenGausStr = "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
TTree *treeLookUp = NULL;
TFile *outFile = NULL;
TF1   *fgengaus = 0;
TStopwatch timer;
//
Double_t fAmpArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBinsInTree]; // fAmpArr[fEtaBin][fCentBin][fMomBin][particleType]
Double_t fMeanArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBinsInTree];
Double_t fSigmaArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBinsInTree];
Double_t fSkewArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBinsInTree];
Double_t fKurtosisArr[fnEtaBins][fnCentBins][fnMomBins][fnParticleBinsInTree];

Int_t wrongCount = 0;
Float_t meanBackground[fnEtaBins][fnCentBins][fnMomBins];
Float_t meanProton[fnEtaBins][fnCentBins][fnMomBins];
Float_t meanAntiProton[fnEtaBins][fnCentBins][fnMomBins];
Float_t meanElectron[fnEtaBins][fnCentBins][fnMomBins];
Int_t fUsedBins[fnEtaBins][fnCentBins][fnMomBins];

enum momentType
{
  kEl=0,
  kPi=1,
  kKa=2,
  kPr=3,
  kBEl=4,
  kBPi=5,
  kBKa=6,
  kBPr=7,
};

enum parType
{
  kPR=0,
  kAP=1,
  kBG=2,
};


enum trackCutBit
{
  kCRows60=0,
  kCRows80=1,
  kCRows100=2,
  kChi2TPC3=3,
  kChi2TPC4=4,
  kChi2TPC5=5,
  kDCAXYSmall=6,
  kDCAXY=7,
  kDCAXYLarge=8,
  kVZSmall=9,
  kVZ=10,
  kVZLarge=11,
  kEventVertexZSmall=12,
  kEventVertexZ=13,
  kEventVertexZLarge=14,
  kClusterRequirementITS=15,
  kNewITSCut=16,
};

Double_t fitFunctionGenGaus(Double_t *x, Double_t *par)
{
  //
  // Generalised gauss function --> "[0]*exp(-(TMath::Abs(x-[1])/[2])**[3])*(1+TMath::Erf([4]*(x-[1])/[2]/TMath::Sqrt(2)))";
  // Skew-normal distribution   --> "[0]*exp(-0.5*((x-[1])/[2])**2)*(1+TMath::Erf([3]*(x-[1])/[2]/TMath::Sqrt(2)))";
  // par[0] --> Amplitude
  // par[1] --> Mean
  // par[2] --> Sigma
  // par[3] --> Kurtosis
  // par[4] --> Skewness
  //
  static Double_t sqrt2 = TMath::Sqrt(2);
  Double_t fun = par[0]*exp(-TMath::Power((TMath::Abs(x[0]-par[1])/par[2]),par[3]))
  *(1+TMath::Erf(par[4]*(x[0]-par[1])/par[2]/sqrt2));
  return fun;
}

void InitializeObjects()
{
  //
  fChi2   = new TH1D("fChi2" ,"fChi2",200 ,0, 10 );
  fcRows  = new TH1D("fcRows","fcRows",200 ,0, 200 );
  fhEta   = new TH1D("fhEta" ,"Eta Bins"       ,fnEtaBins ,fEtaRangeDown, fEtaRangeUp );
  fhPtot  = new TH1D("fhPtot","Momentum Bins"  ,fnMomBins ,fMomRangeDown, fMomRangeUp );
  fhCent  = new TH1D("fhCent","Centrality Bins",fnCentBins ,xCentBins );
  //
  Int_t min = fhPtot -> FindBin(fpDownInput   + 0.0000001) - 1;
  Int_t max = fhPtot -> FindBin(fpUpInput     - 0.0000001) - 1;
  Int_t nbins = max-min;
  hDedxDebug = new TH1D *[nbins];
  for (Int_t i = 0; i < nbins; i++) {
    hDedxDebug[i] = new TH1D(Form("hDedx_%d", i),Form("hDedx_%d", i),nBinsLineShape,fMindEdx,fMaxdEdx);
  }
  //
  fHistWs = new TH1D *[3];
  fHistOmegas = new TH1D *[3];
  for (Int_t i = 0; i < 3; i++)
  {
    fHistWs[i] = new TH1D(Form("hW_%d", i), Form("hW_%d", i), nBinsLineShape, 0., nBinsLineShape);
    fHistOmegas[i] = new TH1D(Form("hOmega_%d", i), Form("hOmega_%d", i), 100, 0., 1.);
  }
  //
  fParticles = new TF1 ****[fnEtaBins];
  for (Int_t i = 0; i<fnEtaBins; i++){
    fParticles[i] = new TF1***[fnCentBins];
    for (Int_t j = 0; j<fnCentBins; j++){
      fParticles[i][j] = new TF1**[fnMomBins];
      for (Int_t k = 0; k<fnMomBins; k++){
        fParticles[i][j][k] = new TF1*[fnParticleBins];
        for (Int_t ii = 0; ii<fnParticleBins; ii++){
          fParticles[i][j][k][ii] = NULL;
        }
      }
    }
  }
  hParticles = new TH1D ****[fnEtaBins];
  for (Int_t i = 0; i<fnEtaBins; i++){
    hParticles[i] = new TH1D***[fnCentBins];
    for (Int_t j = 0; j<fnCentBins; j++){
      hParticles[i][j] = new TH1D**[fnMomBins];
      for (Int_t k = 0; k<fnMomBins; k++){
        hParticles[i][j][k] = new TH1D*[fnParticleBins];
        for (Int_t ii = 0; ii<fnParticleBins; ii++){
          hParticles[i][j][k][ii] = NULL;
        }
      }
    }
  }
  //
  //
  for( Int_t i = 0; i < fnEtaBins; i++ ){
    for( Int_t j = 0; j < fnCentBins; j++ ){
      for( Int_t k = 0; k < fnMomBins; k++ ){
        fUsedBins[i][j][k] = -1;
      }
    }
  }
  //
  // initialize counters
  fIntegrals    = new TVectorF(3); // 0, 1, 2
  fMoments1st   = new TVectorF(3); // 0, 1, 2
  fMoments2nd   = new TVectorF(6); // 00, 11, 22, 01, 02, 1,2
  fMoments3rd   = new TVectorF(10); // 000, 111, 222, 001, 002, 110, 112, 220, 221, 123
  //
  for(Int_t i=0;i<3; i++){
    (*fIntegrals)[i]=0.;
    (*fMoments1st)[i]=0.;
  }
  //
  for(Int_t i=0;i<6; i++) (*fMoments2nd)[i]=0.;
  for(Int_t i=0;i<10; i++) (*fMoments3rd)[i]=0.;
  //
  // initialize output tree
  //
  momTree = new TTree("momTree","momTree");
  momTree -> Branch("nEvents",&nEvents);
  momTree -> Branch("nnorm",&nnorm);
  momTree -> Branch("subsample",&fSubsample);
  momTree -> Branch("pDown",&fpDownInput);
  momTree -> Branch("pUp",&fpUpInput);
  momTree -> Branch("etaDown",&fEtaDownInput);
  momTree -> Branch("etaUp",&fEtaUpInput);
  momTree -> Branch("centBinIndex",&fCentInputBin);
  momTree -> Branch("centBin",&fCentBinCenter);
  momTree -> Branch("fIntegrals",&fIntegrals);
  momTree -> Branch("fMoments1st",&fMoments1st);
  momTree -> Branch("fMoments2nd",&fMoments2nd);
  momTree -> Branch("fMoments3rd",&fMoments3rd);

}

TF1 *MergeOtherTF1s(Int_t i, Int_t j, Int_t k)
{

  TF1 *fEl = NULL,  *fPi = NULL,  *fKa = NULL;
  TF1 *fBEl = NULL, *fBPi = NULL, *fBKa = NULL;

  fEl = new TF1("fEl",fitFunctionGenGausStr,fMindEdx,fMaxdEdx);
  fEl->FixParameter(0,fAmpArr[i][j][k][kEl]);
  fEl->FixParameter(1,fMeanArr[i][j][k][kEl]);
  fEl->FixParameter(2,fSigmaArr[i][j][k][kEl]);
  fEl->FixParameter(3,fKurtosisArr[i][j][k][kEl]);
  fEl->FixParameter(4,fSkewArr[i][j][k][kEl]);
  fBEl = new TF1("fBEl",fitFunctionGenGausStr,fMindEdx,fMaxdEdx);
  fBEl->FixParameter(0,fAmpArr[i][j][k][kBEl]);
  fBEl->FixParameter(1,fMeanArr[i][j][k][kBEl]);
  fBEl->FixParameter(2,fSigmaArr[i][j][k][kBEl]);
  fBEl->FixParameter(3,fKurtosisArr[i][j][k][kBEl]);
  fBEl->FixParameter(4,fSkewArr[i][j][k][kBEl]);
  //
  fPi = new TF1("fPi",fitFunctionGenGausStr,fMindEdx,fMaxdEdx);
  fPi->FixParameter(0,fAmpArr[i][j][k][kPi]);
  fPi->FixParameter(1,fMeanArr[i][j][k][kPi]);
  fPi->FixParameter(2,fSigmaArr[i][j][k][kPi]);
  fPi->FixParameter(3,fKurtosisArr[i][j][k][kPi]);
  fPi->FixParameter(4,fSkewArr[i][j][k][kPi]);
  fBPi = new TF1("fBPi",fitFunctionGenGausStr,fMindEdx,fMaxdEdx);
  fBPi->FixParameter(0,fAmpArr[i][j][k][kBPi]);
  fBPi->FixParameter(1,fMeanArr[i][j][k][kBPi]);
  fBPi->FixParameter(2,fSigmaArr[i][j][k][kBPi]);
  fBPi->FixParameter(3,fKurtosisArr[i][j][k][kBPi]);
  fBPi->FixParameter(4,fSkewArr[i][j][k][kBPi]);
  //
  fKa = new TF1("fKa",fitFunctionGenGausStr,fMindEdx,fMaxdEdx);
  fKa->FixParameter(0,fAmpArr[i][j][k][kKa]);
  fKa->FixParameter(1,fMeanArr[i][j][k][kKa]);
  fKa->FixParameter(2,fSigmaArr[i][j][k][kKa]);
  fKa->FixParameter(3,fKurtosisArr[i][j][k][kKa]);
  fKa->FixParameter(4,fSkewArr[i][j][k][kKa]);
  fBKa = new TF1("fBKa",fitFunctionGenGausStr,fMindEdx,fMaxdEdx);
  fBKa->FixParameter(0,fAmpArr[i][j][k][kBKa]);
  fBKa->FixParameter(1,fMeanArr[i][j][k][kBKa]);
  fBKa->FixParameter(2,fSigmaArr[i][j][k][kBKa]);
  fBKa->FixParameter(3,fKurtosisArr[i][j][k][kBKa]);
  fBKa->FixParameter(4,fSkewArr[i][j][k][kBKa]);

  TF1 *fSumAll = new TF1(Form("particle_2_bin_%d_bin_%d_bin_%d",i,j,k),"fPi+fBPi+fEl+fBEl+fKa+fBKa",fMindEdx,fMaxdEdx);
  return fSumAll;

}

void ReadFitParamsFromTree(TString paramTreeName, Int_t fitIter)
{
  //
  // Read the fit parameters from the paramtree file "ParamTree.root" which comes from PID FITS
  //
  // Open the lookup table and initialize array
  fLineShapesLookUpTable = new TFile(paramTreeName);
  cout << " ReadFitParamsFromTree.Info: Read line shapes from ttree " << endl;
  treeLookUp = (TTree*)fLineShapesLookUpTable->Get(treeLineShapes);
  Int_t nent = treeLookUp -> GetEntries();
  cout << nent <<  " ReadFitParamsFromTree.Info: tree is fine go ahead " << endl;

  Int_t sign       = 0;
  Int_t it         = 0;
  Int_t sl         = 0;
  Int_t myBin[3]  = {0};  // 0; eta, 1;cent, 2;ptot

  Double_t elA  = 0., elM  = 0., elSi = 0., elK = 0., elSk = 0.;
  Double_t piA  = 0., piM  = 0., piSi = 0., piK = 0., piSk = 0.;
  Double_t kaA  = 0., kaM  = 0., kaSi = 0., kaK = 0., kaSk = 0.;
  Double_t prA  = 0., prM  = 0., prSi = 0., prK = 0., prSk = 0.;

  treeLookUp->SetBranchAddress("sign"   ,&sign);
  treeLookUp->SetBranchAddress("it"     ,&it);
  treeLookUp->SetBranchAddress("sl"     ,&sl);
  treeLookUp->SetBranchAddress("myBin"  ,myBin);

  treeLookUp->SetBranchAddress("elM" ,&elM);
  treeLookUp->SetBranchAddress("piM" ,&piM);
  treeLookUp->SetBranchAddress("kaM" ,&kaM);
  treeLookUp->SetBranchAddress("prM" ,&prM);

  treeLookUp->SetBranchAddress("elSi" ,&elSi);
  treeLookUp->SetBranchAddress("piSi" ,&piSi);
  treeLookUp->SetBranchAddress("kaSi" ,&kaSi);
  treeLookUp->SetBranchAddress("prSi" ,&prSi);

  treeLookUp->SetBranchAddress("elA" ,&elA);
  treeLookUp->SetBranchAddress("piA" ,&piA);
  treeLookUp->SetBranchAddress("kaA" ,&kaA);
  treeLookUp->SetBranchAddress("prA" ,&prA);

  treeLookUp->SetBranchAddress("elSk" ,&elSk);
  treeLookUp->SetBranchAddress("piSk" ,&piSk);
  treeLookUp->SetBranchAddress("kaSk" ,&kaSk);
  treeLookUp->SetBranchAddress("prSk" ,&prSk);

  treeLookUp->SetBranchAddress("elK" ,&elK);
  treeLookUp->SetBranchAddress("piK" ,&piK);
  treeLookUp->SetBranchAddress("kaK" ,&kaK);
  treeLookUp->SetBranchAddress("prK" ,&prK);

  for(Int_t i = 0; i < nent; ++i)
  {
    // myBin[0] --> Eta, myBin[1]--> Centrality, myBin[2]-->Momentum
    treeLookUp -> GetEntry(i);
    if (it != fitIter) continue;
    //
    if (sign==1){
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prA;

      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prM;

      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prSi;

      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prSk;

      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kEl] = elK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kPi] = piK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kKa] = kaK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kPr] = prK;
    }

    if (sign==-1){
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = elA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = piA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = kaA;
      fAmpArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = prA;

      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = -elM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = -piM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = -kaM;
      fMeanArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = -prM;

      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = elSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = piSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = kaSi;
      fSigmaArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = prSi;

      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = -elSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = -piSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = -kaSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = -prSk;

      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = elK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = piK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = kaK;
      fKurtosisArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = prK;
    }

  }

  Int_t centBinRange[2] = {TMath::Max(fCentInputBin-1,0), TMath::Min(fCentInputBin+1,fnCentBins)};
  Int_t etaBinRange[2]  = {TMath::Max(fEtaDownBin-1,0), TMath::Min(fEtaUpBin+1,fnEtaBins)};
  Int_t momBinRange[2]  = {TMath::Max(fpDownBin-1,0), TMath::Min(fpUpBin+1,fnMomBins)};
  cout << "==================================" << endl;
  cout << "   Reading the file is started " << endl;
  cout << "==================================" << endl;
  cout << "cent Bin window = " << centBinRange[0] << "  " << centBinRange[1] << endl;
  cout << "eta  Bin window = " << etaBinRange[0]  << "  " << etaBinRange[1] << endl;
  cout << "mom  Bin window = " << momBinRange[0]  << "  " << momBinRange[1] << endl;
  cout << "==================================" << endl;
  //
  // Fill TF1s into TClonesArray
  outFile->cd();
  // for (Int_t i = 0; i<fnEtaBins; i++){
  //   for (Int_t j = 0; j< fnCentBins; j++){
  //     for (Int_t k = 0; k<fnMomBins; k++){
  for (Int_t i = etaBinRange[0]; i<etaBinRange[1]; i++){
    for (Int_t j = centBinRange[0]; j<centBinRange[1]; j++){
      for (Int_t k = momBinRange[0];  k<momBinRange[1]; k++){

        TString objPtotonName     = Form("particle_0_bin_%d_bin_%d_bin_%d",i,j,k);
        TString objAntiProtonName = Form("particle_1_bin_%d_bin_%d_bin_%d",i,j,k);
        TString objOthersName     = Form("particle_2_bin_%d_bin_%d_bin_%d",i,j,k);
        //
        // proton
        fParticles[i][j][k][0] = new TF1(objPtotonName,fitFunctionGenGausStr,fMindEdx,fMaxdEdx); // protons
        fParticles[i][j][k][0]->FixParameter(0,fAmpArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->FixParameter(1,fMeanArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->FixParameter(2,fSigmaArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->FixParameter(3,fKurtosisArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->FixParameter(4,fSkewArr[i][j][k][kPr]);
        fParticles[i][j][k][0]->SetLineColor(kGreen+2);
        fParticles[i][j][k][0]->SetNpx(nBinsLineShape);
        //
        hParticles[i][j][k][0] = (TH1D*)fParticles[i][j][k][0]->GetHistogram();
        if (fSubsample>0) hParticles[i][j][k][0]->Scale(1./nSubSample);
        hParticles[i][j][k][0]->SetName(objPtotonName);
        hParticles[i][j][k][0]->SetLineColor(kGreen+2);
        //
        // antiprotons
        fParticles[i][j][k][1] = new TF1(objAntiProtonName,fitFunctionGenGausStr,fMindEdx,fMaxdEdx); // antiprotons
        fParticles[i][j][k][1]->FixParameter(0,fAmpArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->FixParameter(1,fMeanArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->FixParameter(2,fSigmaArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->FixParameter(3,fKurtosisArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->FixParameter(4,fSkewArr[i][j][k][kBPr]);
        fParticles[i][j][k][1]->SetLineColor(kBlue+2);
        fParticles[i][j][k][1]->SetNpx(nBinsLineShape);
        //
        hParticles[i][j][k][1] = (TH1D*)fParticles[i][j][k][1]->GetHistogram();
        if (fSubsample>0) hParticles[i][j][k][1]->Scale(1./nSubSample);
        hParticles[i][j][k][1]->SetName(objAntiProtonName);
        hParticles[i][j][k][1]->SetLineColor(kBlue+2);

        //
        // Other particles as backgtound
        fParticles[i][j][k][2] = MergeOtherTF1s(i,j,k);
        fParticles[i][j][k][2]->SetName(objOthersName);
        fParticles[i][j][k][2]->SetLineColor(kBlack);
        fParticles[i][j][k][2]->SetNpx(nBinsLineShape);
        //
        hParticles[i][j][k][2] = (TH1D*)fParticles[i][j][k][2]->GetHistogram();
        if (fSubsample>0) hParticles[i][j][k][2]->Scale(1./nSubSample);
        hParticles[i][j][k][2]->SetName(objOthersName);
        hParticles[i][j][k][2]->SetLineColor(kBlack);

        // if (i==etaBinRange[0] && j==centBinRange[0]) {
        //   hParticles[i][j][k][0]->Write();
        //   hParticles[i][j][k][1]->Write();
        //   hParticles[i][j][k][2]->Write();
        // }


      }
    }
  }

}

Bool_t ApplyTreeSelection(UInt_t cut, Int_t syst)
{
  /*
  syst: 0 -->  Reference
  1 -->  CRows60
  2 -->  CRows100
  3 -->  Chi2TPC3
  4 -->  Chi2TPC5
  5 -->  DCAXYSmall
  6 -->  DCAXYLarge
  7 -->  VZSmall
  8 -->  VZLarge
  9 -->  EventVertexZSmall
  10 --> EventVertexZLarge
  11 --> ClusterRequirementITS
  12 --> NewITSCut
  */

  const Int_t fnCutBins=7;
  Int_t fCutArr[fnCutBins]={0};
  if(syst==0) {  // 0 -->  Reference
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==1) {  // 1 -->  CRows60
    Int_t fCutArrTmp[fnCutBins] = {kCRows60,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==2) {  // 2 -->  CRows100
    Int_t fCutArrTmp[fnCutBins] = {kCRows100, kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==3) {  // 3 -->  Chi2TPC3
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC3, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==4) {  // 4 -->  Chi2TPC5
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC5, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==5) {  // 5 -->  kDCAXYSmall
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXYSmall, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==6) {  // 6 -->  kDCAXYLarge
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXYLarge, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==7) {  // 7 -->  kVZSmall
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZSmall, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==8) {  // 8 -->  kVZLarge
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZLarge, kEventVertexZ, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==9) {  // 9 -->  kEventVertexZSmall
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZSmall, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==10) { // 10 -->  kEventVertexZLarge
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZLarge, kClusterRequirementITS, kNewITSCut};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==11) { // 11 -->  kNewITSCut
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kNewITSCut, kCRows80};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  if(syst==12) { // 12 -->  kClusterRequirementITS
    Int_t fCutArrTmp[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kCRows80};
    for(Int_t i=0;i<fnCutBins;i++) fCutArr[i] = fCutArrTmp[i];
  }
  //
  //
  //  Apply conditions
  for (Int_t i=0;i<fnCutBins;i++){
    if( ((cut >> fCutArr[i]) & 1) == 0 ) return kFALSE;
  }

  return kTRUE;

}

// =======================================================================================================
// =======================================================================================================
// =======================================================================================================

void setpars(TF1 *fun, Double_t a, Double_t b, Double_t c)
{
  fun->SetParameter(0, a);
  fun->SetParameter(1, b);
  fun->SetParameter(2, c);
}

Double_t getW2(Int_t mix, Int_t a, Int_t b, Int_t c, Int_t i, Int_t l)
{
  return (wmix[mix][i] * wmean[c][l] - wmean[a][i] * wmean[b][i] * wmean[c][l]);
}

Double_t getValueH(Double_t *xval, Double_t *par)
{
  Double_t xx = xval[0];
  Int_t pp = (Int_t)par[0];
  //
  //
  Double_t piValue=0., kValue=0., prValue=0.;
  if (lookUpTableLineMode==0) {
    Int_t bin = hParticles[fEtaBin][fCentBin][fMomBin][pp]->FindBin(xx);
    prValue = hParticles[fEtaBin][fCentBin][fMomBin][0]->GetBinContent(bin);
    piValue = hParticles[fEtaBin][fCentBin][fMomBin][1]->GetBinContent(bin);
    kValue  = hParticles[fEtaBin][fCentBin][fMomBin][2]->GetBinContent(bin);
  }
  if (lookUpTableLineMode==1) {
    prValue = fParticles[fEtaBin][fCentBin][fMomBin][0]->Eval(xx);
    piValue = fParticles[fEtaBin][fCentBin][fMomBin][1]->Eval(xx);
    kValue  = fParticles[fEtaBin][fCentBin][fMomBin][2]->Eval(xx);
  }

  if (pp == 0) return prValue;
  else if (pp == 1) return piValue;
  else if (pp == 2) return kValue;
  else return (piValue + kValue + prValue);

}

Double_t myIntegral(TF1 *fun)
{
  Double_t xx, sum = 0;
  Double_t step = (fMaxdEdx - fMindEdx) / Double_t(2 * size_size);
  for (Int_t i = 1; i < 2 * size_size + 1; i++)
  {
    xx = fMindEdx + step * i;
    sum += fun->Eval(xx);
  }
  return sum * step;
}

void initFunctions_M()
{

  if (funProton) delete funProton;
  funProton = new TF1("funProton", getValueH, fMindEdx, fMaxdEdx, 1);
  funProton->SetParameter(0, 0);

  if (funAProton)  delete funAProton;
  funAProton = new TF1("funAProton", getValueH, fMindEdx, fMaxdEdx, 1);
  funAProton->SetParameter(0, 1);

  if (funBackground)  delete funBackground;
  funBackground = new TF1("funBackground", getValueH, fMindEdx, fMaxdEdx, 1);
  funBackground->SetParameter(0, 2);

}

Double_t getFunctions(Double_t *xx, Double_t *par)
{
  Double_t x = xx[0];
  Int_t j = (Int_t)par[0];
  Int_t i = (Int_t)par[1];
  Int_t k = (Int_t)par[2];
  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funBackground->Eval(x);
  val[2] = funAProton->Eval(x);

  Double_t sumVal = val[0] + val[1] + val[2];
  Double_t relVal[3];
  if (sumVal < 1e-15) return 0.;
  Double_t power = k;
  for (Int_t m = 0; m < 3; m++)
  {
    relVal[m] = val[m] / sumVal;
    relVal[m] = TMath::Power(relVal[m], power);
  }
  return relVal[j] * val[i];
}

Double_t getFunctionsMix(Double_t *xx, Double_t *par)
{
  Double_t x = xx[0];
  Int_t j = (Int_t)par[0];
  Int_t k = (Int_t)par[1];

  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funBackground->Eval(x);
  val[2] = funAProton->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15) return 0.;

  Double_t myVal[3];
  myVal[0] = val[0] * val[1] / sumVal / sumVal;
  myVal[1] = val[0] * val[2] / sumVal / sumVal;
  myVal[2] = val[1] * val[2] / sumVal / sumVal;

  return myVal[j] * val[k];
}

Double_t getFunctionsMix2(Double_t *xx, Double_t *par)
{

  Double_t x = xx[0];
  Int_t m = (Int_t)par[0];
  Int_t j = (Int_t)par[1];
  Int_t k = (Int_t)par[2];

  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funBackground->Eval(x);
  val[2] = funAProton->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15) return 0.;

  Double_t myVal[3];
  myVal[0] = val[m] * val[m] * val[0] / sumVal / sumVal / sumVal;
  myVal[1] = val[m] * val[m] * val[1] / sumVal / sumVal / sumVal;
  myVal[2] = val[m] * val[m] * val[2] / sumVal / sumVal / sumVal;

  return myVal[j] * val[k];
}

Double_t getFunctionsMix3(Double_t *xx, Double_t *par)
{
  Double_t x = xx[0];
  Int_t k = (Int_t)par[0];
  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funBackground->Eval(x);
  val[2] = funAProton->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15) return 0.;

  Double_t myVal;
  myVal = val[0] * val[1] * val[2] / sumVal / sumVal / sumVal;
  return myVal * val[k];
}

void myDelete(TF1 *fun)
{
  if (fun)
  {
    delete fun;
    fun = NULL;
  }
}

void initFunctions()
{

  funP_P = new TF1("funP_P", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_P, 0, 0, 1);
  funK_P = new TF1("funK_P", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_P, 1, 0, 1);
  funAP_P = new TF1("funAP_P", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_P, 2, 0, 1);
  funP_P2 = new TF1("funP_P2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_P2, 0, 0, 2);
  funK_P2 = new TF1("funK_P2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_P2, 1, 0, 2);
  funAP_P2 = new TF1("funAP_P2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_P2, 2, 0, 2);
  funP_P3 = new TF1("funP_P3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_P3, 0, 0, 3);
  funK_P3 = new TF1("funK_P3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_P3, 1, 0, 3);
  funAP_P3 = new TF1("funAP_P3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_P3, 2, 0, 3);
  funP_K = new TF1("funP_K", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_K, 0, 1, 1);
  funK_K = new TF1("funK_K", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_K, 1, 1, 1);
  funAP_K = new TF1("funAP_K", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_K, 2, 1, 1);
  funP_K2 = new TF1("funP_K2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_K2, 0, 1, 2);
  funK_K2 = new TF1("funK_K2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_K2, 1, 1, 2);
  funAP_K2 = new TF1("funAP_K2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_K2, 2, 1, 2);
  funP_K3 = new TF1("funP_K3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_K3, 0, 1, 3);
  funK_K3 = new TF1("funK_K3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_K3, 1, 1, 3);
  funAP_K3 = new TF1("funAP_K3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_K3, 2, 1, 3);
  funP_AP = new TF1("funP_AP", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_AP, 0, 2, 1);
  funK_AP = new TF1("funK_AP", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_AP, 1, 2, 1);
  funAP_AP = new TF1("funAP_AP", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_AP, 2, 2, 1);
  funP_AP2 = new TF1("funP_AP2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_AP2, 0, 2, 2);
  funK_AP2 = new TF1("funK_AP2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_AP2, 1, 2, 2);
  funAP_AP2 = new TF1("funAP_AP2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_AP2, 2, 2, 2);
  funP_AP3 = new TF1("funP_AP3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_AP3, 0, 2, 3);
  funK_AP3 = new TF1("funK_AP3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_AP3, 1, 2, 3);
  funAP_AP3 = new TF1("funAP_AP3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funAP_AP3, 2, 2, 3);
  ////////////////// mixed ones
  funPK_P = new TF1("funPK_P", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPK_P, 0, 0, -10);
  funPAP_P = new TF1("funPAP_P", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPAP_P, 1, 0, -10);
  funAPK_P = new TF1("funAPK_P", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funAPK_P, 2, 0, -10);
  funPK_K = new TF1("funPK_K", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPK_K, 0, 1, -10);
  funPAP_K = new TF1("funPAP_K", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPAP_K, 1, 1, -10);
  funAPK_K = new TF1("funAPK_K", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funAPK_K, 2, 1, -10);
  funPK_AP = new TF1("funPK_AP", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPK_AP, 0, 2, -10);
  funPAP_AP = new TF1("funPAP_AP", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPAP_AP, 1, 2, -10);
  funAPK_AP = new TF1("funAPK_AP", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funAPK_AP, 2, 2, -10);

  funP2AP_P = new TF1("funP2AP_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2AP_P, 0, 2, 0);
  funP2AP_AP = new TF1("funP2AP_AP", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2AP_AP, 0, 2, 2);
  funP2AP_K = new TF1("funP2AP_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2AP_K, 0, 2, 1);
  funAP2P_P = new TF1("funAP2P_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funAP2P_P, 2, 0, 0);
  funAP2P_AP = new TF1("funAP2P_AP", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funAP2P_AP, 2, 0, 2);
  funAP2P_K = new TF1("funAP2P_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funAP2P_K, 2, 0, 1);
  funP2K_P = new TF1("funP2K_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2K_P, 0, 1, 0);
  funP2K_AP = new TF1("funP2K_AP", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2K_AP, 0, 1, 2);
  funP2K_K = new TF1("funP2K_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2K_K, 0, 1, 1);
  funK2P_P = new TF1("funK2P_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2P_P, 1, 0, 0);
  funK2P_AP = new TF1("funK2P_AP", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2P_AP, 1, 0, 2);
  funK2P_K = new TF1("funK2P_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2P_K, 1, 0, 1);
  funAP2K_P = new TF1("funAP2K_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funAP2K_P, 2, 1, 0);
  funAP2K_AP = new TF1("funAP2K_AP", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funAP2K_AP, 2, 1, 2);
  funAP2K_K = new TF1("funAP2K_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funAP2K_K, 2, 1, 1);
  funK2AP_P = new TF1("funK2AP_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2AP_P, 1, 2, 0);
  funK2AP_AP = new TF1("funK2AP_AP", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2AP_AP, 1, 2, 2);
  funK2AP_K = new TF1("funK2AP_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2AP_K, 1, 2, 1);
  funAPPrK_P = new TF1("funAPPrK_P", getFunctionsMix3, fMindEdx, fMaxdEdx, 3);
  setpars(funAPPrK_P, 0, -100, -100);
  funAPPrK_K = new TF1("funAPPrK_K", getFunctionsMix3, fMindEdx, fMaxdEdx, 3);
  setpars(funAPPrK_K, 1, -100, -100);
  funAPPrK_AP = new TF1("funAPPrK_AP", getFunctionsMix3, fMindEdx, fMaxdEdx, 3);
  setpars(funAPPrK_AP, 2, -100, -100);
}

void calcIntegrals(Float_t *normme)
{
  Double_t npNorm = normme[0];  // * 0.01;
  Double_t npiNorm = normme[1]; /// * 0.01;
  Double_t nkNorm = normme[2];  /// * 0.01;

  if (nkNorm == 0)
    nkNorm = 1;
  if (npNorm == 0)
    npNorm = 1;
  if (npiNorm == 0)
    npiNorm = 1;

  wP_P += myIntegral(funP_P) / npNorm;
  wK_P += myIntegral(funK_P) / npNorm;
  wAP_P += myIntegral(funAP_P) / npNorm;

  wP_P2 += myIntegral(funP_P2) / npNorm;
  wK_P2 += myIntegral(funK_P2) / npNorm;
  wAP_P2 += myIntegral(funAP_P2) / npNorm;

  wP_P3 += myIntegral(funP_P3) / npNorm;
  wK_P3 += myIntegral(funK_P3) / npNorm;
  wAP_P3 += myIntegral(funAP_P3) / npNorm;

  wP_K += myIntegral(funP_K) / nkNorm;
  wK_K += myIntegral(funK_K) / nkNorm;
  wAP_K += myIntegral(funAP_K) / nkNorm;

  wP_K2 += myIntegral(funP_K2) / nkNorm;
  wK_K2 += myIntegral(funK_K2) / nkNorm;
  wAP_K2 += myIntegral(funAP_K2) / nkNorm;

  wP_K3 += myIntegral(funP_K3) / nkNorm;
  wK_K3 += myIntegral(funK_K3) / nkNorm;
  wAP_K3 += myIntegral(funAP_K3) / nkNorm;

  wP_AP += myIntegral(funP_AP) / npiNorm;
  wK_AP += myIntegral(funK_AP) / npiNorm;
  wAP_AP += myIntegral(funAP_AP) / npiNorm;

  wP_AP2 += myIntegral(funP_AP2) / npiNorm;
  wK_AP2 += myIntegral(funK_AP2) / npiNorm;
  wAP_AP2 += myIntegral(funAP_AP2) / npiNorm;

  wP_AP3 += myIntegral(funP_AP3) / npiNorm;
  wK_AP3 += myIntegral(funK_AP3) / npiNorm;
  wAP_AP3 += myIntegral(funAP_AP3) / npiNorm;

  ///mixed
  wPK_P += myIntegral(funPK_P) / npNorm;
  wPAP_P += myIntegral(funPAP_P) / npNorm;
  wAPK_P += myIntegral(funAPK_P) / npNorm;

  wPK_K += myIntegral(funPK_K) / nkNorm;
  wPAP_K += myIntegral(funPAP_K) / nkNorm;
  wAPK_K += myIntegral(funAPK_K) / nkNorm;

  wPK_AP += myIntegral(funPK_AP) / npiNorm;
  wPAP_AP += myIntegral(funPAP_AP) / npiNorm;
  wAPK_AP += myIntegral(funAPK_AP) / npiNorm;

  ///////////////////////////////////////////////////////////////////////

  wP2AP_P += myIntegral(funP2AP_P) / npNorm;
  wP2AP_AP += myIntegral(funP2AP_AP) / npiNorm;
  wP2AP_K += myIntegral(funP2AP_K) / nkNorm;

  wAP2P_P += myIntegral(funAP2P_P) / npNorm;
  wAP2P_AP += myIntegral(funAP2P_AP) / npiNorm;
  wAP2P_K += myIntegral(funAP2P_K) / nkNorm;

  wP2K_P += myIntegral(funP2K_P) / npNorm;
  wP2K_AP += myIntegral(funP2K_AP) / npiNorm;
  wP2K_K += myIntegral(funP2K_K) / nkNorm;

  wK2P_P += myIntegral(funK2P_P) / npNorm;
  wK2P_AP += myIntegral(funK2P_AP) / npiNorm;
  wK2P_K += myIntegral(funK2P_K) / nkNorm;

  wAP2K_P += myIntegral(funAP2K_P) / npNorm;
  wAP2K_AP += myIntegral(funAP2K_AP) / npiNorm;
  wAP2K_K += myIntegral(funAP2K_K) / nkNorm;

  wK2AP_P += myIntegral(funK2AP_P) / npNorm;
  wK2AP_AP += myIntegral(funK2AP_AP) / npiNorm;
  wK2AP_K += myIntegral(funK2AP_K) / nkNorm;

  wAPPrK_Pr += myIntegral(funAPPrK_P) / npNorm;
  wAPPrK_K += myIntegral(funAPPrK_K) / nkNorm;
  wAPPrK_AP += myIntegral(funAPPrK_AP) / npiNorm;
}

/////////////////////////////

template <class T>
Float_t Sum(vector<T> &v)
{
  return accumulate(v.begin(), v.end(), 0.0);
}

void resetValues()
{
  WP_sum = WK_sum = WAP_sum = WMult_sum = 0.;
}

void addParticles(Float_t mVal, Int_t &count)
{
  Double_t kValue, prValue, piValue;
  Double_t WP1, WK1, WAP1;
  //
  prValue = funProton->Eval(mVal);
  kValue  = funBackground->Eval(mVal);
  piValue = funAProton->Eval(mVal);

  Double_t sumValue = prValue + kValue + piValue;
  if (sumValue == 0)
  {
    cout << "why:   " << mVal << endl;
    wrongCount++;
    WP1 = WK1 = WAP1 = 0;
  }
  else
  {
    WP1 = prValue / sumValue;
    WK1 = kValue / sumValue;
    WAP1 = piValue / sumValue;

    count++;

    WP_sum += WP1;
    WK_sum += WK1;
    WAP_sum += WAP1;
    WMult_sum += 1;
  }
}

int identity3Particle()
{

  cout << " identity3Particle.Info: burada " << endl;
  TFile *dataFile = new TFile(fileNameDataTree);
  wrongCount = 0;
  TTree *myTree = (TTree *)dataFile->Get(treeIdentity);
  cout << " identity3Particle.Info: Funtions are inited" << endl;
  Int_t prevEvt = -1000;
  cout << " identity3Particle.Info: set address" << endl;
  //
  Int_t fEventNum;  // ULong64_t fEventNum;
  Double_t fDEdx;   // Float_t fDEdx;
  UInt_t fCutBit;
  Int_t fSign;
  TBranch *fMyBinBrach = (TBranch *)myTree->FindBranch("myBin");
  if (!fMyBinBrach)
  { // new version of tree format
    myTree->SetBranchAddress("gid",    &fEventNum);
    myTree->SetBranchAddress("dEdx",   &fDEdx);
    myTree->SetBranchAddress("sign",   &fSign);
    myTree->SetBranchAddress("cutBit", &fCutBit);
  }
  else
  { // old version of tree format
    myTree->SetBranchAddress("evtNum", &fEventNum);
    myTree->SetBranchAddress("myDeDx", &fDEdx);
    myTree->SetBranchAddress("sign",   &fSign);
    myTree->SetBranchAddress("myBin",  fMyBin);
  }
  cout << "setted address" << endl;
  resetValues();
  //
  Int_t count = 0;
  for (Int_t i = 0; i < fnEtaBins; i++){
    for (Int_t j = 0; j < fnCentBins; j++){
      for (Int_t k = 0; k < fnMomBins; k++){
        meanBackground[i][j][k] = 0;
        meanProton[i][j][k] = 0;
        meanAntiProton[i][j][k] = 0;
        fUsedBins[i][j][k] = 0;
      }
    }
  }

  numAllEvents = 0;
  numAllCutEvents = 0;

  Int_t runEvt = (Int_t)myTree->GetEntries();
  cout << " identity3Particle.Info: total number " << runEvt << endl;
  Int_t runStat1 = 0;
  Int_t runStat2 = runEvt;

  cout << " identity3Particle.Info: running from " << runStat1 << " to " << runStat2 << endl;
  Int_t countVeto = 0;
  Int_t prevEvtVeto = -1;
  initFunctions_M();

  myTree->GetEntry(runStat2 - 1);
  Int_t remEvent = fEventNum;
  cout << " ==================================" << endl;
  cout << " identity3Particle.Info: loop over tracks starts " <<endl;
  timer.Reset(); timer.Start();
  cout << " ==================================" << endl;
  for (Int_t i = runStat1; i < runStat2; i++)
  {
    if (i % 2000000 == 0) cout << "track " << i << " of " << runStat2 - runStat1 << endl;
    if (i>fTestNtracks && fTestNtracks>0) break;
    myTree->GetEntry(i);
    if ((Int_t)fEventNum == remEvent) continue;
    if ((Int_t)fEventNum != prevEvtVeto && prevEvtVeto > 0) countVeto++;
    fEtaBin  = fMyBin[0];
    fCentBin = fMyBin[1];
    fMomBin  = fMyBin[2];
    prevEvtVeto = fEventNum;
    fDEdx = fSign*fDEdx;

    if( fDEdx < fMindEdx || fDEdx > fMaxdEdx) continue;
    if( fCentBin != fCentInputBin ) continue;
    if( fMomBin < fpDownBin   || fMomBin > fpUpBin ) continue;
    if( fEtaBin < fEtaDownBin || fEtaBin > fEtaUpBin ) continue;
    //
    // Fill debug historam
    if (fEtaBin==fEtaDownBin && fCentBin==fCentInputBin && fMomBin < fpUpBin){
      hDedxDebug[fMomBin-fpDownBin]->Fill(fDEdx);
    }
    //
    fUsedBins[fEtaBin][fCentBin][fMomBin] = 1;
    // cout << " aaaa " << hParticles[fEtaBin][fCentBin][fMomBin][0]->GetMean() << "   " <<  fEtaBin << "  " << fCentBin << "  " << fMomBin << endl;
    if (!meanProton[fEtaBin][fCentBin][fMomBin]) meanProton[fEtaBin][fCentBin][fMomBin] = myIntegral(funProton);
    if (!meanAntiProton[fEtaBin][fCentBin][fMomBin])   meanAntiProton[fEtaBin][fCentBin][fMomBin]   = myIntegral(funAProton);
    if (!meanBackground[fEtaBin][fCentBin][fMomBin])   meanBackground[fEtaBin][fCentBin][fMomBin]   = myIntegral(funBackground);

    numAllEvents += multEv;

    if ((Int_t)fEventNum == prevEvt)
    {
      addParticles(fDEdx, count);
    }
    else
    {
      if (count != 0)
      {
        if (WMult_sum)
        {
          sumMult.push_back(count);

          W3P_sum_vec.push_back(TMath::Power(WP_sum, 3));
          W3AP_sum_vec.push_back(TMath::Power(WAP_sum, 3));
          W3K_sum_vec.push_back(TMath::Power(WK_sum, 3));

          W2P_sum_vec.push_back(WP_sum * WP_sum);
          W2K_sum_vec.push_back(WK_sum * WK_sum);
          W2AP_sum_vec.push_back(WAP_sum * WAP_sum);

          WPM_sum_vec.push_back(WP_sum);
          WKM_sum_vec.push_back(WK_sum);
          WAPM_sum_vec.push_back(WAP_sum);

          WPK_sum_vec.push_back(WP_sum * WK_sum);
          WPAP_sum_vec.push_back(WP_sum * WAP_sum);
          WAPK_sum_vec.push_back(WAP_sum * WK_sum);

          WAPPrK_sum_vec.push_back(WAP_sum * WP_sum * WK_sum);

          WPr2AP_sum_vec.push_back(TMath::Power(WP_sum, 2) * WAP_sum);
          WPr2K_sum_vec.push_back(TMath::Power(WP_sum, 2) * WK_sum);
          WAP2K_sum_vec.push_back(TMath::Power(WAP_sum, 2) * WK_sum);

          WAP2Pr_sum_vec.push_back(WAP_sum * WAP_sum * WP_sum);
          WK2Pr_sum_vec.push_back(WK_sum * WK_sum * WP_sum);
          WK2AP_sum_vec.push_back(WK_sum * WK_sum * WAP_sum);
        }
      }
      resetValues();
      count = 0;
      addParticles(fDEdx, count);
      numAllCutEvents++;
    }
    prevEvt = fEventNum;
  } //end event
  //
  // dump some debug histograms
  outFile->cd();
  if (fSubsample<2){
    for (Int_t i=0; i<fpUpBin-fpDownBin-1; i++) {
      hDedxDebug[i]->Write();
      hParticles[fEtaDownBin][fCentInputBin][fpDownBin+i][0]->Write();
      hParticles[fEtaDownBin][fCentInputBin][fpDownBin+i][1]->Write();
      hParticles[fEtaDownBin][fCentInputBin][fpDownBin+i][2]->Write();
    }
  }


  cout << "====================================" << endl;
  cout << " track loop is over " << endl;
  timer.Stop(); timer.Print();
  cout << "====================================" << endl;
  //
  nEvents = countVeto;
  nEvents = sumMult.size();
  countVeto = sumMult.size();

  cout << "# of Events == " << nEvents << endl;
  Double_t corrFactor = multEv / nEvents;
  multEv = countVeto;
  corrFactor = 1;
  cout << "corrFactor == " << corrFactor << endl;

  ////////////////////

  Double_t meanMult = Sum<int>(sumMult) / nEvents;
  cout << "mean mult   " << meanMult << endl;

  Double_t proton_aver = 0;
  Double_t background_aver = 0;
  Double_t antiP_aver = 0;

  for (Int_t i = 0; i < fnEtaBins; i++){
    for (Int_t j = 0; j < fnCentBins; j++){
      for (Int_t k = 0; k < fnMomBins; k++){
        if(fUsedBins[i][j][k] != 1) continue;
        background_aver += meanBackground[i][j][k];
        proton_aver += meanProton[i][j][k];
        antiP_aver += meanAntiProton[i][j][k];
      }
    }
  }

  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "mean multiplicities : ONLY TRUE FOR ALL STATISTICS !" << endl;
  cout << "proton       " << proton_aver * corrFactor << endl;
  cout << "background   " << background_aver * corrFactor << endl;
  cout << "antiP        " << antiP_aver * corrFactor << endl;

  proton_aver *= corrFactor;
  background_aver *= corrFactor;
  antiP_aver *= corrFactor;

  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "***********" << endl;

  //Double_t averMult = Sum<double>(WMult_sum_vec)/nEvents;
  //Double_t averMult2 = Sum<double>(W2Mult_sum_vec)/nEvents;

  Double_t W2P_aver = Sum<double>(W2P_sum_vec) / nEvents;
  Double_t W2K_aver = Sum<double>(W2K_sum_vec) / nEvents;
  Double_t W2AP_aver = Sum<double>(W2AP_sum_vec) / nEvents;

  Double_t WPM_aver = Sum<double>(WPM_sum_vec) / nEvents;
  Double_t WKM_aver = Sum<double>(WKM_sum_vec) / nEvents;
  Double_t WAPM_aver = Sum<double>(WAPM_sum_vec) / nEvents;

  Double_t W3P_aver = Sum<double>(W3P_sum_vec) / nEvents;
  Double_t W3AP_aver = Sum<double>(W3AP_sum_vec) / nEvents;
  Double_t W3K_aver = Sum<double>(W3K_sum_vec) / nEvents;
  Double_t WAPPrK_aver = Sum<double>(WAPPrK_sum_vec) / nEvents;
  //    Double_t WAP2P_aver     = Sum<double>(WAP2P_sum_vec)/nEvents;

  Double_t WPr2AP_aver = Sum<double>(WPr2AP_sum_vec) / nEvents;
  Double_t WPr2K_aver = Sum<double>(WPr2K_sum_vec) / nEvents;
  Double_t WAP2K_aver = Sum<double>(WAP2K_sum_vec) / nEvents;

  Double_t WAP2Pr_aver = Sum<double>(WAP2Pr_sum_vec) / nEvents;
  Double_t WK2Pr_aver = Sum<double>(WK2Pr_sum_vec) / nEvents;
  Double_t WK2AP_aver = Sum<double>(WK2AP_sum_vec) / nEvents;

  Double_t WPK_aver = Sum<double>(WPK_sum_vec) / nEvents;
  Double_t WPAP_aver = Sum<double>(WPAP_sum_vec) / nEvents;
  Double_t WAPK_aver = Sum<double>(WAPK_sum_vec) / nEvents;

  Double_t pr_aver = 0.;  //Sum<int>(P_mult_vec)/nEvents;
  Double_t k_aver = 0.;   ///Sum<int>(K_mult_vec)/nEvents;
  Double_t apr_aver = 0.;  ///Sum<int>(AP_mult_vec)/nEvents;

  Double_t pr2_aver = 0.; ////Sum<int>(P2_mult_vec)/nEvents;
  Double_t k2_aver = 0.;  ////Sum<int>(K2_mult_vec)/nEvents;
  Double_t pi2_aver = 0.; ///Sum<int>(AP2_mult_vec)/nEvents;

  Double_t prapr_aver = 0.; ///Sum<int>(PAP_mult_vec)/nEvents;
  Double_t prk_aver = 0.;  ///Sum<int>(PK_mult_vec)/nEvents;
  Double_t pik_aver = 0.;  ////Sum<int>(APK_mult_vec)/nEvents;

  cout << " =========================== " << endl;
  cout << " W2P_aver " << W2P_aver << endl;
  cout << " W2K_aver " << W2K_aver << endl;
  cout << " W2AP_aver " << W2AP_aver << endl;
  cout << " WPK_aver " << WPK_aver << endl;
  cout << " WPAP_aver " << WPAP_aver << endl;
  cout << " W2APK_aver " << WAPK_aver << endl;
  cout << " =========================== " << endl;
  cout << " mean from identities " << endl;
  cout << " =========================== " << endl;
  cout << " WPM_aver  " << WPM_aver  << " -> " << proton_aver << endl;
  cout << " WKM_aver  " << WKM_aver  << " -> " << background_aver   << endl;
  cout << " WAPM_aver " << WAPM_aver << " -> " << antiP_aver   <<  endl;
  cout << " =========================== " << endl;
  (*fMoments1st)[kPR]=WPM_aver;
  (*fMoments1st)[kAP]=WAPM_aver;
  (*fMoments1st)[kBG]=WKM_aver;
  /*
  proton_aver   = WPM_aver;
  background_aver     = WKM_aver;
  antiP_aver     = WAPM_aver;
  */
  //
  wP_P = 0., wK_P = 0., wAP_P = 0., wP_P2 = 0., wK_P2 = 0., wAP_P2 = 0., wP_P3 = 0., wK_P3 = 0., wAP_P3 = 0., wP_K = 0., wK_K = 0., wAP_K = 0., wP_K2 = 0., wK_K2 = 0.;
  wAP_K2 = 0., wP_K3 = 0., wK_K3 = 0., wAP_K3 = 0., wP_AP = 0., wK_AP = 0., wAP_AP = 0., wP_AP2 = 0., wK_AP2 = 0., wAP_AP2 = 0., wP_AP3 = 0., wK_AP3 = 0., wAP_AP3 = 0.;
  ///mixed
  wPK_P = 0., wPAP_P = 0., wAPK_P = 0., wPK_K = 0., wPAP_K = 0., wAPK_K = 0., wPK_AP = 0., wPAP_AP = 0., wAPK_AP = 0., wP2AP_P = 0., wP2AP_AP = 0., wP2AP_K = 0.;
  wAP2P_P = 0., wAP2P_AP = 0., wAP2P_K = 0., wP2K_P = 0., wP2K_AP = 0., wP2K_K = 0., wK2P_P = 0., wK2P_AP = 0., wK2P_K = 0., wAP2K_P = 0., wAP2K_AP = 0., wAP2K_K = 0.;
  wK2AP_P = 0., wK2AP_AP = 0., wK2AP_K = 0., wAPPrK_Pr = 0., wAPPrK_K = 0., wAPPrK_AP = 0.;
  Float_t normme[] = {(Float_t)proton_aver, (Float_t)antiP_aver, (Float_t)background_aver};
  ///////////
  cout << " ==================================" << endl;
  cout << " main calculation starts " <<endl;
  timer.Reset(); timer.Start();
  cout << " ==================================" << endl;
  cout << " initFunctions  " << endl;
  //
  //
  initFunctions();
  cout << " calcIntegrals  " << endl;
  cout << " ==================================" << endl;
  for (Int_t i = 0; i < fnEtaBins; i++){
    for (Int_t j = 0; j < fnCentBins; j++){
      for (Int_t k = 0; k < fnMomBins; k++){
        if(fUsedBins[i][j][k] != 1) continue;
        fEtaBin   =  i;
        fCentBin  =  j;
        fMomBin   =  k;
        calcIntegrals(normme);
      }
    }
  }
  pr_aver = proton_aver;
  apr_aver = antiP_aver;
  k_aver = background_aver;

  cout << " ==================================" << endl;
  cout << " test " << wK_P << "  " << wK_K << "  " << wK_AP << endl;
  cout << " 1 " << W2K_aver << endl;
  cout << " 2 " << pr_aver * (wK_P2 - wK_P * wK_P) << endl;
  cout << " 3 " << apr_aver * (wK_AP2 - wK_AP * wK_AP) << endl;
  cout << " 4 " << k_aver * (wK_K2 - wK_K * wK_K) << endl;
  cout << " ==================================" << endl;
  cout << " ==================================" << endl;
  cout << "  Unfolding      " << endl;
  cout << " ==================================" << endl;
  cout << " ==================================" << endl;
  TMatrixD A2(6, 6);
  A2(0, 0) = wP_P * wP_P;
  A2(0, 1) = wP_AP * wP_AP;
  A2(0, 2) = wP_K * wP_K;
  A2(0, 3) = 2. * wP_P * wP_AP;
  A2(0, 4) = 2. * wP_P * wP_K;
  A2(0, 5) = 2. * wP_AP * wP_K;

  A2(1, 0) = wAP_P * wAP_P;
  A2(1, 1) = wAP_AP * wAP_AP;
  A2(1, 2) = wAP_K * wAP_K;
  A2(1, 3) = 2. * wAP_P * wAP_AP;
  A2(1, 4) = 2. * wAP_P * wAP_K;
  A2(1, 5) = 2. * wAP_AP * wAP_K;

  A2(2, 0) = wK_P * wK_P;
  A2(2, 1) = wK_AP * wK_AP;
  A2(2, 2) = wK_K * wK_K;
  A2(2, 3) = 2. * wK_P * wK_AP;
  A2(2, 4) = 2. * wK_P * wK_K;
  A2(2, 5) = 2. * wK_AP * wK_K;

  A2(3, 0) = wP_P * wAP_P;
  A2(3, 1) = wP_AP * wAP_AP;
  A2(3, 2) = wP_K * wAP_K;
  A2(3, 3) = wP_P * wAP_AP + wP_AP * wAP_P;
  A2(3, 4) = wP_P * wAP_K + wP_K * wAP_P;
  A2(3, 5) = wP_AP * wAP_K + wP_K * wAP_AP;

  A2(4, 0) = wP_P * wK_P;
  A2(4, 1) = wP_AP * wK_AP;
  A2(4, 2) = wP_K * wK_K;
  A2(4, 3) = wP_P * wK_AP + wP_AP * wK_P;
  A2(4, 4) = wP_P * wK_K + wP_K * wK_P;
  A2(4, 5) = wP_AP * wK_K + wP_K * wK_AP;

  A2(5, 0) = wAP_P * wK_P;
  A2(5, 1) = wAP_AP * wK_AP;
  A2(5, 2) = wAP_K * wK_K;
  A2(5, 3) = wAP_P * wK_AP + wAP_AP * wK_P;
  A2(5, 4) = wAP_P * wK_K + wAP_K * wK_P;
  A2(5, 5) = wAP_AP * wK_K + wAP_K * wK_AP;

  proton_aver = WPM_aver;
  background_aver = WKM_aver;
  antiP_aver = WAPM_aver;

  pr_aver = proton_aver;
  apr_aver = antiP_aver;
  k_aver = background_aver;

  Double_t B2[6];
  B2[0] = W2P_aver -
          pr_aver * (wP_P2 - wP_P * wP_P) -
          apr_aver * (wP_AP2 - wP_AP * wP_AP) -
          k_aver * (wP_K2 - wP_K * wP_K);

  B2[1] = W2AP_aver -
          pr_aver * (wAP_P2 - wAP_P * wAP_P) -
          apr_aver * (wAP_AP2 - wAP_AP * wAP_AP) -
          k_aver * (wAP_K2 - wAP_K * wAP_K);

  B2[2] = W2K_aver -
          pr_aver * (wK_P2 - wK_P * wK_P) -
          apr_aver * (wK_AP2 - wK_AP * wK_AP) -
          k_aver * (wK_K2 - wK_K * wK_K);

  B2[3] = WPAP_aver -
          pr_aver * (wPAP_P - wP_P * wAP_P) -
          apr_aver * (wPAP_AP - wP_AP * wAP_AP) -
          k_aver * (wPAP_K - wP_K * wAP_K);

  B2[4] = WPK_aver -
          pr_aver * (wPK_P - wP_P * wK_P) -
          apr_aver * (wPK_AP - wP_AP * wK_AP) -
          k_aver * (wPK_K - wP_K * wK_K);

  B2[5] = WAPK_aver -
          pr_aver * (wAPK_P - wAP_P * wK_P) -
          apr_aver * (wAPK_AP - wAP_AP * wK_AP) -
          k_aver * (wAPK_K - wAP_K * wK_K);

  cout << "A2(0,0) " << A2(0, 0) << "  " << B2[0] << endl;
  A2.Print();

  TMatrixD invA2 = A2.Invert();
  //   invA.Print();
  Double_t recP2_av = 0.;
  Double_t recAP2_av = 0.;
  Double_t recK2_av = 0.;
  Double_t recPAP_av = 0.;
  Double_t recPK_av = 0.;
  Double_t recAPK_av = 0.;

  //recP2_av = recAP2_av = recK2_av  = 0.;
  //recPAP_av = recPK_av = recAPK_av = 0.;
  for (Int_t tt = 0; tt < 6; tt++)
  {
    recP2_av  += invA2(0, tt) * B2[tt];
    recAP2_av += invA2(1, tt) * B2[tt];
    recK2_av  += invA2(2, tt) * B2[tt];
    //cout<< "background test "<<invA2(2,tt)<<"  "<<B2[tt]<<endl;
    recPAP_av += invA2(3, tt) * B2[tt];
    recPK_av  += invA2(4, tt) * B2[tt];
    recAPK_av += invA2(5, tt) * B2[tt];
  }

  cout << " ==================================" << endl;
  cout << " ==================================" << endl;
  cout << " print first moments " << endl;
  cout << " ==================================" << endl;
  cout << " ==================================" << endl;
  cout << " pr  " << WPM_aver   << endl;
  cout << " apr " << WAPM_aver  << endl;
  cout << " bg  " << WKM_aver   << endl;
  cout << " ==================================" << endl;
  cout << " ==================================" << endl;
  cout << " print second moments " << endl;
  cout << " ==================================" << endl;
  cout << " ==================================" << endl;
  cout << " pr2     " << recP2_av  << endl;
  cout << " apr2    " << recAP2_av << endl;
  cout << " bg2     " << recK2_av  << endl;
  cout << " pr_apr  " << recPAP_av << endl;
  cout << " pr_bg   " << recPK_av  << endl;
  cout << " apr_bg  " << recAPK_av << endl;
  cout << " ==================================" << endl;
  cout << " ==================================" << endl;
  (*fMoments2nd)[0]=recP2_av;
  (*fMoments2nd)[1]=recAP2_av;
  (*fMoments2nd)[2]=recK2_av;
  (*fMoments2nd)[3]=recPAP_av;
  (*fMoments2nd)[4]=recPK_av;
  (*fMoments2nd)[5]=recAPK_av;
  //
  pr2_aver = recP2_av;
  pi2_aver = recAP2_av;
  k2_aver = recK2_av;
  prapr_aver = recPAP_av;
  prk_aver = recPK_av;
  pik_aver = recAPK_av;

  TMatrixD A(10, 10);
  TMatrixD BB(10, 9);

  wmean[0][0] = wP_P;
  wmean[0][1] = wP_AP;
  wmean[0][2] = wP_K;
  wmean[1][0] = wAP_P;
  wmean[1][1] = wAP_AP;
  wmean[1][2] = wAP_K;
  wmean[2][0] = wK_P;
  wmean[2][1] = wK_AP;
  wmean[2][2] = wK_K;

  wmean2[0][0] = wP_P2;
  wmean2[0][1] = wP_AP2;
  wmean2[0][2] = wP_K2;

  wmean2[1][0] = wAP_P2;
  wmean2[1][1] = wAP_AP2;
  wmean2[1][2] = wAP_K2;

  wmean2[2][0] = wK_P2;
  wmean2[2][1] = wK_AP2;
  wmean2[2][2] = wK_K2;

  wmean3[0][0] = wP_P3;
  wmean3[0][1] = wP_AP3;
  wmean3[0][2] = wP_K3;

  wmean3[1][0] = wAP_P3;
  wmean3[1][1] = wAP_AP3;
  wmean3[1][2] = wAP_K3;

  wmean3[2][0] = wK_P3;
  wmean3[2][1] = wK_AP3;
  wmean3[2][2] = wK_K3;

  wmix[0][0] = wPAP_P;
  wmix[0][1] = wPAP_AP;
  wmix[0][2] = wPAP_K;

  wmix[1][0] = wPK_P;
  wmix[1][1] = wPK_AP;
  wmix[1][2] = wPK_K;

  wmix[2][0] = wAPK_P;
  wmix[2][1] = wAPK_AP;
  wmix[2][2] = wAPK_K;

  Float_t wmix2A[3][3];

  wmix2A[0][0] = wP2AP_P;
  wmix2A[0][1] = wP2AP_AP;
  wmix2A[0][2] = wP2AP_K;

  wmix2A[1][0] = wP2K_P;
  wmix2A[1][1] = wP2K_AP;
  wmix2A[1][2] = wP2K_K;

  wmix2A[2][0] = wAP2K_P;
  wmix2A[2][1] = wAP2K_AP;
  wmix2A[2][2] = wAP2K_K;

  Float_t wmix2B[3][3];

  wmix2B[0][0] = wAP2P_P;
  wmix2B[0][1] = wAP2P_AP;
  wmix2B[0][2] = wAP2P_K;

  wmix2B[1][0] = wK2P_P;
  wmix2B[1][1] = wK2P_AP;
  wmix2B[1][2] = wK2P_K;

  wmix2B[2][0] = wK2AP_P;
  wmix2B[2][1] = wK2AP_AP;
  wmix2B[2][2] = wK2AP_K;

  Float_t wmeanPrAPK[3];

  wmeanPrAPK[0] = wAPPrK_Pr;
  wmeanPrAPK[1] = wAPPrK_AP;
  wmeanPrAPK[2] = wAPPrK_K;

  Int_t npart = 3;
  Int_t nn = 0;
  for (Int_t p = 0; p < npart; p++)
  {
    nn = 0;
    for (Int_t i = 0; i < npart; i++)
    {
      A(p, nn) = TMath::Power(wmean[p][i], 3);
      nn++;
    }

    for (Int_t i = 0; i < (npart - 1); i++)
      for (Int_t l = i + 1; l < npart; l++)
      {
        A(p, nn) = 3. * TMath::Power(wmean[p][i], 2) * wmean[p][l];
        nn++;
      }

    for (Int_t i = 0; i < (npart - 1); i++)
      for (Int_t l = i + 1; l < npart; l++)
      {
        A(p, nn) = 3. * TMath::Power(wmean[p][l], 2) * wmean[p][i];
        nn++;
      }

    A(p, nn) = 6. * wmean[p][0] * wmean[p][1] * wmean[p][2];
    nn++;
  }

  //cout<<"nn == "<<nn<<endl;
  for (Int_t p = 0; p < npart; p++)
  {
    nn = 0;
    for (Int_t i = 0; i < npart; i++)
    {
      BB(p, nn) = 3. * (wmean2[p][i] * wmean[p][i] - TMath::Power(wmean[p][i], 3));
      nn++;
    }
    for (Int_t i = 0; i < (npart - 1); i++)
      for (Int_t l = i + 1; l < npart; l++)
      {
        BB(p, nn) = 3. * (wmean2[p][i] * wmean[p][l] + wmean2[p][l] * wmean[p][i] - TMath::Power(wmean[p][i], 2) * wmean[p][l] - TMath::Power(wmean[p][l], 2) * wmean[p][i]);
        nn++;
      }
    for (Int_t i = 0; i < npart; i++)
    {
      BB(p, nn) = 2. * TMath::Power(wmean[p][i], 3) + wmean3[p][i] - 3. * wmean2[p][i] * wmean[p][i];
      nn++;
    }
  }
  //cout<<"nn == "<<nn<<endl;

  //////////////////////
  //////////////////////
  //////////////////////

  nn = 0;
  for (Int_t i = 0; i < 3; i++)
  {
    A(3, nn) = wmean[0][i] * wmean[1][i] * wmean[2][i];
    nn++;
  }

  for (Int_t i = 0; i < (npart - 1); i++)
    for (Int_t l = i + 1; l < npart; l++)
    {
      A(3, nn) = wmean[0][i] * wmean[1][i] * wmean[2][l] + wmean[0][i] * wmean[1][l] * wmean[2][i] + wmean[0][l] * wmean[1][i] * wmean[2][i];
      nn++;
    }

  for (Int_t i = 0; i < (npart - 1); i++)
    for (Int_t l = i + 1; l < npart; l++)
    {
      A(3, nn) = wmean[0][l] * wmean[1][l] * wmean[2][i] + wmean[0][l] * wmean[1][i] * wmean[2][l] + wmean[0][i] * wmean[1][l] * wmean[2][l];
      nn++;
    }

  A(3, nn) = wmean[0][0] * wmean[1][1] * wmean[2][2] + wmean[0][1] * wmean[1][2] * wmean[2][0] + wmean[0][2] * wmean[1][0] * wmean[2][1] +
             wmean[0][0] * wmean[1][2] * wmean[2][1] + wmean[0][1] * wmean[1][0] * wmean[2][2] + wmean[0][2] * wmean[1][1] * wmean[2][0];

  nn = 0;
  Int_t p = 3;
  for (Int_t i = 0; i < npart; i++)
  {
    BB(3, nn) = wmix[0][i] * wmean[2][i] + wmix[1][i] * wmean[1][i] + wmix[2][i] * wmean[0][i] - 3. * wmean[0][i] * wmean[1][i] * wmean[2][i];
    nn++;
  }

  for (Int_t i = 0; i < (npart - 1); i++)
    for (Int_t l = i + 1; l < npart; l++)
    {
      BB(3, nn) = getW2(0, 0, 1, 2, i, l) + getW2(0, 0, 1, 2, l, i) + getW2(1, 0, 2, 1, i, l) + getW2(1, 0, 2, 1, l, i) +
                  getW2(2, 1, 2, 0, i, l) + getW2(2, 1, 2, 0, l, i);
      nn++;
    }

  for (Int_t i = 0; i < npart; i++)
  {
    BB(p, nn) = wmeanPrAPK[i] + 2. * wmean[0][i] * wmean[1][i] * wmean[2][i] - wmix[0][i] * wmean[2][i] - wmix[1][i] * wmean[1][i] - wmix[2][i] * wmean[0][i];
    nn++;
  }

  ////////////////////////////////////////////////////
  ////////////////////////////////////////////////////

  Int_t mm = 4;
  for (Int_t p = 0; p < (npart - 1); p++)
    for (Int_t q = p + 1; q < npart; q++)
    {
      nn = 0;
      for (Int_t i = 0; i < 3; i++)
      {
        A(mm, nn) = TMath::Power(wmean[p][i], 2) * wmean[q][i];
        nn++;
        // cout<<"mm "<<mm<<endl;
      }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          A(mm, nn) = TMath::Power(wmean[p][i], 2) * wmean[q][l] + 2. * wmean[p][i] * wmean[p][l] * wmean[q][i];
          nn++;
        }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          A(mm, nn) = TMath::Power(wmean[p][l], 2) * wmean[q][i] + 2. * wmean[p][i] * wmean[p][l] * wmean[q][l];
          nn++;
        }

      A(mm, nn) = 2. * (wmean[p][0] * wmean[p][1] * wmean[q][2] + wmean[p][0] * wmean[p][2] * wmean[q][1] + wmean[p][1] * wmean[p][2] * wmean[q][0]);

      nn = 0;
      for (Int_t i = 0; i < 3; i++)
      {

        BB(mm, nn) = wmean2[p][i] * wmean[q][i] - 3. * TMath::Power(wmean[p][i], 2) * wmean[q][i] + 2. * wmix[mm - 4][i] * wmean[p][i];
        nn++;
      }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          BB(mm, nn) = (2. * getW2(mm - 4, p, q, p, i, l) + 2. * getW2(mm - 4, p, q, p, l, i) + wmean2[p][i] * wmean[q][l] + wmean2[p][l] * wmean[q][i] -
                        TMath::Power(wmean[p][i], 2) * wmean[q][l] - TMath::Power(wmean[p][l], 2) * wmean[q][i]);
          nn++;
        }

      for (Int_t i = 0; i < 3; i++)
      {
        BB(mm, nn) = 2. * TMath::Power(wmean[p][i], 2) * wmean[q][i] + wmix2A[mm - 4][i] - wmean2[p][i] * wmean[q][i] - 2. * wmix[mm - 4][i] * wmean[p][i];
        nn++;
      }
      mm++;
    }
  ////////////////////////////////////////////////////

  mm = 7;
  for (Int_t q = 0; q < (npart - 1); q++)
    for (Int_t p = q + 1; p < npart; p++)
    {
      nn = 0;
      for (Int_t i = 0; i < 3; i++)
      {
        A(mm, nn) = TMath::Power(wmean[p][i], 2) * wmean[q][i];
        nn++;
      }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          A(mm, nn) = TMath::Power(wmean[p][i], 2) * wmean[q][l] + 2. * wmean[p][i] * wmean[p][l] * wmean[q][i];
          nn++;
        }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          A(mm, nn) = TMath::Power(wmean[p][l], 2) * wmean[q][i] + 2. * wmean[p][i] * wmean[p][l] * wmean[q][l];
          nn++;
        }

      A(mm, nn) = 2. * (wmean[p][0] * wmean[p][1] * wmean[q][2] + wmean[p][0] * wmean[p][2] * wmean[q][1] + wmean[p][1] * wmean[p][2] * wmean[q][0]);

      nn = 0;
      for (Int_t i = 0; i < 3; i++)
      {

        BB(mm, nn) = wmean2[p][i] * wmean[q][i] - 3. * TMath::Power(wmean[p][i], 2) * wmean[q][i] + 2. * wmix[mm - 7][i] * wmean[p][i];
        nn++;
      }

      for (Int_t i = 0; i < (npart - 1); i++)
        for (Int_t l = i + 1; l < npart; l++)
        {
          // cout << "nn == " << mm << "  " << nn << endl;
          BB(mm, nn) = (2. * getW2(mm - 7, p, q, p, i, l) + 2. * getW2(mm - 7, p, q, p, l, i) + wmean2[p][i] * wmean[q][l] + wmean2[p][l] * wmean[q][i] - TMath::Power(wmean[p][i], 2) * wmean[q][l] - TMath::Power(wmean[p][l], 2) * wmean[q][i]);
          nn++;
        }

      for (Int_t i = 0; i < 3; i++)
      {
        BB(mm, nn) = 2. * TMath::Power(wmean[p][i], 2) * wmean[q][i] + wmix2B[mm - 7][i] - wmean2[p][i] * wmean[q][i] - 2. * wmix[mm - 7][i] * wmean[p][i];
        nn++;
      }
      mm++;
    }

  Double_t B[10] = {W3P_aver, W3AP_aver, W3K_aver, WAPPrK_aver, WPr2AP_aver,
                    WPr2K_aver, WAP2K_aver, WAP2Pr_aver, WK2Pr_aver, WK2AP_aver};

  for (Int_t i = 0; i < 10; i++)
  {
    B[i] -= (BB(i, 0) * pr2_aver + BB(i, 1) * pi2_aver +
             BB(i, 2) * k2_aver + BB(i, 3) * prapr_aver +
             BB(i, 4) * prk_aver + BB(i, 5) * pik_aver +
             BB(i, 6) * pr_aver + BB(i, 7) * apr_aver + BB(i, 8) * k_aver);
  }

  // TString mom3Names[] = {"P3", "AP3", "K3", "P2AP", "P2K", "AP2K", "AP2P", "K2P", "K2AP", "APPK"};
  TString mom3Names[] = {"pr3", "apr3", "bg3", "pr2_apr", "pr2_bg", "apr2_bg", "apr2_pr", "bg2_pr", "bg2_apr", "apr_pr_bg"};

  TMatrixD invA = A.Invert();
  Double_t rec3mom[10] = {0};
  cout << " ==================================" << endl;
  cout << " ==================================" << endl;
  cout << " print third moments " << endl;
  for (Int_t i = 0; i < 10; i++)
  {
    for (Int_t tt = 0; tt < 10; tt++)
    {
      rec3mom[i] += invA(i, tt) * B[tt];
    }
    // N^3+3N^2+N, (A*A+A)*B --> poisson 3rd moment
    cout << " " << mom3Names[i].Data() << "   :  "  << rec3mom[i] << endl;
    (*fMoments3rd)[i]=rec3mom[i];
  }
  cout << " ================================== " << endl;
  cout << " ================================== " << endl;
  cout << " calculation is over " << endl;
  timer.Stop(); timer.Print();
  cout << " ================================== " << endl;
  cout << " ================================== " << endl;

  //cout<<"rec3mom 0 "<<rec3mom[0]<<endl;
  Double_t skew = rec3mom[0] - 3. * pr2_aver * proton_aver + 2. * TMath::Power(proton_aver, 3);

  // cout<<"third variance "<<skew/proton_aver<<endl;
  cout << " =========================== " << endl;
  cout << "ratio to poisson 3rd moment  " << (TMath::Power(proton_aver, 3) + 3. * TMath::Power(proton_aver, 2) + proton_aver) / rec3mom[0] << endl;
  cout << "ratio to poisson 2nd moment  " << (TMath::Power(proton_aver, 2) + proton_aver) / pr2_aver << endl;
  cout << " =========================== " << endl;

  skew /= (pr2_aver - proton_aver * proton_aver);
  cout << "skew === " << skew << endl;
  momTree->Fill();
  outFile->cd();
  momTree -> Write();
  outFile -> Close();
  delete outFile;
  return 0;
}

int main(int argc, char *argv[])
{
  cout << "      " << endl;
  cout << "***********************************************" << endl;
  cout << "***********************************************" << endl;
  cout << "******************** IDENTITY METHOD **********" << endl;
  cout << "***********************************************" << endl;
  cout << "***********************************************" << endl;
  cout << "      " << endl;
  //
  if(argc == 9)
  {
    sprintf(inputfileNameDataTree,"%s",argv[1]);
    sprintf(inputfileNameLineShapes,"%s",argv[2]);
    fSubsample       = atoi(argv[3]);
    fCentInput       = atof(argv[4]);
    fpDownInput      = atof(argv[5]);
    fpUpInput        = atof(argv[6]);
    fEtaDownInput    = atof(argv[7]);
    fEtaUpInput      = atof(argv[8]);
    cout<<" main.Info: read file names from input "<<endl;
    TString outPutFileNAme = Form("TIMoments3D_sub%d_cent_%3.2f_mom_%3.2f_%3.2f_eta_%3.2f_%3.2f.root",fSubsample,fCentInput,fpDownInput,fpUpInput,fEtaDownInput,fEtaUpInput);
    outFile = new TFile(outPutFileNAme,"recreate");
    cout << " main.Info: write output into:    " << outPutFileNAme << endl;
  }
  else
  {
    cout<<" main.Error: wrong input list"<<endl;
  }
  //
  //
  InitializeObjects();
  fileNameDataTree   = inputfileNameDataTree;
  fileNameLineShapes = inputfileNameLineShapes;
  fpDownBin       = fhPtot -> FindBin(fpDownInput   + 0.0000001) - 1;
  fpUpBin         = fhPtot -> FindBin(fpUpInput     - 0.0000001) - 1;
  fEtaDownBin     = fhEta  -> FindBin(fEtaDownInput + 0.0000001) - 1;
  fEtaUpBin       = fhEta  -> FindBin(fEtaUpInput   - 0.0000001) - 1;
  fCentInputBin   = fhCent -> FindBin(fCentInput    + 0.0000001) - 1;
  fCentBinCenter  = fhCent -> GetXaxis()->GetBinCenter(fCentInputBin+1);
  //
  //
  TROOT IdentityMethod("IdentityMethod", "compiled identity method");
  ReadFitParamsFromTree(fileNameLineShapes,6);
  identity3Particle();
  return 0;
}
