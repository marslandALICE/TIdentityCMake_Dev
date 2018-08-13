
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

TF1 *funProton, *funKaon, *funPion;

Double_t evtNumber;
Double_t WK_sum, WP_sum, WPi_sum, WMult_sum;
Int_t multEv;

vector<int> sumMult;

vector<double> W2P_sum_vec;

vector<double> W3P_sum_vec;
vector<double> W3Pi_sum_vec;
vector<double> W3K_sum_vec;

vector<double> W2K_sum_vec;
vector<double> W2Pi_sum_vec;
vector<double> WPK_sum_vec;
vector<double> WPPi_sum_vec;
vector<double> WPiK_sum_vec;
vector<double> WPiPrK_sum_vec;

vector<double> WPr2Pi_sum_vec;
vector<double> WPr2K_sum_vec;
vector<double> WPi2K_sum_vec;

vector<double> WPi2Pr_sum_vec;
vector<double> WK2Pr_sum_vec;
vector<double> WK2Pi_sum_vec;

Double_t recP2_av, recPi2_av, recK2_av, recPPi_av, recPK_av, recPiK_av;
Double_t wP_P, wK_P, wPi_P, wP_P2, wK_P2, wPi_P2, wP_P3, wK_P3, wPi_P3, wP_K, wK_K, wPi_K, wP_K2, wK_K2, wPi_K2;
Double_t wP_K3, wK_K3, wPi_K3, wP_Pi, wK_Pi, wPi_Pi, wP_Pi2, wK_Pi2, wPi_Pi2, wP_Pi3, wK_Pi3, wPi_Pi3;

///mixed
Double_t wPK_P, wPPi_P, wPiK_P, wPK_K, wPPi_K, wPiK_K, wPK_Pi, wPPi_Pi, wPiK_Pi, wP2Pi_P, wP2Pi_Pi, wP2Pi_K, wPi2P_P, wPi2P_Pi;
Double_t wPi2P_K, wP2K_P, wP2K_Pi, wP2K_K, wK2P_P, wK2P_Pi, wK2P_K, wPi2K_P, wPi2K_Pi, wPi2K_K, wK2Pi_P, wK2Pi_Pi, wK2Pi_K, wPiPrK_Pr, wPiPrK_K, wPiPrK_Pi;
//
TF1 *funP_P = NULL,    *funK_P = NULL,     *funPi_P = NULL,    *funP_P2 = NULL,   *funK_P2 = NULL,    *funPi_P2 = NULL,  *funP_P3 = NULL,   *funK_P3 = NULL,   *funPi_P3 = NULL;
TF1 *funP_K = NULL,    *funK_K = NULL,     *funPi_K = NULL,    *funP_K2 = NULL,   *funK_K2 = NULL,    *funPi_K2 = NULL,  *funP_K3 = NULL,   *funK_K3 = NULL,   *funPi_K3 = NULL;
TF1 *funP_Pi = NULL,   *funK_Pi = NULL,    *funPi_Pi = NULL,   *funP_Pi2 = NULL,  *funK_Pi2 = NULL,   *funPi_Pi2 = NULL, *funP_Pi3 = NULL,  *funK_Pi3 = NULL,  *funPi_Pi3 = NULL;
TF1 *funPK_P = NULL,   *funPPi_P = NULL,   *funPiK_P = NULL,   *funPK_K = NULL,   *funPPi_K = NULL,   *funPiK_K = NULL,  *funPK_Pi = NULL,  *funPPi_Pi = NULL, *funPiK_Pi = NULL;
TF1 *funP2Pi_P = NULL, *funP2Pi_Pi = NULL, *funP2Pi_K = NULL,  *funPi2P_P = NULL, *funPi2P_Pi = NULL, *funPi2P_K = NULL, *funP2K_P = NULL,  *funP2K_Pi = NULL, *funP2K_K = NULL;
TF1 *funK2P_P = NULL,  *funK2P_Pi = NULL,  *funK2P_K = NULL,   *funPi2K_P = NULL, *funPi2K_Pi = NULL, *funPi2K_K = NULL, *funK2Pi_P = NULL, *funK2Pi_Pi = NULL;
TF1 *funK2Pi_K = NULL, *funPiPrK_P = NULL, *funPiPrK_K = NULL, *funPiPrK_Pi = NULL;
////////////////////////////

Int_t size_size = 3000;
Int_t fMyBin[3];

Float_t wmean[3][3];
Float_t wmean2[3][3];
Float_t wmean3[3][3];
Float_t wmix[3][3];

vector<double> WPM_sum_vec;
vector<double> WKM_sum_vec;
vector<double> WPiM_sum_vec;

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
//
const Int_t nBinsLineShape = 4080; // 6120
Int_t fnTestEntries = 0;
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
TH1D *hDedxDebug;

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
TVectorF *fMoments1st,  *fMoments2nd, *fMoments2ndMixed;
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
Float_t meanKaon[fnEtaBins][fnCentBins][fnMomBins];
Float_t meanProton[fnEtaBins][fnCentBins][fnMomBins];
Float_t meanPion[fnEtaBins][fnCentBins][fnMomBins];
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
  fIntegrals    = new TVectorF(fnParticleBins);
  fMoments1st   = new TVectorF(fnParticleBins);
  fMoments2nd   = new TVectorF(fnParticleBins);
  fMoments2ndMixed   = new TVectorF(fnParticleBins*fnParticleBins);

  for(Int_t i=0;i<fnParticleBins*fnParticleBins; i++){
    if (i<fnParticleBins) {
      (*fIntegrals)[i]=0.;
      (*fMoments1st)[i]=0.;
      (*fMoments2nd)[i]=0.;
    }
    (*fMoments2ndMixed)[i]=0.;
  }
  //
  // initialize output tree
  //
  momTree = new TTree("momTree","momTree");
  momTree -> Branch("nEvents",&nEvents);
  momTree -> Branch("nnorm",&nnorm);
  momTree -> Branch("pDown",&fpDownInput);
  momTree -> Branch("pUp",&fpUpInput);
  momTree -> Branch("etaDown",&fEtaDownInput);
  momTree -> Branch("etaUp",&fEtaUpInput);
  momTree -> Branch("cent",&fCentInput);
  momTree -> Branch("centBin",&fCentInputBin);
  momTree -> Branch("fIntegrals",&fIntegrals);
  momTree -> Branch("fMoments1st",&fMoments1st);
  momTree -> Branch("fMoments2nd",&fMoments2nd);
  momTree -> Branch("fMoments2ndMixed",&fMoments2ndMixed);

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

      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBEl] = elSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBPi] = piSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBKa] = kaSk;
      fSkewArr[myBin[0]][myBin[1]][myBin[2]][kBPr] = prSk;

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
        hParticles[i][j][k][2]->SetName(objOthersName);
        hParticles[i][j][k][2]->SetLineColor(kBlack);

        if (i==etaBinRange[0] && j==centBinRange[0]) {
          hParticles[i][j][k][0]->Write();
          hParticles[i][j][k][1]->Write();
          hParticles[i][j][k][2]->Write();
        }

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
  // cout << "  aaaa ->  " <<  xval[0] << "  "  << par[0] << " " << fEtaBin << " " <<  fCentBin << " " << fMomBin << "  " << hParticles[fEtaBin][fCentBin][fMomBin][pp]->GetName() << endl;
  // if (pp==0) cout << xx << "   -------  " <<  hParticles[fEtaBin][fCentBin][fMomBin][0]->FindBin(xx) << "  " << fEtaBin << " " <<  fCentBin << " " << fMomBin << "  hhhhhh  " << endl;
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
    // cout << xx << "  jjj " << step << "  " << i << "  " << fMaxdEdx - fMindEdx << "  " << 2 * size_size << endl;
    sum += fun->Eval(xx);
  }
  return sum * step;
}

void initFunctions_M()
{

  if (funProton) delete funProton;
  funProton = new TF1("funProton", getValueH, fMindEdx, fMaxdEdx, 1);
  funProton->SetParameter(0, 0);

  if (funPion)  delete funPion;
  funPion = new TF1("funPion", getValueH, fMindEdx, fMaxdEdx, 1);
  funPion->SetParameter(0, 1);

  if (funKaon)  delete funKaon;
  funKaon = new TF1("funKaon", getValueH, fMindEdx, fMaxdEdx, 1);
  funKaon->SetParameter(0, 2);

}

Double_t getFunctions(Double_t *xx, Double_t *par)
{
  Double_t x = xx[0];
  Int_t j = (Int_t)par[0];
  Int_t i = (Int_t)par[1];
  Int_t k = (Int_t)par[2];
  Double_t val[3];
  val[0] = funProton->Eval(x);
  val[1] = funKaon->Eval(x);
  val[2] = funPion->Eval(x);
  // cout << "oooo " << val[0] << "  " << val[1] << "  " << val[2] << endl;
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
  val[1] = funKaon->Eval(x);
  val[2] = funPion->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15)
    return 0.;

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
  val[1] = funKaon->Eval(x);
  val[2] = funPion->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15)
    return 0.;

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
  val[1] = funKaon->Eval(x);
  val[2] = funPion->Eval(x);
  Double_t sumVal = val[0] + val[1] + val[2];
  if (sumVal < 1e-15)
    return 0.;
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
  funPi_P = new TF1("funPi_P", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_P, 2, 0, 1);
  funP_P2 = new TF1("funP_P2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_P2, 0, 0, 2);
  funK_P2 = new TF1("funK_P2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_P2, 1, 0, 2);
  funPi_P2 = new TF1("funPi_P2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_P2, 2, 0, 2);
  funP_P3 = new TF1("funP_P3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_P3, 0, 0, 3);
  funK_P3 = new TF1("funK_P3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_P3, 1, 0, 3);
  funPi_P3 = new TF1("funPi_P3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_P3, 2, 0, 3);
  funP_K = new TF1("funP_K", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_K, 0, 1, 1);
  funK_K = new TF1("funK_K", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_K, 1, 1, 1);
  funPi_K = new TF1("funPi_K", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_K, 2, 1, 1);
  funP_K2 = new TF1("funP_K2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_K2, 0, 1, 2);
  funK_K2 = new TF1("funK_K2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_K2, 1, 1, 2);
  funPi_K2 = new TF1("funPi_K2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_K2, 2, 1, 2);
  funP_K3 = new TF1("funP_K3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_K3, 0, 1, 3);
  funK_K3 = new TF1("funK_K3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_K3, 1, 1, 3);
  funPi_K3 = new TF1("funPi_K3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_K3, 2, 1, 3);
  funP_Pi = new TF1("funP_Pi", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_Pi, 0, 2, 1);
  funK_Pi = new TF1("funK_Pi", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_Pi, 1, 2, 1);
  funPi_Pi = new TF1("funPi_Pi", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_Pi, 2, 2, 1);
  funP_Pi2 = new TF1("funP_Pi2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_Pi2, 0, 2, 2);
  funK_Pi2 = new TF1("funK_Pi2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_Pi2, 1, 2, 2);
  funPi_Pi2 = new TF1("funPi_Pi2", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_Pi2, 2, 2, 2);
  funP_Pi3 = new TF1("funP_Pi3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funP_Pi3, 0, 2, 3);
  funK_Pi3 = new TF1("funK_Pi3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funK_Pi3, 1, 2, 3);
  funPi_Pi3 = new TF1("funPi_Pi3", getFunctions, fMindEdx, fMaxdEdx, 3);
  setpars(funPi_Pi3, 2, 2, 3);
  ////////////////// mixed ones
  funPK_P = new TF1("funPK_P", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPK_P, 0, 0, -10);
  funPPi_P = new TF1("funPPi_P", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPPi_P, 1, 0, -10);
  funPiK_P = new TF1("funPiK_P", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPiK_P, 2, 0, -10);
  funPK_K = new TF1("funPK_K", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPK_K, 0, 1, -10);
  funPPi_K = new TF1("funPPi_K", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPPi_K, 1, 1, -10);
  funPiK_K = new TF1("funPiK_K", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPiK_K, 2, 1, -10);
  funPK_Pi = new TF1("funPK_Pi", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPK_Pi, 0, 2, -10);
  funPPi_Pi = new TF1("funPPi_Pi", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPPi_Pi, 1, 2, -10);
  funPiK_Pi = new TF1("funPiK_Pi", getFunctionsMix, fMindEdx, fMaxdEdx, 3);
  setpars(funPiK_Pi, 2, 2, -10);

  funP2Pi_P = new TF1("funP2Pi_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2Pi_P, 0, 2, 0);
  funP2Pi_Pi = new TF1("funP2Pi_Pi", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2Pi_Pi, 0, 2, 2);
  funP2Pi_K = new TF1("funP2Pi_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2Pi_K, 0, 2, 1);
  funPi2P_P = new TF1("funPi2P_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funPi2P_P, 2, 0, 0);
  funPi2P_Pi = new TF1("funPi2P_Pi", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funPi2P_Pi, 2, 0, 2);
  funPi2P_K = new TF1("funPi2P_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funPi2P_K, 2, 0, 1);
  funP2K_P = new TF1("funP2K_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2K_P, 0, 1, 0);
  funP2K_Pi = new TF1("funP2K_Pi", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2K_Pi, 0, 1, 2);
  funP2K_K = new TF1("funP2K_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funP2K_K, 0, 1, 1);
  funK2P_P = new TF1("funK2P_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2P_P, 1, 0, 0);
  funK2P_Pi = new TF1("funK2P_Pi", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2P_Pi, 1, 0, 2);
  funK2P_K = new TF1("funK2P_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2P_K, 1, 0, 1);
  funPi2K_P = new TF1("funPi2K_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funPi2K_P, 2, 1, 0);
  funPi2K_Pi = new TF1("funPi2K_Pi", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funPi2K_Pi, 2, 1, 2);
  funPi2K_K = new TF1("funPi2K_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funPi2K_K, 2, 1, 1);
  funK2Pi_P = new TF1("funK2Pi_P", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2Pi_P, 1, 2, 0);
  funK2Pi_Pi = new TF1("funK2Pi_Pi", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2Pi_Pi, 1, 2, 2);
  funK2Pi_K = new TF1("funK2Pi_K", getFunctionsMix2, fMindEdx, fMaxdEdx, 3);
  setpars(funK2Pi_K, 1, 2, 1);
  funPiPrK_P = new TF1("funPiPrK_P", getFunctionsMix3, fMindEdx, fMaxdEdx, 3);
  setpars(funPiPrK_P, 0, -100, -100);
  funPiPrK_K = new TF1("funPiPrK_K", getFunctionsMix3, fMindEdx, fMaxdEdx, 3);
  setpars(funPiPrK_K, 1, -100, -100);
  funPiPrK_Pi = new TF1("funPiPrK_Pi", getFunctionsMix3, fMindEdx, fMaxdEdx, 3);
  setpars(funPiPrK_Pi, 2, -100, -100);
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
  wPi_P += myIntegral(funPi_P) / npNorm;

  wP_P2 += myIntegral(funP_P2) / npNorm;
  wK_P2 += myIntegral(funK_P2) / npNorm;
  wPi_P2 += myIntegral(funPi_P2) / npNorm;

  wP_P3 += myIntegral(funP_P3) / npNorm;
  wK_P3 += myIntegral(funK_P3) / npNorm;
  wPi_P3 += myIntegral(funPi_P3) / npNorm;

  wP_K += myIntegral(funP_K) / nkNorm;
  wK_K += myIntegral(funK_K) / nkNorm;
  wPi_K += myIntegral(funPi_K) / nkNorm;

  wP_K2 += myIntegral(funP_K2) / nkNorm;
  wK_K2 += myIntegral(funK_K2) / nkNorm;
  wPi_K2 += myIntegral(funPi_K2) / nkNorm;

  wP_K3 += myIntegral(funP_K3) / nkNorm;
  wK_K3 += myIntegral(funK_K3) / nkNorm;
  wPi_K3 += myIntegral(funPi_K3) / nkNorm;

  wP_Pi += myIntegral(funP_Pi) / npiNorm;
  wK_Pi += myIntegral(funK_Pi) / npiNorm;
  wPi_Pi += myIntegral(funPi_Pi) / npiNorm;

  wP_Pi2 += myIntegral(funP_Pi2) / npiNorm;
  wK_Pi2 += myIntegral(funK_Pi2) / npiNorm;
  wPi_Pi2 += myIntegral(funPi_Pi2) / npiNorm;

  wP_Pi3 += myIntegral(funP_Pi3) / npiNorm;
  wK_Pi3 += myIntegral(funK_Pi3) / npiNorm;
  wPi_Pi3 += myIntegral(funPi_Pi3) / npiNorm;

  ///mixed
  wPK_P += myIntegral(funPK_P) / npNorm;
  wPPi_P += myIntegral(funPPi_P) / npNorm;
  wPiK_P += myIntegral(funPiK_P) / npNorm;

  wPK_K += myIntegral(funPK_K) / nkNorm;
  wPPi_K += myIntegral(funPPi_K) / nkNorm;
  wPiK_K += myIntegral(funPiK_K) / nkNorm;

  wPK_Pi += myIntegral(funPK_Pi) / npiNorm;
  wPPi_Pi += myIntegral(funPPi_Pi) / npiNorm;
  wPiK_Pi += myIntegral(funPiK_Pi) / npiNorm;

  ///////////////////////////////////////////////////////////////////////

  wP2Pi_P += myIntegral(funP2Pi_P) / npNorm;
  wP2Pi_Pi += myIntegral(funP2Pi_Pi) / npiNorm;
  wP2Pi_K += myIntegral(funP2Pi_K) / nkNorm;

  wPi2P_P += myIntegral(funPi2P_P) / npNorm;
  wPi2P_Pi += myIntegral(funPi2P_Pi) / npiNorm;
  wPi2P_K += myIntegral(funPi2P_K) / nkNorm;

  wP2K_P += myIntegral(funP2K_P) / npNorm;
  wP2K_Pi += myIntegral(funP2K_Pi) / npiNorm;
  wP2K_K += myIntegral(funP2K_K) / nkNorm;

  wK2P_P += myIntegral(funK2P_P) / npNorm;
  wK2P_Pi += myIntegral(funK2P_Pi) / npiNorm;
  wK2P_K += myIntegral(funK2P_K) / nkNorm;

  wPi2K_P += myIntegral(funPi2K_P) / npNorm;
  wPi2K_Pi += myIntegral(funPi2K_Pi) / npiNorm;
  wPi2K_K += myIntegral(funPi2K_K) / nkNorm;

  wK2Pi_P += myIntegral(funK2Pi_P) / npNorm;
  wK2Pi_Pi += myIntegral(funK2Pi_Pi) / npiNorm;
  wK2Pi_K += myIntegral(funK2Pi_K) / nkNorm;

  wPiPrK_Pr += myIntegral(funPiPrK_P) / npNorm;
  wPiPrK_K += myIntegral(funPiPrK_K) / nkNorm;
  wPiPrK_Pi += myIntegral(funPiPrK_Pi) / npiNorm;
}

/////////////////////////////

template <class T>
Float_t Sum(vector<T> &v)
{
  return accumulate(v.begin(), v.end(), 0.0);
}

void resetValues()
{
  WP_sum = WK_sum = WPi_sum = WMult_sum = 0.;
}

void addParticles(Float_t mVal, Int_t &count)
{
  Double_t kValue, prValue, piValue;
  Double_t WP1, WK1, WPi1;
  //
  prValue = funProton->Eval(mVal);
  kValue  = funKaon->Eval(mVal);
  piValue = funPion->Eval(mVal);

  Double_t sumValue = prValue + kValue + piValue;
  if (sumValue == 0)
  {
    cout << "why:   " << mVal << endl;
    wrongCount++;
    WP1 = WK1 = WPi1 = 0;
  }
  else
  {
    WP1 = prValue / sumValue;
    WK1 = kValue / sumValue;
    WPi1 = piValue / sumValue;

    count++;

    WP_sum += WP1;
    WK_sum += WK1;
    WPi_sum += WPi1;
    WMult_sum += 1;
  }
}

int identity3Particle()
{

  cout << "burada " << endl;
  TFile *dataFile = new TFile(fileNameDataTree);
  wrongCount = 0;
  TTree *myTree = (TTree *)dataFile->Get(treeIdentity);
  cout << "Funtions are inited" << endl;
  Int_t prevEvt = -1000;
  cout << "set address" << endl;
  //
  Int_t fEventNum;  // ULong64_t fEventNum;
  Double_t fDEdx;    // Float_t fDEdx;
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
        meanKaon[i][j][k] = 0;
        meanProton[i][j][k] = 0;
        meanPion[i][j][k] = 0;
        fUsedBins[i][j][k] = 0;
      }
    }
  }

  numAllEvents = 0;
  numAllCutEvents = 0;

  Int_t runEvt = (Int_t)myTree->GetEntries();
  cout << "total number " << runEvt << endl;
  Int_t runStat1 = 0;
  Int_t runStat2 = runEvt;

  cout << "running from " << runStat1 << " to " << runStat2 << endl;
  Int_t countVeto = 0;
  Int_t prevEvtVeto = -1;
  initFunctions_M();

  myTree->GetEntry(runStat2 - 1);
  Int_t remEvent = fEventNum;
  cout << " ==================================" << endl;
  cout << " loop over tracks starts " <<endl;
  timer.Reset(); timer.Start();
  cout << " ==================================" << endl;
  for (Int_t i = runStat1; i < runStat2; i++)
  {
    if (i % 2000000 == 0) cout << "track " << i << " of " << runStat2 - runStat1 << endl;
    if (i>100000) break;
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

    fUsedBins[fEtaBin][fCentBin][fMomBin] = 1;
    // cout << " aaaa " << hParticles[fEtaBin][fCentBin][fMomBin][0]->GetMean() << "   " <<  fEtaBin << "  " << fCentBin << "  " << fMomBin << endl;
    if (!meanProton[fEtaBin][fCentBin][fMomBin]) meanProton[fEtaBin][fCentBin][fMomBin] = myIntegral(funProton);
    if (!meanPion[fEtaBin][fCentBin][fMomBin])   meanPion[fEtaBin][fCentBin][fMomBin]   = myIntegral(funPion);
    if (!meanKaon[fEtaBin][fCentBin][fMomBin])   meanKaon[fEtaBin][fCentBin][fMomBin]   = myIntegral(funKaon);

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
          W3Pi_sum_vec.push_back(TMath::Power(WPi_sum, 3));
          W3K_sum_vec.push_back(TMath::Power(WK_sum, 3));

          W2P_sum_vec.push_back(WP_sum * WP_sum);
          W2K_sum_vec.push_back(WK_sum * WK_sum);
          W2Pi_sum_vec.push_back(WPi_sum * WPi_sum);

          WPM_sum_vec.push_back(WP_sum);
          WKM_sum_vec.push_back(WK_sum);
          WPiM_sum_vec.push_back(WPi_sum);

          WPK_sum_vec.push_back(WP_sum * WK_sum);
          WPPi_sum_vec.push_back(WP_sum * WPi_sum);
          WPiK_sum_vec.push_back(WPi_sum * WK_sum);

          WPiPrK_sum_vec.push_back(WPi_sum * WP_sum * WK_sum);

          WPr2Pi_sum_vec.push_back(TMath::Power(WP_sum, 2) * WPi_sum);
          WPr2K_sum_vec.push_back(TMath::Power(WP_sum, 2) * WK_sum);
          WPi2K_sum_vec.push_back(TMath::Power(WPi_sum, 2) * WK_sum);

          WPi2Pr_sum_vec.push_back(WPi_sum * WPi_sum * WP_sum);
          WK2Pr_sum_vec.push_back(WK_sum * WK_sum * WP_sum);
          WK2Pi_sum_vec.push_back(WK_sum * WK_sum * WPi_sum);
        }
      }
      resetValues();
      count = 0;
      addParticles(fDEdx, count);
      numAllCutEvents++;
    }
    prevEvt = fEventNum;
  } //end event
  cout << "====================================" << endl;
  cout << " track loop is over " << endl;
  timer.Stop(); timer.Print();
  cout << "====================================" << endl;
  //
  Double_t nEvents = countVeto;
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
  Double_t kaon_aver = 0;
  Double_t pion_aver = 0;

  for (Int_t i = 0; i < fnEtaBins; i++){
    for (Int_t j = 0; j < fnCentBins; j++){
      for (Int_t k = 0; k < fnMomBins; k++){
        if(fUsedBins[i][j][k] != 1) continue;
        // cout << i << "  " << j << "  " << k << "  " << meanKaon[i][j][k] << "  " << meanProton[i][j][k] << "  " <<  meanPion[i][j][k] << endl;
        kaon_aver += meanKaon[i][j][k];
        proton_aver += meanProton[i][j][k];
        pion_aver += meanPion[i][j][k];
      }
    }
  }

  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "mean multiplicities : ONLY TRUE FOR ALL STATISTICS !" << endl;
  cout << "proton " << proton_aver * corrFactor << endl;
  cout << "pkaon  " << kaon_aver * corrFactor << endl;
  cout << "pion   " << pion_aver * corrFactor << endl;

  proton_aver *= corrFactor;
  kaon_aver *= corrFactor;
  pion_aver *= corrFactor;

  cout << "***********" << endl;
  cout << "***********" << endl;
  cout << "***********" << endl;

  //Double_t averMult = Sum<double>(WMult_sum_vec)/nEvents;
  //Double_t averMult2 = Sum<double>(W2Mult_sum_vec)/nEvents;

  Double_t W2P_aver = Sum<double>(W2P_sum_vec) / nEvents;
  Double_t W2K_aver = Sum<double>(W2K_sum_vec) / nEvents;
  Double_t W2Pi_aver = Sum<double>(W2Pi_sum_vec) / nEvents;

  Double_t WPM_aver = Sum<double>(WPM_sum_vec) / nEvents;
  Double_t WKM_aver = Sum<double>(WKM_sum_vec) / nEvents;
  Double_t WPiM_aver = Sum<double>(WPiM_sum_vec) / nEvents;

  Double_t W3P_aver = Sum<double>(W3P_sum_vec) / nEvents;
  Double_t W3Pi_aver = Sum<double>(W3Pi_sum_vec) / nEvents;
  Double_t W3K_aver = Sum<double>(W3K_sum_vec) / nEvents;
  Double_t WPiPrK_aver = Sum<double>(WPiPrK_sum_vec) / nEvents;
  //    Double_t WPi2P_aver     = Sum<double>(WPi2P_sum_vec)/nEvents;

  Double_t WPr2Pi_aver = Sum<double>(WPr2Pi_sum_vec) / nEvents;
  Double_t WPr2K_aver = Sum<double>(WPr2K_sum_vec) / nEvents;
  Double_t WPi2K_aver = Sum<double>(WPi2K_sum_vec) / nEvents;

  Double_t WPi2Pr_aver = Sum<double>(WPi2Pr_sum_vec) / nEvents;
  Double_t WK2Pr_aver = Sum<double>(WK2Pr_sum_vec) / nEvents;
  Double_t WK2Pi_aver = Sum<double>(WK2Pi_sum_vec) / nEvents;

  Double_t WPK_aver = Sum<double>(WPK_sum_vec) / nEvents;
  Double_t WPPi_aver = Sum<double>(WPPi_sum_vec) / nEvents;
  Double_t WPiK_aver = Sum<double>(WPiK_sum_vec) / nEvents;

  Double_t pr_aver = 0.;  //Sum<int>(P_mult_vec)/nEvents;
  Double_t k_aver = 0.;   ///Sum<int>(K_mult_vec)/nEvents;
  Double_t pi_aver = 0.;  ///Sum<int>(Pi_mult_vec)/nEvents;

  Double_t pr2_aver = 0.; ////Sum<int>(P2_mult_vec)/nEvents;
  Double_t k2_aver = 0.;  ////Sum<int>(K2_mult_vec)/nEvents;
  Double_t pi2_aver = 0.; ///Sum<int>(Pi2_mult_vec)/nEvents;

  Double_t prpi_aver = 0.; ///Sum<int>(PPi_mult_vec)/nEvents;
  Double_t prk_aver = 0.;  ///Sum<int>(PK_mult_vec)/nEvents;
  Double_t pik_aver = 0.;  ////Sum<int>(PiK_mult_vec)/nEvents;

  cout << "W2P_aver " << W2P_aver << endl;
  cout << "W2K_aver " << W2K_aver << endl;
  cout << "W2Pi_aver " << W2Pi_aver << endl;
  cout << "WPK_aver " << WPK_aver << endl;
  cout << "WPPi_aver " << WPPi_aver << endl;
  cout << "W2PiK_aver " << WPiK_aver << endl;

  cout << "mean from identities " << endl;
  cout << "WPM_aver  " << WPM_aver  << " -> " << proton_aver << endl;
  cout << "WKM_aver  " << WKM_aver  << " -> " << kaon_aver   << endl;
  cout << "WPiM_aver " << WPiM_aver << " -> " << pion_aver   <<  endl;

  /*
   proton_aver   = WPM_aver;
   kaon_aver     = WKM_aver;
   pion_aver     = WPiM_aver;
   */

  cout << "start to calculate integrals " << endl;
  wP_P = 0., wK_P = 0., wPi_P = 0., wP_P2 = 0., wK_P2 = 0., wPi_P2 = 0., wP_P3 = 0., wK_P3 = 0., wPi_P3 = 0., wP_K = 0., wK_K = 0., wPi_K = 0., wP_K2 = 0., wK_K2 = 0.;
  wPi_K2 = 0., wP_K3 = 0., wK_K3 = 0., wPi_K3 = 0., wP_Pi = 0., wK_Pi = 0., wPi_Pi = 0., wP_Pi2 = 0., wK_Pi2 = 0., wPi_Pi2 = 0., wP_Pi3 = 0., wK_Pi3 = 0., wPi_Pi3 = 0.;
  ///mixed
  wPK_P = 0., wPPi_P = 0., wPiK_P = 0., wPK_K = 0., wPPi_K = 0., wPiK_K = 0., wPK_Pi = 0., wPPi_Pi = 0., wPiK_Pi = 0., wP2Pi_P = 0., wP2Pi_Pi = 0., wP2Pi_K = 0.;
  wPi2P_P = 0., wPi2P_Pi = 0., wPi2P_K = 0., wP2K_P = 0., wP2K_Pi = 0., wP2K_K = 0., wK2P_P = 0., wK2P_Pi = 0., wK2P_K = 0., wPi2K_P = 0., wPi2K_Pi = 0., wPi2K_K = 0.;
  wK2Pi_P = 0., wK2Pi_Pi = 0., wK2Pi_K = 0., wPiPrK_Pr = 0., wPiPrK_K = 0., wPiPrK_Pi = 0.;
  Float_t normme[] = {(Float_t)proton_aver, (Float_t)pion_aver, (Float_t)kaon_aver};
  ///////////
  cout << " ==================================" << endl;
  cout << " main calculation starts " <<endl;
  timer.Reset(); timer.Start();
  cout << " ==================================" << endl;
  cout << "  initFunctions  " << endl;
  //
  //
  initFunctions();
  cout << "  calcIntegrals  " << endl;
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
  pi_aver = pion_aver;
  k_aver = kaon_aver;

  cout << "test " << wK_P << "  " << wK_K << "  " << wK_Pi << endl;
  cout << "1 " << W2K_aver << endl;
  cout << "2 " << pr_aver * (wK_P2 - wK_P * wK_P) << endl;
  cout << "3 " << pi_aver * (wK_Pi2 - wK_Pi * wK_Pi) << endl;
  cout << "4 " << k_aver * (wK_K2 - wK_K * wK_K) << endl;
  //
  //
  cout << "  Unfolding      " << endl;
  TMatrixD A2(6, 6);
  A2(0, 0) = wP_P * wP_P;
  A2(0, 1) = wP_Pi * wP_Pi;
  A2(0, 2) = wP_K * wP_K;
  A2(0, 3) = 2. * wP_P * wP_Pi;
  A2(0, 4) = 2. * wP_P * wP_K;
  A2(0, 5) = 2. * wP_Pi * wP_K;

  A2(1, 0) = wPi_P * wPi_P;
  A2(1, 1) = wPi_Pi * wPi_Pi;
  A2(1, 2) = wPi_K * wPi_K;
  A2(1, 3) = 2. * wPi_P * wPi_Pi;
  A2(1, 4) = 2. * wPi_P * wPi_K;
  A2(1, 5) = 2. * wPi_Pi * wPi_K;

  A2(2, 0) = wK_P * wK_P;
  A2(2, 1) = wK_Pi * wK_Pi;
  A2(2, 2) = wK_K * wK_K;
  A2(2, 3) = 2. * wK_P * wK_Pi;
  A2(2, 4) = 2. * wK_P * wK_K;
  A2(2, 5) = 2. * wK_Pi * wK_K;

  A2(3, 0) = wP_P * wPi_P;
  A2(3, 1) = wP_Pi * wPi_Pi;
  A2(3, 2) = wP_K * wPi_K;
  A2(3, 3) = wP_P * wPi_Pi + wP_Pi * wPi_P;
  A2(3, 4) = wP_P * wPi_K + wP_K * wPi_P;
  A2(3, 5) = wP_Pi * wPi_K + wP_K * wPi_Pi;

  A2(4, 0) = wP_P * wK_P;
  A2(4, 1) = wP_Pi * wK_Pi;
  A2(4, 2) = wP_K * wK_K;
  A2(4, 3) = wP_P * wK_Pi + wP_Pi * wK_P;
  A2(4, 4) = wP_P * wK_K + wP_K * wK_P;
  A2(4, 5) = wP_Pi * wK_K + wP_K * wK_Pi;

  A2(5, 0) = wPi_P * wK_P;
  A2(5, 1) = wPi_Pi * wK_Pi;
  A2(5, 2) = wPi_K * wK_K;
  A2(5, 3) = wPi_P * wK_Pi + wPi_Pi * wK_P;
  A2(5, 4) = wPi_P * wK_K + wPi_K * wK_P;
  A2(5, 5) = wPi_Pi * wK_K + wPi_K * wK_Pi;

  proton_aver = WPM_aver;
  kaon_aver = WKM_aver;
  pion_aver = WPiM_aver;

  pr_aver = proton_aver;
  pi_aver = pion_aver;
  k_aver = kaon_aver;

  Double_t B2[6];
  B2[0] = W2P_aver -
          pr_aver * (wP_P2 - wP_P * wP_P) -
          pi_aver * (wP_Pi2 - wP_Pi * wP_Pi) -
          k_aver * (wP_K2 - wP_K * wP_K);

  B2[1] = W2Pi_aver -
          pr_aver * (wPi_P2 - wPi_P * wPi_P) -
          pi_aver * (wPi_Pi2 - wPi_Pi * wPi_Pi) -
          k_aver * (wPi_K2 - wPi_K * wPi_K);

  B2[2] = W2K_aver -
          pr_aver * (wK_P2 - wK_P * wK_P) -
          pi_aver * (wK_Pi2 - wK_Pi * wK_Pi) -
          k_aver * (wK_K2 - wK_K * wK_K);

  B2[3] = WPPi_aver -
          pr_aver * (wPPi_P - wP_P * wPi_P) -
          pi_aver * (wPPi_Pi - wP_Pi * wPi_Pi) -
          k_aver * (wPPi_K - wP_K * wPi_K);

  B2[4] = WPK_aver -
          pr_aver * (wPK_P - wP_P * wK_P) -
          pi_aver * (wPK_Pi - wP_Pi * wK_Pi) -
          k_aver * (wPK_K - wP_K * wK_K);

  B2[5] = WPiK_aver -
          pr_aver * (wPiK_P - wPi_P * wK_P) -
          pi_aver * (wPiK_Pi - wPi_Pi * wK_Pi) -
          k_aver * (wPiK_K - wPi_K * wK_K);

  cout << "A2(0,0) " << A2(0, 0) << "  " << B2[0] << endl;
  A2.Print();

  TMatrixD invA2 = A2.Invert();
  //   invA.Print();
  Double_t recP2_av = 0.;
  Double_t recPi2_av = 0.;
  Double_t recK2_av = 0.;
  Double_t recPPi_av = 0.;
  Double_t recPK_av = 0.;
  Double_t recPiK_av = 0.;

  //recP2_av = recPi2_av = recK2_av  = 0.;
  //recPPi_av = recPK_av = recPiK_av = 0.;
  for (Int_t tt = 0; tt < 6; tt++)
  {
    recP2_av += invA2(0, tt) * B2[tt];
    recPi2_av += invA2(1, tt) * B2[tt];
    recK2_av += invA2(2, tt) * B2[tt];
    //cout<< "kaon test "<<invA2(2,tt)<<"  "<<B2[tt]<<endl;
    recPPi_av += invA2(3, tt) * B2[tt];
    recPK_av += invA2(4, tt) * B2[tt];
    recPiK_av += invA2(5, tt) * B2[tt];
  }

  cout << "print second moments " << endl;
  cout << "pr2  " << recP2_av  << endl;
  cout << "pi2  " << recPi2_av << endl;
  cout << "k2   " << recK2_av  << endl;
  cout << "prpi " << recPPi_av << endl;
  cout << "prk  " << recPK_av  << endl;
  cout << "pik  " << recPiK_av << endl;

  pr2_aver = recP2_av;
  pi2_aver = recPi2_av;
  k2_aver = recK2_av;
  prpi_aver = recPPi_av;
  prk_aver = recPK_av;
  pik_aver = recPiK_av;

  TMatrixD A(10, 10);
  TMatrixD BB(10, 9);

  wmean[0][0] = wP_P;
  wmean[0][1] = wP_Pi;
  wmean[0][2] = wP_K;
  wmean[1][0] = wPi_P;
  wmean[1][1] = wPi_Pi;
  wmean[1][2] = wPi_K;
  wmean[2][0] = wK_P;
  wmean[2][1] = wK_Pi;
  wmean[2][2] = wK_K;

  wmean2[0][0] = wP_P2;
  wmean2[0][1] = wP_Pi2;
  wmean2[0][2] = wP_K2;

  wmean2[1][0] = wPi_P2;
  wmean2[1][1] = wPi_Pi2;
  wmean2[1][2] = wPi_K2;

  wmean2[2][0] = wK_P2;
  wmean2[2][1] = wK_Pi2;
  wmean2[2][2] = wK_K2;

  wmean3[0][0] = wP_P3;
  wmean3[0][1] = wP_Pi3;
  wmean3[0][2] = wP_K3;

  wmean3[1][0] = wPi_P3;
  wmean3[1][1] = wPi_Pi3;
  wmean3[1][2] = wPi_K3;

  wmean3[2][0] = wK_P3;
  wmean3[2][1] = wK_Pi3;
  wmean3[2][2] = wK_K3;

  wmix[0][0] = wPPi_P;
  wmix[0][1] = wPPi_Pi;
  wmix[0][2] = wPPi_K;

  wmix[1][0] = wPK_P;
  wmix[1][1] = wPK_Pi;
  wmix[1][2] = wPK_K;

  wmix[2][0] = wPiK_P;
  wmix[2][1] = wPiK_Pi;
  wmix[2][2] = wPiK_K;

  Float_t wmix2A[3][3];

  wmix2A[0][0] = wP2Pi_P;
  wmix2A[0][1] = wP2Pi_Pi;
  wmix2A[0][2] = wP2Pi_K;

  wmix2A[1][0] = wP2K_P;
  wmix2A[1][1] = wP2K_Pi;
  wmix2A[1][2] = wP2K_K;

  wmix2A[2][0] = wPi2K_P;
  wmix2A[2][1] = wPi2K_Pi;
  wmix2A[2][2] = wPi2K_K;

  Float_t wmix2B[3][3];

  wmix2B[0][0] = wPi2P_P;
  wmix2B[0][1] = wPi2P_Pi;
  wmix2B[0][2] = wPi2P_K;

  wmix2B[1][0] = wK2P_P;
  wmix2B[1][1] = wK2P_Pi;
  wmix2B[1][2] = wK2P_K;

  wmix2B[2][0] = wK2Pi_P;
  wmix2B[2][1] = wK2Pi_Pi;
  wmix2B[2][2] = wK2Pi_K;

  Float_t wmeanPrPiK[3];

  wmeanPrPiK[0] = wPiPrK_Pr;
  wmeanPrPiK[1] = wPiPrK_Pi;
  wmeanPrPiK[2] = wPiPrK_K;

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
    BB(p, nn) = wmeanPrPiK[i] + 2. * wmean[0][i] * wmean[1][i] * wmean[2][i] - wmix[0][i] * wmean[2][i] - wmix[1][i] * wmean[1][i] - wmix[2][i] * wmean[0][i];
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

  Double_t B[10] = {W3P_aver, W3Pi_aver, W3K_aver, WPiPrK_aver, WPr2Pi_aver,
                    WPr2K_aver, WPi2K_aver, WPi2Pr_aver, WK2Pr_aver, WK2Pi_aver};

  for (Int_t i = 0; i < 10; i++)
  {
    B[i] -= (BB(i, 0) * pr2_aver + BB(i, 1) * pi2_aver +
             BB(i, 2) * k2_aver + BB(i, 3) * prpi_aver +
             BB(i, 4) * prk_aver + BB(i, 5) * pik_aver +
             BB(i, 6) * pr_aver + BB(i, 7) * pi_aver + BB(i, 8) * k_aver);
  }

  TString mom3Names[] = {"NP3", "NPi3", "NK3", "Pr2Pi", "Pr2K",
                         "Pi2K", "Pi2Pr", "K2Pr", "K2Pi", "PiPrK"};

  TMatrixD invA = A.Invert();
  Double_t rec3mom[10] = {0};
  cout << "print third moments " << endl;
  for (Int_t i = 0; i < 10; i++)
  {
    for (Int_t tt = 0; tt < 10; tt++)
    {
      rec3mom[i] += invA(i, tt) * B[tt];
    }
    // N^3+3N^2+N --> poisson 3rd moment
    // (A*A+A)*B
    cout << mom3Names[i].Data() << "   :  "  << rec3mom[i] << endl;
  }
  cout << "====================================" << endl;
  cout << " calculation is over " << endl;
  timer.Stop(); timer.Print();
  cout << "====================================" << endl;

  //cout<<"rec3mom 0 "<<rec3mom[0]<<endl;
  Double_t skew = rec3mom[0] - 3. * pr2_aver * proton_aver + 2. * TMath::Power(proton_aver, 3);

  // cout<<"third variance "<<skew/proton_aver<<endl;
  cout << " =========================== " << endl;
  cout << "ratio to poisson 3rd moment  " << (TMath::Power(proton_aver, 3) + 3. * TMath::Power(proton_aver, 2) + proton_aver) / rec3mom[0] << endl;
  cout << "ratio to poisson 2nd moment  " << (TMath::Power(proton_aver, 2) + proton_aver) / pr2_aver << endl;
  cout << " =========================== " << endl;

  skew /= (pr2_aver - proton_aver * proton_aver);
  cout << "skew === " << skew << endl;

  return 0;
}

int main(int argc, char *argv[])
{
  cout << "***********************************************" << endl;
  cout << "***********************************************" << endl;
  cout << "******************** IDENTITY METHOD **********" << endl;
  cout << "***********************************************" << endl;
  cout << "***********************************************" << endl;
  cout << "      " << endl;
  cout << "      " << endl;
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
    TString outPutFileNAme = Form("TIMoments_sub%d_cent_%3.2f_mom_%3.2f_%3.2f_eta_%3.2f_%3.2f.root",fSubsample,fCentInput,fpDownInput,fpUpInput,fEtaDownInput,fEtaUpInput);
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
  //
  //
  TROOT IdentityMethod("IdentityMethod", "compiled identity method");
  ReadFitParamsFromTree(fileNameLineShapes,6);
  identity3Particle();
  return 0;
}
