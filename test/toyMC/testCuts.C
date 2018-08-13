#include "TSystem.h"
#include "TROOT.h"
#include "TCut.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1D.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;
using std::setw;


enum trackCutBit {
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


const Int_t fnCutBins=7;
Int_t fCutArr[fnCutBins] = {kCRows80,  kChi2TPC4, kDCAXY, kVZ, kEventVertexZ, kClusterRequirementITS, kNewITSCut};





Bool_t ApplyTreeSelection(UInt_t cut);
void   testCopyTree();
TString PrintNumInBinary(UInt_t num);




void testCuts(){
    
    /*
    cd /home/marsland/Desktop/schleching/data/trees
    aliroot
    .L /home/marsland/Desktop/schleching/TIdentity/TIdentity/test/testCuts.C+
    
    */
    
    TFile *f = new TFile("AnalysisResults.root");
    TTree *tree = (TTree*)f->Get("tracks");
    Double_t allEntries = tree->GetEntries();
  
    Float_t eta,cent,ptot,dEdx;
    Int_t event,sign;
    UInt_t cutBit;
    Double_t cRows,tpcchi2;
    tree -> SetBranchAddress("eta"   ,&eta);
    tree -> SetBranchAddress("cent"  ,&cent);
    tree -> SetBranchAddress("ptot"  ,&ptot);
    tree -> SetBranchAddress("dEdx"  ,&dEdx);
    tree -> SetBranchAddress("event" ,&event);
    tree -> SetBranchAddress("cutBit" ,&cutBit);
    tree -> SetBranchAddress("cRows" ,&cRows);
    tree -> SetBranchAddress("tpcchi2" ,&tpcchi2);    
    tree -> SetBranchAddress("sign",&sign);
    
    TH1D *fChi2   = new TH1D("fChi2" ,"fChi2",200 ,0, 10 );
    TH1D *fcRows  = new TH1D("fcRows","fcRows",200 ,0, 200 );

    
    for ( int ient = 0; ient < allEntries/2; ient++ ) {
        tree -> GetEntry(ient);
        
        
        // if (!ApplyTreeSelection(cutBit)) continue;
        
        Bool_t select = (((cutBit >> kCRows80) & 1) == 1) && (((cutBit >> kChi2TPC4) & 1) == 1);
        if (!select) continue;
        
        TString cutBinary = PrintNumInBinary(cutBit);
        if (ient%50000==0 && tpcchi2>3) 
            cout << cutBinary << "  " << cutBit << "  " << cRows << "  " << tpcchi2 << endl;
      
        fChi2->Fill(tpcchi2);
        fcRows->Fill(cRows);
        
        
    }
    
    TFile * outfile = new TFile("testCuts.root","recreate");
    fChi2->Write("fChi2"); 
    fcRows->Write("fcRows"); 
    delete outfile;
    
    
        

}

void testCopyTree(){
    
    
    TFile *f = new TFile("AnalysisResults.root");
    TTree *tree = (TTree*)f->Get("tracks");
    tree->GetEntries();
    
    TString mm = Form("((cutBit >> %d) & 1) == 1 && ((cutBit >> %d) & 1) == 1",kCRows80,kChi2TPC4);
    tree->SetAlias("nn",mm);
    Double_t allEntries = tree->GetEntries();
    TFile * tmpFile = new TFile("tmp.root","recreate");
    
    TStopwatch timer; 
    timer.Start();
    TTree * t = tree->CopyTree("nn","",allEntries);
    timer.Stop(); timer.Print();
    
    tree->Draw("cRows");
    t->Draw("cRows","","same");
    
}

Bool_t ApplyTreeSelection(UInt_t cut)
{
    UInt_t arr[fnCutBins];
    for (Int_t i=0;i<fnCutBins;i++){
        arr[i] = ((cut >> fCutArr[i]) & 1);
    }
    
    return (arr[0]&&arr[1]&&arr[2]&&arr[3]&&arr[4]&&arr[5]&&arr[6]);
}

TString PrintNumInBinary(UInt_t num)
{
    TString bin="";
    Int_t numberOfBits = sizeof(UInt_t)*8;
    for (Int_t i=numberOfBits-1; i>=0; i--) {
        Bool_t isBitSet = (num & (1<<i));
        if (isBitSet) {
            bin+="1";
        } else {
            bin+="0";
        }
    }
    return bin;
}