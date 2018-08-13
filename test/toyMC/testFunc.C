#include "TF1.h"
#include "TClonesArray.h"
#include "TFile.h"
// Double_t myfunction(Double_t *x, Double_t *par)
// {
//   Double_t f = par[0]*exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
//   return f;
// }

// void testFunc(){
//
//   TClonesArray arr("TF1",50);
//   arr.SetOwner(kTRUE);
//
//   // TF1 f1("f1","sin(x)",0,10);
//   // TF1 f2("f2","cos(x)",0,10);
//   // TF1 fsum("f1","[&](double *x, double *p){ return p[0]*f1(x) + p[1]*f2(x); }",0,10,2);
//   // fsum.SetParameters(1,2);
//   // fsum.Draw();
//
//   // TF1 *func = new TF1("func",myfunction,0,50,3);
//   // func->SetParameters(10,10,1);
//
//   TF1 * func = new TF1("func",[&](double*x, double *p){ return p[0]*exp(-0.5*TMath::Power((x[0]-p[1])/p[2],2)); }, 0, 50, 3);
//   func->SetParameters(10,13,2);
//
//   TF1 *func1 = new TF1("func1","[0]*exp(-0.5*TMath::Power((x-[1])/[2],2))",0,50);
//   func1->SetParameters(10,10,2);
//
//   func->Draw();
//   func1->Draw("same");
//   // TF1 fsum("f1","[&](double *x, double *p){ return par[0]*exp(-0.5*TMath::Power((x-par[1])/par[2],2)); }",0,10,2);
//
//
//   // new (arr[objcounter++]) TF1(*fParticles[j][i]);
//   // arr[0] = (TF1*)func;
//   new(arr[0]) TF1(*func1);
//   new(arr[1]) TF1(*func);
//   // new(arr[1]) TF1(fsum);
//   TFile *f = new TFile("aaa.root","recreate");
//   arr . Write("arr",TObject::kSingleKey);
//   //arr . Write();
//   // arr . Clear("C");
//   f->Close();
//   delete f;
//
// }



double fun1(double *xVec, double *Param)
{
  double x, fVal;
  x=xVec[0];
  fVal=Param[0]*x*x + Param[1]*x + Param[0];
  return fVal;
}

double fun2(double *xVec, double *Param)
{
  double x, fVal;
  x=xVec[0];
  fVal=Param[0]*x*x + Param[1]*x + Param[0];
  return fVal;
}

double fun3(double *xVec, double *Param)
{
  return fun1(xVec, Param) + fun2(xVec, &Param[3]);
}

void balic()
{

  TClonesArray arr("TF1",50);
  arr.SetOwner(kTRUE);

  double Params[6]={1, 1, 1, 1, 1, 1};

  TF1 *f1 = new TF1("fun1", fun1, -10, 10, 3);
  f1->SetParameters(Params);
  TF1 *f2 = new TF1("fun2", fun2, -10, 10, 3);
  f2->SetParameters(Params);
  printf(" f1(2)=%f\n", f1->Eval(2));
  printf(" f2(2)=%f\n", f2->Eval(2));

  // the following statement  produces an error "ERROR 4 : Empty String"
  TF1 *f3 = new TF1("fun3", fun3, -10, 10, 6);
  TF1 *fClone = new TF1(); f3->Copy(*fClone);
  f3->SetParameters(Params);
  f3->Draw();
  printf(" f3(2)=%f\n", f3->Eval(2));

  new(arr[0]) TF1(*fClone);
  TFile *f = new TFile("bbb.root","recreate");
  arr . Write("arr",TObject::kSingleKey);
  f->Close();
  delete f;






}
