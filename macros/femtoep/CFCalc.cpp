#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "TF1.h"
#include "TH1.h"
#include "TTree.h"
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TMath.h>

using namespace std;

// enum EventType { kSameEvent, kMixedEvent };
// enum Pairphibin { In_Plane, Out_of_plane };

bool FileExists_bool(const string &filename) {
  ifstream file(filename);
  return file.good();
}

void FileExists_warn(const string &filename) {
  if (!FileExists_bool(filename)) {
    cerr << "Error: File: " << filename << " doesn't exist!" << endl;
  }
}

unique_ptr<TFile> GetROOT(const string Readpath, const string ReadFilename,
                          const string Option) {
  string lowerOption = Option;
  transform(lowerOption.begin(), lowerOption.end(), lowerOption.begin(),
            [](unsigned char c) { return tolower(c); }); // 全转小写
  if (!(lowerOption == "read" || lowerOption == "write" ||
        lowerOption == "recreate")) {
    cout << "Invalid Option: " << Option
         << ", should be one of: read, write, recreate" << endl;
  }
  string rfilename = Readpath + "/" + ReadFilename + ".root"; // 加上root后缀
  FileExists_warn(rfilename);
  cout << "Opening file: " << rfilename << " with option: " << Option << endl;
  return make_unique<TFile>(rfilename.c_str(), Option.c_str());
}

std::string doubleToString(double x, int precision = 2) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << x;
  return ss.str();
}

TH1D *GetCFfromSM(TH1D *h_se, TH1D *h_me, double normLow, double normHigh,
                  double kstarMax) {
  TF1 *constant = new TF1("constant", "1", 0, 10);
  cout << "Normalizing SE and ME histograms between k* = " << normLow << " and "
       << normHigh << " GeV/c." << endl;
  TH1D *h_se_c = (TH1D *)h_se->Clone();
  TH1D *h_me_c = (TH1D *)h_me->Clone();
  cout << "Cloned SE and ME histograms." << endl;

  int binNorm[2];
  binNorm[0] = h_se_c->FindBin(normLow);
  binNorm[1] = h_se_c->FindBin(normHigh);
  double factorN = h_me_c->Integral(binNorm[0], binNorm[1]) /
                   h_se_c->Integral(binNorm[0], binNorm[1]);
  h_me_c->Divide(constant, factorN);
  h_se_c->Divide(h_me_c);
  TH1D *h_cf_re = new TH1D();
  int Upperbin = h_se_c->FindBin(kstarMax);

  for (int nBin = 1; nBin <= Upperbin; nBin++) {
    double content = h_se_c->GetBinContent(nBin);
    double error = h_se_c->GetBinError(nBin);
    h_cf_re->SetBinContent(nBin, content);
    h_cf_re->SetBinError(nBin, error);
  }
  return h_cf_re;
}

void onefromndHisto(string rpath, string rfilename, string wpath,
                    string wfilename) {
  // string fullpath = filepath + "/" +filename;
  auto rfFemtoep = GetROOT(rpath, rfilename, "read");
  auto wfFemtoep = GetROOT(wpath, wfilename, "recreate");
  auto h4_s = (THnSparseF *)rfFemtoep->Get(
      "femto-dream-pair-task-track-track/SameEvent_EP/"
      "relPairkstarmTMultMultPercentileQn");
  auto h4_m = (THnSparseF *)rfFemtoep->Get(
      "femto-dream-pair-task-track-track/MixedEvent_EP/"
      "relPairkstarmTMultMultPercentileQn");

  // 获取轴
  auto ax_kstar_s = h4_s->GetAxis(0);
  auto ax_kstar_m = h4_m->GetAxis(0);

  TH1D *h1s = h4_s->Projection(3);
  TH1D *h1sc = (TH1D *)h1s->Clone();
  TH1D *h1m = h4_m->Projection(3);
  TH1D *h1mc = (TH1D *)h1m->Clone();

  // do cf
  TF1 *constant = new TF1("constant", "1", 0, 10);
  // int binNorm[2];
  // binNorm[0] = h1sc->FindBin(0.5);
  // binNorm[1] = h1sc->FindBin(0.8);
  // double factorN = h1mc->Integral(binNorm[0], binNorm[1]) /
  //                  h1sc->Integral(binNorm[0], binNorm[1]);
  // h1mc->Divide(constant, factorN);
  h1mc->Divide(h1sc);
  TH1D *h1sc_re = (TH1D *)h1mc->Clone();

  //   int binRange[2];
  // //   binRange[0] = h1sc->FindBin(0.);
  //   binRange[0] = h1sc->FindBin(0.);
  //   binRange[1] = h1sc->FindBin(0.25);

  // h1sc_re->Reset(); // 清空内容

  // for (int i = binRange[0]; i <= binRange[1];i++) {
  //   double content = h1sc->GetBinContent(i);
  //   double error = h1sc->GetBinError(i);
  //   int newBin = i - binRange[0] + 1; // 从1开始填
  //   h1sc_re->SetBinContent(newBin, content);
  //   h1sc_re->SetBinError(newBin, error);
  // }

  // 重新设置x轴范围
  // h1sc_re->GetXaxis()->SetRangeUser(0., 0.25);

  // h1sc->SetTitle("CF from ndHisto");
  // h1sc->GetXaxis()->SetTitle("k* (GeV/c)");
  // h1sc->GetYaxis()->SetTitle("C(k*)");
  h1sc_re->SetTitle("M/S Ratio");
  h1sc_re->GetXaxis()->SetTitle("Ratio");
  h1sc_re->GetYaxis()->SetTitle("Phi-Psi");
  TCanvas *c1 = new TCanvas("c1", "CF from ndHisto", 800, 600);
  c1->SetMargin(0.13, 0.05, 0.12, 0.07);
  c1->SetGrid();
  h1sc->Draw("hist");
  c1->Update();

  wfFemtoep->cd();
  h1sc_re->Write("M/S ratio");
  wfFemtoep->Close();
  rfFemtoep->Close();
}

void CFCalcWith_Cent_Mt_pairphi(string rpath_se, string rfilename_se,
                                string rpath_me, string rfilename_me,
                                string wpath, string wfilename,
                                std::vector<std::pair<double, double>> centBins,
                                std::vector<std::pair<double, double>> mTBins,
                                std::vector<std::string> pairphiBins) {
  auto rf_se = GetROOT(rpath_se, rfilename_se, "read");
  auto rf_me = GetROOT(rpath_me, rfilename_me, "read");
  auto wfFemtoep = GetROOT(wpath, wfilename, "recreate");

  // 重归一化的区间
  int binNorm[2];
  // 用于计算的fake histo
  double factorN;
  // 最后要显示的区间
  int binRange[2];
  //   binRange[0] = h1sc->FindBin(0.);

  TF1 *constant = new TF1("constant", "1", 0, 10);

  cout << "Starting CF calculation with cent, mT, pairphi bins..." << endl;

  for (auto [centLow, centHigh] : centBins) {
    for (auto [mTLow, mTHigh] : mTBins) {
      string hname = "cent=" + to_string(static_cast<int>(centLow)) + "-" +
                     to_string(static_cast<int>(centHigh)) +
                     " mT=" + doubleToString(mTLow, 1) + "-" +
                     doubleToString(mTHigh, 1);
      string hname_me = hname + " " + "Min bias EP";
      TH1D *h1_me = rf_me->Get<TH1D>(hname_me.c_str());
      if (!h1_me) {
        cout << "ERROR: ME histo not found: " << hname_me << endl;
        continue;
      }
      TH1D *h1mc = (TH1D *)(h1_me->Clone());
      for (auto phibin : pairphiBins) {
        string hname_se = hname + " " + phibin;
        auto h1_se = (TH1D *)rf_se->Get(hname_se.c_str());
        if (!h1_se) {
          cout << "ERROR: ME histo not found: " << hname_me << endl;
          continue;
        }
        TH1D *h1sc = (TH1D *)h1_se->Clone();
        cout << "Processing: " << hname_se << endl;
        // TH1D *h1cf_re = GetCFfromSM(h1_se, h1_me, 0.5, 0.8, 0.25);

        binNorm[0] = h1sc->FindBin(0.5);
        binNorm[1] = h1sc->FindBin(0.8);
        factorN = h1mc->Integral(binNorm[0], binNorm[1]) /
                  h1sc->Integral(binNorm[0], binNorm[1]);
        h1mc->Divide(constant, factorN);
        h1sc->Divide(h1mc);

        binRange[0] = h1sc->FindBin(0.);
        binRange[1] = h1sc->FindBin(0.25);

        int nBins = binRange[1] - binRange[0] + 1;
        TH1D *h1cf_re = new TH1D("", "", nBins, 0, 0.25);
        // TH1D *h1cf_re = new TH1D();

        for (int i = binRange[0]; i <= binRange[1]; i++) {
          double content = h1sc->GetBinContent(i);
          double error = h1sc->GetBinError(i);
          int newBin = i - binRange[0] + 1; // 从1开始填
          cout << "Bin " << newBin << " content: " << content
               << ", error: " << error << endl;
          h1cf_re->SetBinContent(newBin, content);
          h1cf_re->SetBinError(newBin, error);
        }
        string wHistoname = "CF " + hname + " " + phibin;
        h1cf_re->SetTitle("reranged CF from ndHisto");
        h1cf_re->GetXaxis()->SetTitle("k* (GeV/c)");
        h1cf_re->GetYaxis()->SetTitle("C(k*)");
        wfFemtoep->WriteObject(h1cf_re, wHistoname.c_str());
      }
    }
  }
}

void CFCalc() {
  string rpath1 = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  // string rpath2 = "/Users/allenzhou/ALICE/scripts/femtoep";
  string rpath2 = rpath1;
  string rfilename1 = "Femtoep_Same";
  string rfilename2 = "Femtoep_Mixed";
  // string wpath = "/Users/allenzhou/ALICE/scripts/femtoep";
  string wpath = rpath1;
  string wfilename1 = "Check EP Femto";
  // string wfilename2 = "Check EP mixing 2";

  std::vector<std::pair<double, double>> centBins = {
      {0, 30}, {30, 70}, {70, 100}};
  std::vector<std::pair<double, double>> mTBins = {
      {0.5, 0.7}, {0.7, 1.0}, {1.0, 1.5}};
  std::vector<std::string> pairphiBins = {"In_plane", "Out_of_plane"};

  CFCalcWith_Cent_Mt_pairphi(rpath1, rfilename1, rpath1, rfilename2, wpath,
                             wfilename1, centBins, mTBins, pairphiBins);
  // ndHistoRead(path, filename);
  // onefromndHisto(rpath1, rfilename1, wpath, wfilename1);
  // onefromndHisto(rpath2, rfilename2, wpath, wfilename2);
}