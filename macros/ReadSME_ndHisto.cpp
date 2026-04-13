#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TMath.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TF1.h"
#include "TH1.h"
#include "TTree.h"

using namespace std;

enum EventType { kSameEvent, kMixedEvent };

bool FileExists_bool(const string &filename) {
  ifstream file(filename);
  return file.good();
}

void FileExists_warn(const string &filename) {
  if (!FileExists_bool(filename)) {
    cerr << "Error: File: " << filename << " doesn't exist!" << endl;
  }
}

unique_ptr<TFile> GetROOT(const string Readpath, const string ReadFilename, const string Option) {
  string lowerOption = Option;
  transform(lowerOption.begin(), lowerOption.end(), lowerOption.begin(), [](unsigned char c) {
    return tolower(c);
  });  // 全转小写
  if (!(lowerOption == "read" || lowerOption == "write" || lowerOption == "recreate")) {
    cout << "Invalid Option: " << Option << ", should be one of: read, write, recreate" << endl;
  }
  string rfilename = Readpath + "/" + ReadFilename + ".root";  // 加上root后缀
  FileExists_warn(rfilename);

  return make_unique<TFile>(rfilename.c_str(), Option.c_str());
}

std::string doubleToString(double x, int precision = 2) {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(precision) << x;
  return ss.str();
}

void onefromndHisto(string rpath, string rfilename, string wpath, string wfilename) {
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

  TH1D *h1s = h4_s->Projection(0);
  TH1D *h1sc = (TH1D *)h1s->Clone();
  TH1D *h1m = h4_m->Projection(0);
  TH1D *h1mc = (TH1D *)h1m->Clone();

  // do cf
  TF1 *constant = new TF1("constant", "1", 0, 10);
  int binNorm[2];
  binNorm[0] = h1sc->FindBin(0.5);
  binNorm[1] = h1sc->FindBin(0.8);
  double factorN = h1mc->Integral(binNorm[0], binNorm[1]) / h1sc->Integral(binNorm[0], binNorm[1]);
  h1mc->Divide(constant, factorN);
  h1sc->Divide(h1mc);
  TH1D *h1sc_re = (TH1D *)h1sc->Clone();

  int binRange[2];
  //   binRange[0] = h1sc->FindBin(0.);
  binRange[0] = h1sc->FindBin(0.);
  binRange[1] = h1sc->FindBin(0.25);

  h1sc_re->Reset();  // 清空内容

  for (int i = binRange[0]; i <= binRange[1]; i++) {
    double content = h1sc->GetBinContent(i);
    double error = h1sc->GetBinError(i);
    int newBin = i - binRange[0] + 1;  // 从1开始填
    h1sc_re->SetBinContent(newBin, content);
    h1sc_re->SetBinError(newBin, error);
  }

  // 重新设置x轴范围
  h1sc_re->GetXaxis()->SetRangeUser(0., 0.25);

  h1sc->SetTitle("CF from ndHisto");
  h1sc->GetXaxis()->SetTitle("k* (GeV/c)");
  h1sc->GetYaxis()->SetTitle("C(k*)");
  h1sc_re->SetTitle("reranged CF from ndHisto");
  h1sc_re->GetXaxis()->SetTitle("k* (GeV/c)");
  h1sc_re->GetYaxis()->SetTitle("C(k*)");
  TCanvas *c1 = new TCanvas("c1", "CF from ndHisto", 800, 600);
  c1->SetMargin(0.13, 0.05, 0.12, 0.07);
  c1->SetGrid();
  h1sc->Draw("hist");
  c1->Update();

  wfFemtoep->cd();
  h1sc->Write("CF_from_ndHisto");
  h1sc_re->Write("CF_reranged_from_ndHisto");
  wfFemtoep->Close();
  rfFemtoep->Close();
}

void ndHistoRead(string rpath, string rfilename, string wpath, string wfilename, EventType eventType) {
  // string fullpath = filepath + "/" +filename;
  auto rfFemtoep = GetROOT(rpath, rfilename, "read");
  auto wfFemtoep = GetROOT(wpath, wfilename, "recreate");
  THnSparseF *h4;
  if (eventType == kSameEvent) {
    cout << "Processing Same Event ndHisto..." << endl;
    h4 = (THnSparseF *)rfFemtoep->Get(
        "femto-dream-pair-task-track-track/SameEvent_EP/"
        "relPairkstarmTMultMultPercentileQn");
  } else {
    cout << "Processing Mixed Event ndHisto..." << endl;
    h4 = (THnSparseF *)rfFemtoep->Get(
        "femto-dream-pair-task-track-track/MixedEvent_EP/"
        "relPairkstarmTMultMultPercentileQn");
  }

  // 获取轴
  auto ax_kstar = h4->GetAxis(0);
  auto ax_mT = h4->GetAxis(1);
  auto ax_cent = h4->GetAxis(2);
  auto ax_EP = h4->GetAxis(3);

  // 定义中心度和 mT 区间
  std::vector<std::pair<double, double>> centBins = {{0, 30}, {30, 70}, {70, 100}};
  std::vector<std::pair<double, double>> mTBins = {{0.5, 0.7}, {0.7, 1.0}, {1.0, 1.5}};

  // 定义 EP 对称区间
  struct EPBin {
    double low1, high1;
    double low2, high2;
    const char *name;
  };
  std::vector<EPBin> epBins;
  if (eventType == kSameEvent) {
    std::cout << "Using Same Event EP bins." << std::endl;
    epBins = {{0, TMath::Pi() / 4, 3 * TMath::Pi() / 4, TMath::Pi(), "In_plane"},
              {TMath::Pi() / 4, 3 * TMath::Pi() / 4, -1, -1, "Out_of_plane"}};
  } else {
    std::cout << "Using Mixed Event EP bins." << std::endl;
    epBins = {{0, TMath::Pi(), -1, -1, "Min bias EP"}};
  }

  int colorList[] = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7};

  wfFemtoep->cd();

  // 主循环：每个中心度和 mT 各一张图
  for (auto [centLow, centHigh] : centBins) {
    for (auto [mTLow, mTHigh] : mTBins) {
      // 创建新的画布
      string cname = "cent=" + to_string(static_cast<int>(centLow)) + "-" + to_string(static_cast<int>(centHigh))
                     + " mT=" + doubleToString(mTLow, 1) + "-" + doubleToString(mTHigh, 1);
      // TCanvas *c = new TCanvas(cname, cname, 800,600);
      //  c->SetMargin(0.13, 0.05, 0.12, 0.07); c->SetGrid();

      TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.04);

      std::cout << cname << std::endl;

      // 限制中心度与 mT
      ax_cent->SetRangeUser(centLow, centHigh);
      ax_mT->SetRangeUser(mTLow, mTHigh);

      int iColor = 0;

      for (auto &ep : epBins) {
        // 如果是双区间（in-plane），要合并两部分
        TH1D *h_kstar_total = nullptr;

        if (ep.low2 > 0) {  // 有镜像部分
          ax_EP->SetRangeUser(ep.low1, ep.high1);
          auto h1 = (TH1D *)h4->Projection(0);
          ax_EP->SetRangeUser(ep.low2, ep.high2);
          auto h2 = (TH1D *)h4->Projection(0);
          h1->Add(h2);
          h_kstar_total = h1;
        } else {
          ax_EP->SetRangeUser(ep.low1, ep.high1);
          h_kstar_total = (TH1D *)h4->Projection(0);
        }
        // 改到这
        //  设置样式
        h_kstar_total->SetLineColor(colorList[iColor++ % 4]);
        h_kstar_total->SetLineWidth(2);
        h_kstar_total->SetTitle(cname.c_str());
        h_kstar_total->GetXaxis()->SetTitle("k* (GeV/c)");
        h_kstar_total->GetYaxis()->SetTitle("Counts");

        // 归一化
        // if (h_kstar_total->Integral() > 0)
        //     h_kstar_total->Scale(1.0 / h_kstar_total->Integral());

        string Histname = cname + " " + ep.name;
        if (iColor == 1)
          h_kstar_total->Draw("hist");
        else
          h_kstar_total->Draw("hist same");

        h_kstar_total->Write(Histname.c_str());
        leg->AddEntry(h_kstar_total, ep.name, "l");
      }
      // leg->Write();
      // c->Update();
    }
  }
}

void ReadSME_ndHisto() {
  // string rpath = "/Users/allenzhou/ALICE/scripts/femtoep";
  // string rpath = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  string rpath = "/Users/allenzhou/ALICE/alidata/hyperloopres/femtoep";
  // string rfilename = "AnalysisResults_femtoep_sm";
  // string rfilename = "AnalysisResults_femtoep_rec";
  string rfilename = "AnalysisResults";
  string wpath = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  // string wpath = "/Users/allenzhou/ALICE/scripts/femtoep";
  string wfilename1 = "Femtoep_Same";
  string wfilename2 = "Femtoep_Mixed";
  // ndHistoRead(path, filename);
  // onefromndHisto(rpath, rfilename, wpath, wfilename);
  ndHistoRead(rpath, rfilename, wpath, wfilename1, kSameEvent);
  ndHistoRead(rpath, rfilename, wpath, wfilename2, kMixedEvent);
}