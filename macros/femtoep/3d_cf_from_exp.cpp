#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "TF1.h"
#include "TF3.h"
#include "TH1.h"
#include "TH3.h"
#include "TTree.h"
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <THnSparse.h>
#include <TKey.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveText.h>

using namespace std;

enum EventType { kSameEvent, kMixedEvent };

struct EPBin {
  double low1, high1;
  double low2, high2;
  string name;
};

constexpr double kProjection1DWindow = 0.04;

bool FileExists_bool(const string &filename) {
  ifstream file(filename);
  return file.good();
}

void FileExists_warn(const string &filename) {
  if (!FileExists_bool(filename)) {
    cerr << "Error: File: " << filename << " doesn't exist!" << endl;
  }
}

unique_ptr<TFile> GetROOT_unique(const string Readpath,
                                 const string ReadFilename,
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

  return make_unique<TFile>(rfilename.c_str(), Option.c_str());
}

TFile *GetROOT(const string Readpath, const string ReadFilename,
               const string Option) {
  string lowerOption = Option;
  transform(lowerOption.begin(), lowerOption.end(), lowerOption.begin(),
            [](unsigned char c) { return tolower(c); }); // 全转小写
  if (!(lowerOption == "read" || lowerOption == "write" ||
        lowerOption == "recreate" || lowerOption == "update")) {
    cout << "Invalid Option: " << Option
         << ", should be one of: read, write, recreate or update" << endl;
  }
  string rfilename = Readpath + "/" + ReadFilename + ".root"; // 加上root后缀
  FileExists_warn(rfilename);

  return new TFile(rfilename.c_str(), Option.c_str());
}

bool InitializeROOTFile(const string &path, const string &filename) {
  auto wf = GetROOT(path, filename, "recreate");
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot initialize ROOT file " << path << "/" << filename
         << ".root" << endl;
    delete wf;
    return false;
  }
  wf->Close();
  delete wf;
  return true;
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

TH1D *onefromndHisto(unique_ptr<TFile> rfile, string hpath, bool isReranged,
                     double rerangeMax) {
  auto h4 = (THnSparseF *)rfile->Get(hpath.c_str());

  // 获取轴
  auto ax_kstar = h4->GetAxis(0);

  TH1D *h1 = h4->Projection(0);
  TH1D *h1c = (TH1D *)h1->Clone();

  h1c->SetTitle("CF from ndHisto");
  h1c->GetXaxis()->SetTitle("k* (GeV/c)");
  h1c->GetYaxis()->SetTitle("C(k*)");

  if (!isReranged) {
    return h1c;
  }

  // do cf
  //   TF1 *constant = new TF1("constant", "1", 0, 10);
  //   int binNorm[2];
  //   binNorm[0] = h1sc->FindBin(0.5);
  //   binNorm[1] = h1sc->FindBin(0.8);
  //   double factorN = h1mc->Integral(binNorm[0], binNorm[1]) /
  //                    h1sc->Integral(binNorm[0], binNorm[1]);
  //   h1mc->Divide(constant, factorN);
  //   h1sc->Divide(h1mc);
  //   TH1D *h1sc_re = (TH1D *)h1sc->Clone();

  int binRange[2];
  //   binRange[0] = h1sc->FindBin(0.);
  binRange[0] = h1c->FindBin(0.);
  binRange[1] = h1c->FindBin(rerangeMax);

  TH1D *h1c_re = (TH1D *)h1c->Clone();

  h1c_re->Reset(); // 清空内容

  for (int i = binRange[0]; i <= binRange[1]; i++) {
    double content = h1c->GetBinContent(i);
    double error = h1c->GetBinError(i);
    int newBin = i - binRange[0] + 1; // 从1开始填
    h1c_re->SetBinContent(newBin, content);
    h1c_re->SetBinError(newBin, error);
  }

  // 重新设置x轴范围
  h1c_re->GetXaxis()->SetRangeUser(0., 0.25);

  h1c_re->SetTitle("reranged CF from ndHisto");
  h1c_re->GetXaxis()->SetTitle("k* (GeV/c)");
  h1c_re->GetYaxis()->SetTitle("C(k*)");
  return h1c_re;
  //   TCanvas *c1 = new TCanvas("c1", "CF from ndHisto", 800, 600);
  //   c1->SetMargin(0.13, 0.05, 0.12, 0.07);
  //   c1->SetGrid();
  //   h1sc->Draw("hist");
  //   c1->Update();

  //   wfFemtoep->cd();
  //   h1sc->Write("CF_from_ndHisto");
  //   h1sc_re->Write("CF_reranged_from_ndHisto");
  //   wfFemtoep->Close();
  //   rfFemtoep->Close();
}

TH1D *ndHistoRead(TFile *rf, string histoname, int centLow, int centHigh,
                  double mTLow, double mTHigh, EPBin epbin) {
  // string fullpath = filepath + "/" +filename;
  THnSparseF *h4 = (THnSparseF *)rf->Get(histoname.c_str());
  //                              "relPairkstarmTMultMultPercentileQn");;
  // if (eventType == kSameEvent) {
  //   cout << "Processing Same Event ndHisto..." << endl;
  //   h4 = (THnSparseF
  //   *)rf->Get("femto-dream-pair-task-track-track/SameEvent_EP/"
  //                              "relPairkstarmTMultMultPercentileQn");
  // } else {
  //   cout << "Processing Mixed Event ndHisto..." << endl;
  //   h4 =
  //       (THnSparseF
  //       *)rf->Get("femto-dream-pair-task-track-track/MixedEvent_EP/"
  //                             "relPairkstarmTMultMultPercentileQn");
  // }

  // 获取轴
  auto ax_kstar = h4->GetAxis(0);
  auto ax_mT = h4->GetAxis(1);
  auto ax_cent = h4->GetAxis(2);
  auto ax_EP = h4->GetAxis(3);

  //   std::vector<EPBin> epBins;
  //   if (eventType == kSameEvent) {
  //     std::cout << "Using Same Event EP bins." << std::endl;
  //     epBins = {
  //         {0, TMath::Pi() / 4, 3 * TMath::Pi() / 4, TMath::Pi(), "In_plane"},
  //         {TMath::Pi() / 4, 3 * TMath::Pi() / 4, -1, -1, "Out_of_plane"}};
  //   } else {
  //     std::cout << "Using Mixed Event EP bins." << std::endl;
  //     epBins = {{0, TMath::Pi(), -1, -1, "Min bias EP"}};
  //   }

  int colorList[] = {kRed + 1, kBlue + 1, kGreen + 2, kOrange + 7};

  // 创建新的画布
  string cname = "cent=" + to_string(static_cast<int>(centLow)) + "-" +
                 to_string(static_cast<int>(centHigh)) +
                 "_mT=" + doubleToString(mTLow, 1) + "-" +
                 doubleToString(mTHigh, 1);
  // TCanvas *c = new TCanvas(cname, cname, 800,600);
  //  c->SetMargin(0.13, 0.05, 0.12, 0.07); c->SetGrid();

  TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);

  // std::cout << cname << std::endl;

  // 限制中心度与 mT
  ax_cent->SetRangeUser(centLow, centHigh);
  ax_mT->SetRangeUser(mTLow, mTHigh);

  int iColor = 0;

  // 如果是双区间（in-plane），要合并两部分
  TH1D *h_kstar_total = nullptr;

  if (epbin.low2 > 0) { // 有镜像部分
    ax_EP->SetRangeUser(epbin.low1, epbin.high1);
    auto h1 = (TH1D *)h4->Projection(0);
    ax_EP->SetRangeUser(epbin.low2, epbin.high2);
    auto h2 = (TH1D *)h4->Projection(0);
    h1->Add(h2);
    h_kstar_total = h1;
  } else {
    ax_EP->SetRangeUser(epbin.low1, epbin.high1);
    h_kstar_total = (TH1D *)h4->Projection(0);
  }
  // 改到这
  //  设置样式
  h_kstar_total->SetLineColor(colorList[iColor++ % 4]);
  h_kstar_total->SetLineWidth(2);
  h_kstar_total->SetTitle(cname.c_str());
  h_kstar_total->GetXaxis()->SetTitle("k* (GeV/c)");
  h_kstar_total->GetYaxis()->SetTitle("Counts");
  string Histname = cname + "_" + epbin.name;
  h_kstar_total->SetName(Histname.c_str());

  return h_kstar_total;

  // 归一化
  // if (h_kstar_total->Integral() > 0)
  //     h_kstar_total->Scale(1.0 / h_kstar_total->Integral());

  //   if (iColor == 1)
  //     h_kstar_total->Draw("hist");
  //   else
  //     h_kstar_total->Draw("hist same");

  // leg->AddEntry(h_kstar_total, ep.name, "l");
}

TH1D *calc_cf_from_sme_rerange(TH1D *h_se, TH1D *h_me,
                               std::pair<double, double> normrange,
                               std::pair<double, double> kstarRange) {
  TF1 *constant = new TF1("constant", "1", 0, 10);

  h_se->Sumw2();
  h_me->Sumw2();

  TH1D *h_se_c = (TH1D *)h_se->Clone();
  TH1D *h_me_c = (TH1D *)h_me->Clone();

  int binNorm[2];
  binNorm[0] = h_se_c->FindBin(normrange.first);
  binNorm[1] = h_se_c->FindBin(normrange.second);
  double factorN = h_me_c->Integral(binNorm[0], binNorm[1]) /
                   h_se_c->Integral(binNorm[0], binNorm[1]);
  h_me_c->Divide(constant, factorN);
  h_se_c->Divide(h_me_c);

  int binRange[2];
  binRange[0] = h_se_c->FindBin(kstarRange.first);
  binRange[1] = h_se_c->FindBin(kstarRange.second);

  int nBins = binRange[1] - binRange[0] + 1;
  TH1D *h_cf_re = new TH1D("", "", nBins, kstarRange.first, kstarRange.second);

  for (int i = binRange[0]; i <= binRange[1]; i++) {
    double content = h_se_c->GetBinContent(i);
    double error = h_se_c->GetBinError(i);
    int newBin = i - binRange[0] + 1; // 从1开始填
    h_cf_re->SetBinContent(newBin, content);
    h_cf_re->SetBinError(newBin, error);
  }
  return h_cf_re;
}

TH1D *BuildProjectionXWithinWindow(TH3D *hist, const string &name, double qMax);
TH1D *BuildProjectionYWithinWindow(TH3D *hist, const string &name, double qMax);
TH1D *BuildProjectionZWithinWindow(TH3D *hist, const string &name, double qMax);

double IntegralVisibleRange(TH3D *hist, bool useWidth = false) {
  if (!hist) {
    return 0.0;
  }
  return hist->Integral(1, hist->GetNbinsX(), 1, hist->GetNbinsY(), 1,
                        hist->GetNbinsZ(), useWidth ? "width" : "");
}

void Write1DProjections(TH3D *hCF, TDirectory *dir, const string &baseName,
                        const string &yTitle = "C(q)", bool useWindow = true) {
  TH1D *hProjX = nullptr;
  TH1D *hProjY = nullptr;
  TH1D *hProjZ = nullptr;

  if (useWindow) {
    hProjX = BuildProjectionXWithinWindow(hCF, baseName + "_ProjX",
                                          kProjection1DWindow);
    hProjY = BuildProjectionYWithinWindow(hCF, baseName + "_ProjY",
                                          kProjection1DWindow);
    hProjZ = BuildProjectionZWithinWindow(hCF, baseName + "_ProjZ",
                                          kProjection1DWindow);
  } else {
    hProjX = hCF->ProjectionX((baseName + "_ProjX").c_str(), 1,
                              hCF->GetNbinsY(), 1, hCF->GetNbinsZ());
    hProjY = hCF->ProjectionY((baseName + "_ProjY").c_str(), 1,
                              hCF->GetNbinsX(), 1, hCF->GetNbinsZ());
    hProjZ = hCF->ProjectionZ((baseName + "_ProjZ").c_str(), 1,
                              hCF->GetNbinsX(), 1, hCF->GetNbinsY());
  }

  string hProjXName = baseName + "_ProjX";
  hProjX->SetName(hProjXName.c_str());
  hProjX->SetTitle(hProjXName.c_str());
  hProjX->GetXaxis()->SetTitle("q_{out} (GeV/c)");
  hProjX->GetYaxis()->SetTitle(yTitle.c_str());
  dir->WriteObject(hProjX, hProjXName.c_str());

  string hProjYName = baseName + "_ProjY";
  hProjY->SetName(hProjYName.c_str());
  hProjY->SetTitle(hProjYName.c_str());
  hProjY->GetXaxis()->SetTitle("q_{side} (GeV/c)");
  hProjY->GetYaxis()->SetTitle(yTitle.c_str());
  dir->WriteObject(hProjY, hProjYName.c_str());

  string hProjZName = baseName + "_ProjZ";
  hProjZ->SetName(hProjZName.c_str());
  hProjZ->SetTitle(hProjZName.c_str());
  hProjZ->GetXaxis()->SetTitle("q_{long} (GeV/c)");
  hProjZ->GetYaxis()->SetTitle(yTitle.c_str());
  dir->WriteObject(hProjZ, hProjZName.c_str());

  delete hProjX;
  delete hProjY;
  delete hProjZ;
}

void CFCalc3D(string rpath, string rfilename, string taskname,
              string subtaskname_se, string subtaskname_me, string wpath,
              string wfilename, std::vector<std::pair<double, double>> centBins,
              std::vector<std::pair<double, double>> mTBins) {
  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!InitializeROOTFile(wpath, wfilename)) {
    if (rf) {
      rf->Close();
      delete rf;
    }
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  string hpath_se =
      taskname + "/" + subtaskname_se + "/relPair3dRmTMultPercentileQnPairphi";
  string hpath_me =
      taskname + "/" + subtaskname_me + "/relPair3dRmTMultPercentileQnPairphi";

  auto hSE_sparse_origin = (THnSparseF *)rf->Get(hpath_se.c_str());
  auto hME_sparse_origin = (THnSparseF *)rf->Get(hpath_me.c_str());

  if (!hSE_sparse_origin || !hME_sparse_origin) {
    cout << "ERROR: Cannot find THnSparse!" << endl;
    if (rf) {
      rf->Close();
      delete rf;
    }
    TH1::AddDirectory(oldAddDirectory);
    return;
  }
  auto hSE_sparse = (THnSparseF *)hSE_sparse_origin->Clone();
  auto hSE_sparse_allphi = (THnSparseF *)hSE_sparse_origin->Clone();
  auto hME_sparse = (THnSparseF *)hME_sparse_origin->Clone();
  hSE_sparse->Sumw2();
  hSE_sparse_allphi->Sumw2();
  hME_sparse->Sumw2();

  auto ax_phi = hSE_sparse->GetAxis(6);
  int nPhiBins = ax_phi->GetNbins();

  cout << "Starting 3D CF calculation..." << endl;

  for (auto [centLow, centHigh] : centBins) {
    for (auto [mTLow, mTHigh] : mTBins) {

      string baseName =
          "cent=" + to_string((int)centLow) + "-" + to_string((int)centHigh) +
          "_mT=" + doubleToString(mTLow, 2) + "-" + doubleToString(mTHigh, 2);

      cout << "Processing " << baseName << endl;

      // 设置 cut（只设置一次）
      hSE_sparse->GetAxis(4)->SetRangeUser(centLow, centHigh);
      hSE_sparse->GetAxis(3)->SetRangeUser(mTLow, mTHigh);

      hSE_sparse_allphi->GetAxis(4)->SetRangeUser(centLow, centHigh);
      hSE_sparse_allphi->GetAxis(3)->SetRangeUser(mTLow, mTHigh);

      hME_sparse->GetAxis(4)->SetRangeUser(centLow, centHigh);
      hME_sparse->GetAxis(3)->SetRangeUser(mTLow, mTHigh);

      auto hME_norm = (TH3D *)hME_sparse->Projection(0, 1, 2);
      hME_norm->SetDirectory(nullptr);
      const double intME = IntegralVisibleRange(hME_norm, true);
      if (intME == 0.0) {
        cout << "WARNING: zero mixed-event integral for " << baseName << endl;
        delete hME_norm;
        continue;
      }
      hME_norm->Scale(1.0 / intME);

      auto writeStoredSlice = [&](TH3D *hSE_norm, const string &hname) {
        string hSEName = hname + "_SE_norm3d";
        string hMEName = hname + "_ME_norm3d";

        hSE_norm->SetName(hSEName.c_str());
        hSE_norm->SetTitle(hSEName.c_str());

        auto hMEStored = (TH3D *)hME_norm->Clone(hMEName.c_str());
        hMEStored->SetDirectory(nullptr);
        hMEStored->SetTitle(hMEName.c_str());

        auto hCF = (TH3D *)hSE_norm->Clone(hname.c_str());
        hCF->SetDirectory(nullptr);
        hCF->Divide(hMEStored);

        hCF->SetName(hname.c_str());
        hCF->SetTitle(hname.c_str());
        hCF->GetXaxis()->SetTitle("q_{out} (GeV/c)");
        hCF->GetYaxis()->SetTitle("q_{side} (GeV/c)");
        hCF->GetZaxis()->SetTitle("q_{long} (GeV/c)");

        auto wf = GetROOT(wpath, wfilename, "update");
        if (!wf || wf->IsZombie()) {
          cout << "ERROR: cannot update output ROOT file " << wpath << "/"
               << wfilename << ".root" << endl;
          delete wf;
          delete hCF;
          delete hMEStored;
          return;
        }

        auto dir = wf->GetDirectory(hname.c_str());
        if (!dir) {
          dir = wf->mkdir(hname.c_str());
        }
        dir->WriteObject(hSE_norm, hSEName.c_str());
        Write1DProjections(hSE_norm, dir, hSEName, "Normalized density", true);
        dir->WriteObject(hMEStored, hMEName.c_str());
        Write1DProjections(hMEStored, dir, hMEName, "Normalized density", true);
        dir->WriteObject(hCF, hname.c_str());
        Write1DProjections(hCF, dir, hname, "C(q)", true);
        wf->Close();
        delete wf;
        delete hCF;
        delete hMEStored;
      };

      auto hSE_all = (TH3D *)hSE_sparse_allphi->Projection(0, 1, 2);
      hSE_all->SetDirectory(nullptr);
      const double intSEAll = IntegralVisibleRange(hSE_all, true);
      if (intSEAll == 0.0) {
        cout << "WARNING: zero same-event integral for " << baseName
             << " phi=all" << endl;
      } else {
        hSE_all->Scale(1.0 / intSEAll);
        writeStoredSlice(hSE_all, baseName + "_phi=all_CF3D");
      }
      delete hSE_all;

      // 遍历 pair-phi bin，并把 (pi/2, pi) 映射到 (-pi/2, 0)
      for (int i = 1; i <= nPhiBins; ++i) {

        double phi = ax_phi->GetBinCenter(i);

        // -------- SE --------
        hSE_sparse->GetAxis(6)->SetRange(i, i);
        auto hSE_a = (TH3D *)hSE_sparse->Projection(0, 1, 2);
        hSE_a->SetDirectory(nullptr);

        // -------- 归一化（关键：分别归一）--------
        double intSE = IntegralVisibleRange(hSE_a, true);

        if (intSE == 0.0) {
          cout << "WARNING: zero integral, skip bin " << i << endl;
          delete hSE_a;
          continue;
        }

        hSE_a->Scale(1.0 / intSE);

        // -------- 命名并写出 --------
        double phi_mapped = phi;
        if (phi > TMath::Pi() / 2.) {
          phi_mapped -= TMath::Pi();
        }
        string hname =
            baseName + "_phi=" + doubleToString(phi_mapped, 2) + "_CF3D";
        writeStoredSlice(hSE_a, hname);

        delete hSE_a;
      }

      delete hME_norm;
    }
  }

  rf->Close();
  delete hSE_sparse;
  delete hSE_sparse_allphi;
  delete hME_sparse;
  delete rf;
  cout << "3D CF results have been written to " << wpath << "/" << wfilename
       << ".root" << endl;
  TH1::AddDirectory(oldAddDirectory);
}

struct Levy3DFitResult {
  string fitModel = "diag";
  string histName;
  string groupKey;
  int centLow = -1;
  int centHigh = -1;
  double mTLow = 0.;
  double mTHigh = 0.;
  double phi = 0.;
  bool isPhiIntegrated = false;
  double norm = 0.;
  double normErr = 0.;
  double lambda = 0.;
  double lambdaErr = 0.;
  double rout2 = 0.;
  double rout2Err = 0.;
  double rside2 = 0.;
  double rside2Err = 0.;
  double rlong2 = 0.;
  double rlong2Err = 0.;
  double routside2 = 0.;
  double routside2Err = 0.;
  double routlong2 = 0.;
  double routlong2Err = 0.;
  double rsidelong2 = 0.;
  double rsidelong2Err = 0.;
  double alpha = 0.;
  double alphaErr = 0.;
  double chi2 = 0.;
  int ndf = 0;
  int status = -1;
  bool hasOffDiagonal = false;
  bool usesCoulomb = false;
  bool usesCoreHaloLambda = true;
};

struct Levy3DFitOptions {
  bool useCoulomb = false;
  bool useCoreHaloLambda = true;
  double fitQMax = 0.15;
};

constexpr double kHbarC = 0.1973269804;
constexpr double kPiPiLikeSignBohrRadiusFm = 387.5;

bool EndsWith(const string &value, const string &suffix) {
  if (value.size() < suffix.size()) {
    return false;
  }
  return value.compare(value.size() - suffix.size(), suffix.size(), suffix) ==
         0;
}

bool ParseCF3DHistogramName(const string &histName, Levy3DFitResult &result) {
  int centLow = 0;
  int centHigh = 0;
  double mTLow = 0.;
  double mTHigh = 0.;
  double phi = 0.;

  int matched = sscanf(histName.c_str(), "cent=%d-%d_mT=%lf-%lf_phi=%lf_CF3D",
                       &centLow, &centHigh, &mTLow, &mTHigh, &phi);
  if (matched == 5) {
    result.histName = histName;
    result.centLow = centLow;
    result.centHigh = centHigh;
    result.mTLow = mTLow;
    result.mTHigh = mTHigh;
    result.phi = phi;
    result.isPhiIntegrated = false;
    result.groupKey = "cent=" + to_string(centLow) + "-" + to_string(centHigh) +
                      "_mT=" + doubleToString(mTLow, 2) + "-" +
                      doubleToString(mTHigh, 2);
    return true;
  }

  int matchedAll =
      sscanf(histName.c_str(), "cent=%d-%d_mT=%lf-%lf_phi=all_CF3D", &centLow,
             &centHigh, &mTLow, &mTHigh);
  if (matchedAll != 4) {
    return false;
  }

  result.histName = histName;
  result.centLow = centLow;
  result.centHigh = centHigh;
  result.mTLow = mTLow;
  result.mTHigh = mTHigh;
  result.phi = -1.0;
  result.isPhiIntegrated = true;
  result.groupKey = "cent=" + to_string(centLow) + "-" + to_string(centHigh) +
                    "_mT=" + doubleToString(mTLow, 2) + "-" +
                    doubleToString(mTHigh, 2);
  return true;
}

TH3D *LoadStoredCFHistogram(TFile *rf, const string &objectName) {
  if (!rf) {
    return nullptr;
  }

  if (auto hCF = dynamic_cast<TH3D *>(rf->Get(objectName.c_str()))) {
    auto hClone = (TH3D *)hCF->Clone((objectName + "_data3d").c_str());
    hClone->SetDirectory(nullptr);
    return hClone;
  }

  auto dir = dynamic_cast<TDirectory *>(rf->Get(objectName.c_str()));
  if (!dir) {
    return nullptr;
  }

  auto hCF = dynamic_cast<TH3D *>(dir->Get(objectName.c_str()));
  if (!hCF) {
    return nullptr;
  }

  auto hClone = (TH3D *)hCF->Clone((objectName + "_data3d").c_str());
  hClone->SetDirectory(nullptr);
  return hClone;
}

string BuildFitOptionTag(const Levy3DFitOptions &options) {
  return string(options.useCoreHaloLambda ? "lam" : "nolam") + "_" +
         (options.useCoulomb ? "coul" : "nocoul");
}

double ComputeLikeSignPiPiGamowFactor(double qOut, double qSide, double qLong) {
  const double qMagnitude = std::sqrt(qOut * qOut + qSide * qSide +
                                      qLong * qLong); // Approximate q_inv.
  const double kStarFm = 0.5 * qMagnitude / kHbarC;
  if (kStarFm <= 1e-12) {
    return 0.0;
  }

  const double eta = 1.0 / (kStarFm * kPiPiLikeSignBohrRadiusFm);
  const double twoPiEta = 2.0 * TMath::Pi() * eta;
  if (twoPiEta > 700.0) {
    return 0.0;
  }

  const double denominator = std::exp(twoPiEta) - 1.0;
  if (denominator <= 0.0) {
    return 0.0;
  }

  const double gamow = twoPiEta / denominator;
  return std::max(0.0, std::min(gamow, 1.0));
}

double ComputeBowlerSinyukovLikeSignPiPiValue(
    double norm, double lambda, double levyExponent, bool useCoulomb,
    bool useCoreHaloLambda, double qOut, double qSide, double qLong) {
  const double lambdaEff = useCoreHaloLambda ? lambda : 1.0;
  const double coulombFactor =
      useCoulomb ? ComputeLikeSignPiPiGamowFactor(qOut, qSide, qLong) : 1.0;
  const double quantumStatTerm = std::exp(-levyExponent);
  return norm * ((1.0 - lambdaEff) +
                 lambdaEff * coulombFactor * (1.0 + quantumStatTerm));
}

double Levy3DModel(double *x, double *par) {
  const bool useCoulomb = par[6] > 0.5;
  const bool useCoreHaloLambda = par[7] > 0.5;
  const double qOut = x[0];
  const double qSide = x[1];
  const double qLong = x[2];
  const double qOut2 = x[0] * x[0];
  const double qSide2 = x[1] * x[1];
  const double qLong2 = x[2] * x[2];

  const double argument =
      (par[2] * qOut2 + par[3] * qSide2 + par[4] * qLong2) / (kHbarC * kHbarC);
  const double levyExponent = std::pow(std::max(argument, 0.0), par[5] / 2.0);

  return ComputeBowlerSinyukovLikeSignPiPiValue(par[0], par[1], levyExponent,
                                                useCoulomb, useCoreHaloLambda,
                                                qOut, qSide, qLong);
}

double Levy3DFullModel(double *x, double *par) {
  const bool useCoulomb = par[9] > 0.5;
  const bool useCoreHaloLambda = par[10] > 0.5;
  const double qOut = x[0];
  const double qSide = x[1];
  const double qLong = x[2];

  const double argument =
      (par[2] * qOut * qOut + par[3] * qSide * qSide + par[4] * qLong * qLong +
       2.0 * par[5] * qOut * qSide + 2.0 * par[6] * qOut * qLong +
       2.0 * par[7] * qSide * qLong) /
      (kHbarC * kHbarC);
  const double levyExponent = std::pow(std::max(argument, 0.0), par[8] / 2.0);

  return ComputeBowlerSinyukovLikeSignPiPiValue(par[0], par[1], levyExponent,
                                                useCoulomb, useCoreHaloLambda,
                                                qOut, qSide, qLong);
}

TF3 *BuildLevyFitFunction(const string &funcName,
                          const Levy3DFitOptions &fitOptions) {
  auto fitFunc =
      new TF3(funcName.c_str(), Levy3DModel, -fitOptions.fitQMax,
              fitOptions.fitQMax, -fitOptions.fitQMax, fitOptions.fitQMax,
              -fitOptions.fitQMax, fitOptions.fitQMax, 8);
  fitFunc->SetParName(0, "Norm");
  fitFunc->SetParName(1, "lambda");
  fitFunc->SetParName(2, "Rout2");
  fitFunc->SetParName(3, "Rside2");
  fitFunc->SetParName(4, "Rlong2");
  fitFunc->SetParName(5, "alpha");
  fitFunc->SetParName(6, "UseCoulomb");
  fitFunc->SetParName(7, "UseCoreHaloLambda");

  fitFunc->SetParameters(1.0, 0.5, 25.0, 25.0, 25.0, 1.5, 0.0, 1.0);
  fitFunc->SetParLimits(0, 0.5, 1.5);
  fitFunc->SetParLimits(1, 0.0, 1.0);
  fitFunc->SetParLimits(2, 0.01, 400.0);
  fitFunc->SetParLimits(3, 0.01, 400.0);
  fitFunc->SetParLimits(4, 0.01, 400.0);
  fitFunc->SetParLimits(5, 0.5, 2.0);
  fitFunc->FixParameter(6, fitOptions.useCoulomb ? 1.0 : 0.0);
  fitFunc->FixParameter(7, fitOptions.useCoreHaloLambda ? 1.0 : 0.0);
  if (!fitOptions.useCoreHaloLambda) {
    fitFunc->FixParameter(1, 1.0);
  }
  fitFunc->SetNpx(60);
  fitFunc->SetNpy(60);
  fitFunc->SetNpz(60);
  return fitFunc;
}

TF3 *BuildFullLevyFitFunction(const string &funcName,
                              const Levy3DFitOptions &fitOptions) {
  auto fitFunc =
      new TF3(funcName.c_str(), Levy3DFullModel, -fitOptions.fitQMax,
              fitOptions.fitQMax, -fitOptions.fitQMax, fitOptions.fitQMax,
              -fitOptions.fitQMax, fitOptions.fitQMax, 11);
  fitFunc->SetParName(0, "Norm");
  fitFunc->SetParName(1, "lambda");
  fitFunc->SetParName(2, "Rout2");
  fitFunc->SetParName(3, "Rside2");
  fitFunc->SetParName(4, "Rlong2");
  fitFunc->SetParName(5, "Routside2");
  fitFunc->SetParName(6, "Routlong2");
  fitFunc->SetParName(7, "Rsidelong2");
  fitFunc->SetParName(8, "alpha");
  fitFunc->SetParName(9, "UseCoulomb");
  fitFunc->SetParName(10, "UseCoreHaloLambda");

  fitFunc->SetParameters(1.0, 0.5, 25.0, 25.0, 25.0, 0.0, 0.0, 0.0, 1.5, 0.0,
                         1.0);
  fitFunc->SetParLimits(0, 0.5, 1.5);
  fitFunc->SetParLimits(1, 0.0, 1.0);
  fitFunc->SetParLimits(2, 0.01, 400.0);
  fitFunc->SetParLimits(3, 0.01, 400.0);
  fitFunc->SetParLimits(4, 0.01, 400.0);
  fitFunc->SetParLimits(5, -200.0, 200.0);
  fitFunc->SetParLimits(6, -200.0, 200.0);
  fitFunc->SetParLimits(7, -200.0, 200.0);
  fitFunc->SetParLimits(8, 0.5, 2.0);
  fitFunc->FixParameter(9, fitOptions.useCoulomb ? 1.0 : 0.0);
  fitFunc->FixParameter(10, fitOptions.useCoreHaloLambda ? 1.0 : 0.0);
  if (!fitOptions.useCoreHaloLambda) {
    fitFunc->FixParameter(1, 1.0);
  }
  fitFunc->SetNpx(60);
  fitFunc->SetNpy(60);
  fitFunc->SetNpz(60);
  return fitFunc;
}

TH3D *BuildFittedHistogramLike(const TH3D *hData, TF3 *fitFunc,
                               const string &histName) {
  auto hFit = (TH3D *)hData->Clone(histName.c_str());
  hFit->Reset("ICES");
  hFit->SetTitle(histName.c_str());

  for (int ix = 1; ix <= hFit->GetNbinsX(); ++ix) {
    const double x = hFit->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= hFit->GetNbinsY(); ++iy) {
      const double y = hFit->GetYaxis()->GetBinCenter(iy);
      for (int iz = 1; iz <= hFit->GetNbinsZ(); ++iz) {
        const double z = hFit->GetZaxis()->GetBinCenter(iz);
        const double value = fitFunc->Eval(x, y, z);
        hFit->SetBinContent(ix, iy, iz, value);
        hFit->SetBinError(ix, iy, iz, 0.);
      }
    }
  }

  return hFit;
}

std::pair<int, int> GetAxisRangeForWindow(TAxis *axis, double qMax) {
  const int firstBin = axis->FindBin(-qMax + 1e-9);
  const int lastBin = axis->FindBin(qMax - 1e-9);
  return {std::max(firstBin, 1), std::min(lastBin, axis->GetNbins())};
}

TH1D *BuildProjectionXWithinWindow(TH3D *hist, const string &name,
                                   double qMax) {
  auto [yMin, yMax] = GetAxisRangeForWindow(hist->GetYaxis(), qMax);
  auto [zMin, zMax] = GetAxisRangeForWindow(hist->GetZaxis(), qMax);
  auto hProj = hist->ProjectionX(name.c_str(), yMin, yMax, zMin, zMax);
  const int nWindowBins =
      std::max(0, yMax - yMin + 1) * std::max(0, zMax - zMin + 1);
  if (nWindowBins > 0) {
    hProj->Scale(1.0 / static_cast<double>(nWindowBins));
  }
  return hProj;
}

TH1D *BuildProjectionYWithinWindow(TH3D *hist, const string &name,
                                   double qMax) {
  auto [xMin, xMax] = GetAxisRangeForWindow(hist->GetXaxis(), qMax);
  auto [zMin, zMax] = GetAxisRangeForWindow(hist->GetZaxis(), qMax);
  auto hProj = hist->ProjectionY(name.c_str(), xMin, xMax, zMin, zMax);
  const int nWindowBins =
      std::max(0, xMax - xMin + 1) * std::max(0, zMax - zMin + 1);
  if (nWindowBins > 0) {
    hProj->Scale(1.0 / static_cast<double>(nWindowBins));
  }
  return hProj;
}

TH1D *BuildProjectionZWithinWindow(TH3D *hist, const string &name,
                                   double qMax) {
  auto [xMin, xMax] = GetAxisRangeForWindow(hist->GetXaxis(), qMax);
  auto [yMin, yMax] = GetAxisRangeForWindow(hist->GetYaxis(), qMax);
  auto hProj = hist->ProjectionZ(name.c_str(), xMin, xMax, yMin, yMax);
  const int nWindowBins =
      std::max(0, xMax - xMin + 1) * std::max(0, yMax - yMin + 1);
  if (nWindowBins > 0) {
    hProj->Scale(1.0 / static_cast<double>(nWindowBins));
  }
  return hProj;
}

void StyleProjectionHistogram(TH1D *hist, int color, int markerStyle,
                              const string &xTitle) {
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(0.9);
  hist->SetLineWidth(2);
  hist->GetXaxis()->SetTitle(xTitle.c_str());
  hist->GetYaxis()->SetTitle("C(q)");
}

string FormatParameterLine(const string &label, double value, double error,
                           const string &unit = "") {
  std::stringstream ss;
  ss << std::fixed << std::setprecision(3) << label << " = " << value << " #pm "
     << error;
  if (!unit.empty()) {
    ss << " " << unit;
  }
  return ss.str();
}

string BuildFitModeTitle(const Levy3DFitResult &fitResult) {
  return fitResult.hasOffDiagonal ? "Full Levy fit" : "Diagonal Levy fit";
}

string BuildFitSwitchLine(const Levy3DFitResult &fitResult) {
  return string("Coulomb(#pi^{#pm}#pi^{#pm}): ") +
         (fitResult.usesCoulomb ? "on" : "off") +
         ", core-halo: " + (fitResult.usesCoreHaloLambda ? "on" : "off");
}

TPaveText *BuildFitParameterBox(const Levy3DFitResult &fitResult, double x1,
                                double y1, double x2, double y2) {
  auto box = new TPaveText(x1, y1, x2, y2, "NDC");
  box->SetFillColor(0);
  box->SetFillStyle(1001);
  box->SetBorderSize(1);
  box->SetTextAlign(12);
  box->SetTextFont(42);
  box->SetTextSize(fitResult.hasOffDiagonal ? 0.024 : 0.028);

  box->AddText(BuildFitModeTitle(fitResult).c_str());
  box->AddText(BuildFitSwitchLine(fitResult).c_str());
  box->AddText(
      FormatParameterLine("N", fitResult.norm, fitResult.normErr).c_str());
  if (fitResult.usesCoreHaloLambda) {
    box->AddText(
        FormatParameterLine("#lambda", fitResult.lambda, fitResult.lambdaErr)
            .c_str());
  } else {
    box->AddText("#lambda fixed = 1.000");
  }
  box->AddText(
      FormatParameterLine("#alpha", fitResult.alpha, fitResult.alphaErr)
          .c_str());
  box->AddText(FormatParameterLine("R_{out}^{2}", fitResult.rout2,
                                   fitResult.rout2Err, "fm^{2}")
                   .c_str());
  box->AddText(FormatParameterLine("R_{side}^{2}", fitResult.rside2,
                                   fitResult.rside2Err, "fm^{2}")
                   .c_str());
  box->AddText(FormatParameterLine("R_{long}^{2}", fitResult.rlong2,
                                   fitResult.rlong2Err, "fm^{2}")
                   .c_str());

  if (fitResult.hasOffDiagonal) {
    box->AddText(FormatParameterLine("R_{outside}^{2}", fitResult.routside2,
                                     fitResult.routside2Err, "fm^{2}")
                     .c_str());
    box->AddText(FormatParameterLine("R_{outlong}^{2}", fitResult.routlong2,
                                     fitResult.routlong2Err, "fm^{2}")
                     .c_str());
    box->AddText(FormatParameterLine("R_{sidelong}^{2}", fitResult.rsidelong2,
                                     fitResult.rsidelong2Err, "fm^{2}")
                     .c_str());
  }

  std::stringstream chi2Line;
  chi2Line << std::fixed << std::setprecision(2)
           << "#chi^{2}/NDF = " << fitResult.chi2 << "/" << fitResult.ndf;
  box->AddText(chi2Line.str().c_str());

  return box;
}

TCanvas *BuildProjectionCanvas(const string &canvasName, TH1D *hData,
                               TH1D *hFit, const string &xTitle,
                               const Levy3DFitResult &fitResult) {
  StyleProjectionHistogram(hData, kBlack, 20, xTitle);
  StyleProjectionHistogram(hFit, kRed + 1, 24, xTitle);
  hFit->SetMarkerSize(0.);

  const double maxValue = std::max(hData->GetMaximum(), hFit->GetMaximum());
  hData->SetMaximum(maxValue * 1.15);

  auto canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
  canvas->SetMargin(0.13, 0.05, 0.12, 0.07);
  canvas->SetGrid();
  hData->Draw("E1");
  hFit->Draw("HIST SAME");

  auto legend = new TLegend(0.62, 0.72, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->AddEntry(hData, "Data projection", "lep");
  legend->AddEntry(hFit, "Levy fit projection", "l");
  legend->Draw();

  auto parameterBox = BuildFitParameterBox(
      fitResult, 0.16, fitResult.hasOffDiagonal ? 0.40 : 0.50, 0.58, 0.88);
  parameterBox->Draw();

  canvas->Update();
  return canvas;
}

TCanvas *Build3DComparisonCanvas(const string &canvasName, TH3D *hData,
                                 TH3D *hFit, const Levy3DFitResult &fitResult) {
  auto canvas = new TCanvas(canvasName.c_str(), canvasName.c_str(), 1400, 600);
  canvas->Divide(2, 1);
  canvas->cd(1);
  gPad->SetTheta(24);
  gPad->SetPhi(32);
  hData->Draw("BOX2Z");
  canvas->cd(2);
  gPad->SetTheta(24);
  gPad->SetPhi(32);
  hFit->Draw("BOX2Z");
  auto parameterBox = BuildFitParameterBox(
      fitResult, 0.12, fitResult.hasOffDiagonal ? 0.28 : 0.40, 0.58, 0.88);
  parameterBox->Draw();
  canvas->Update();
  return canvas;
}

void WriteFitResultsSummary(const string &txtPath,
                            const vector<Levy3DFitResult> &results);
void WriteR2Graphs(TFile *wf, const vector<Levy3DFitResult> &results);

TH3D *BuildCFHistogramFromSparse(THnSparseF *hSE_sparse, THnSparseF *hME_sparse,
                                 int phiBin, int phiBinSym,
                                 const string &histName) {
  (void)phiBinSym;
  hSE_sparse->GetAxis(6)->SetRange(phiBin, phiBin);
  auto hSE_a = (TH3D *)hSE_sparse->Projection(0, 1, 2);
  hSE_a->SetDirectory(nullptr);

  auto hME_a = (TH3D *)hME_sparse->Projection(0, 1, 2);
  hME_a->SetDirectory(nullptr);

  const double intSE = IntegralVisibleRange(hSE_a, true);
  const double intME = IntegralVisibleRange(hME_a, true);
  if (intSE == 0.0 || intME == 0.0) {
    delete hSE_a;
    delete hME_a;
    return nullptr;
  }

  hSE_a->Scale(1.0 / intSE);
  hME_a->Scale(1.0 / intME);

  auto hCF = (TH3D *)hSE_a->Clone(histName.c_str());
  hCF->SetDirectory(nullptr);
  hCF->SetTitle(histName.c_str());
  hCF->Divide(hME_a);
  hCF->GetXaxis()->SetTitle("q_{out} (GeV/c)");
  hCF->GetYaxis()->SetTitle("q_{side} (GeV/c)");
  hCF->GetZaxis()->SetTitle("q_{long} (GeV/c)");

  delete hSE_a;
  delete hME_a;
  return hCF;
}

bool FitAndWriteSingleCFHistogram(TH3D *hCF, Levy3DFitResult &fitResult,
                                  const string &wpath, const string &wfilename,
                                  bool useFullModel,
                                  const Levy3DFitOptions &fitOptions) {
  auto fillResultFromFunction = [&](TF3 *fitFunc) {
    fitResult.fitModel = useFullModel ? "full" : "diag";
    fitResult.hasOffDiagonal = useFullModel;
    fitResult.usesCoulomb = fitOptions.useCoulomb;
    fitResult.usesCoreHaloLambda = fitOptions.useCoreHaloLambda;
    fitResult.norm = fitFunc->GetParameter(0);
    fitResult.normErr = fitFunc->GetParError(0);
    fitResult.lambda =
        fitOptions.useCoreHaloLambda ? fitFunc->GetParameter(1) : 1.0;
    fitResult.lambdaErr =
        fitOptions.useCoreHaloLambda ? fitFunc->GetParError(1) : 0.0;
    fitResult.rout2 = fitFunc->GetParameter(2);
    fitResult.rout2Err = fitFunc->GetParError(2);
    fitResult.rside2 = fitFunc->GetParameter(3);
    fitResult.rside2Err = fitFunc->GetParError(3);
    fitResult.rlong2 = fitFunc->GetParameter(4);
    fitResult.rlong2Err = fitFunc->GetParError(4);
    if (useFullModel) {
      fitResult.routside2 = fitFunc->GetParameter(5);
      fitResult.routside2Err = fitFunc->GetParError(5);
      fitResult.routlong2 = fitFunc->GetParameter(6);
      fitResult.routlong2Err = fitFunc->GetParError(6);
      fitResult.rsidelong2 = fitFunc->GetParameter(7);
      fitResult.rsidelong2Err = fitFunc->GetParError(7);
      fitResult.alpha = fitFunc->GetParameter(8);
      fitResult.alphaErr = fitFunc->GetParError(8);
    } else {
      fitResult.alpha = fitFunc->GetParameter(5);
      fitResult.alphaErr = fitFunc->GetParError(5);
    }
    fitResult.chi2 = fitFunc->GetChisquare();
    fitResult.ndf = fitFunc->GetNDF();
  };

  const string objectName = fitResult.histName;

  TF3 *fitFunc =
      useFullModel
          ? BuildFullLevyFitFunction(objectName + "_levy3d_full_fit",
                                     fitOptions)
          : BuildLevyFitFunction(objectName + "_levy3d_fit", fitOptions);
  auto fitStatus = hCF->Fit(fitFunc, "RSMNQ0");
  fillResultFromFunction(fitFunc);
  fitResult.status = static_cast<int>(fitStatus);

  const string fitHistName =
      objectName + (useFullModel ? "_full_fit3d" : "_fit3d");
  auto hFit3D = BuildFittedHistogramLike(hCF, fitFunc, fitHistName);

  auto hProjXData = BuildProjectionXWithinWindow(
      hCF, objectName + "_data_ProjX", kProjection1DWindow);
  auto hProjYData = BuildProjectionYWithinWindow(
      hCF, objectName + "_data_ProjY", kProjection1DWindow);
  auto hProjZData = BuildProjectionZWithinWindow(
      hCF, objectName + "_data_ProjZ", kProjection1DWindow);

  auto hProjXFit = BuildProjectionXWithinWindow(
      hFit3D, objectName + (useFullModel ? "_full_fit_ProjX" : "_fit_ProjX"),
      kProjection1DWindow);
  auto hProjYFit = BuildProjectionYWithinWindow(
      hFit3D, objectName + (useFullModel ? "_full_fit_ProjY" : "_fit_ProjY"),
      kProjection1DWindow);
  auto hProjZFit = BuildProjectionZWithinWindow(
      hFit3D, objectName + (useFullModel ? "_full_fit_ProjZ" : "_fit_ProjZ"),
      kProjection1DWindow);

  auto cProjX = BuildProjectionCanvas(
      objectName + (useFullModel ? "_canvas_full_ProjX" : "_canvas_ProjX"),
      hProjXData, hProjXFit, "q_{out} (GeV/c)", fitResult);
  auto cProjY = BuildProjectionCanvas(
      objectName + (useFullModel ? "_canvas_full_ProjY" : "_canvas_ProjY"),
      hProjYData, hProjYFit, "q_{side} (GeV/c)", fitResult);
  auto cProjZ = BuildProjectionCanvas(
      objectName + (useFullModel ? "_canvas_full_ProjZ" : "_canvas_ProjZ"),
      hProjZData, hProjZFit, "q_{long} (GeV/c)", fitResult);
  auto c3D = Build3DComparisonCanvas(
      objectName + (useFullModel ? "_canvas_full_3D" : "_canvas_3D"), hCF,
      hFit3D, fitResult);

  auto wf = GetROOT(wpath, wfilename, "update");
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot write output ROOT file " << wpath << "/" << wfilename
         << ".root" << endl;
    delete wf;
    delete cProjX;
    delete cProjY;
    delete cProjZ;
    delete c3D;
    delete hProjXData;
    delete hProjYData;
    delete hProjZData;
    delete hProjXFit;
    delete hProjYFit;
    delete hProjZFit;
    delete hFit3D;
    delete fitFunc;
    return false;
  }

  auto dir = wf->GetDirectory(objectName.c_str());
  if (!dir) {
    dir = wf->mkdir(objectName.c_str());
  }
  dir->cd();

  hCF->Write();
  fitFunc->Write();
  hFit3D->Write();
  hProjXData->Write();
  hProjYData->Write();
  hProjZData->Write();
  hProjXFit->Write();
  hProjYFit->Write();
  hProjZFit->Write();
  cProjX->Write();
  cProjY->Write();
  cProjZ->Write();
  c3D->Write();

  wf->Close();
  delete wf;
  delete cProjX;
  delete cProjY;
  delete cProjZ;
  delete c3D;
  delete hProjXData;
  delete hProjYData;
  delete hProjZData;
  delete hProjXFit;
  delete hProjYFit;
  delete hProjZFit;
  delete hFit3D;
  delete fitFunc;
  return true;
}

bool MatchSelectedBin(const Levy3DFitResult &fitResult,
                      const vector<pair<double, double>> &centBins,
                      const vector<pair<double, double>> &mTBins) {
  const bool centMatched = std::any_of(
      centBins.begin(), centBins.end(), [&](const pair<double, double> &bin) {
        return fitResult.centLow == static_cast<int>(bin.first) &&
               fitResult.centHigh == static_cast<int>(bin.second);
      });
  if (!centMatched) {
    return false;
  }

  return std::any_of(mTBins.begin(), mTBins.end(),
                     [&](const pair<double, double> &bin) {
                       return std::abs(fitResult.mTLow - bin.first) < 1e-6 &&
                              std::abs(fitResult.mTHigh - bin.second) < 1e-6;
                     });
}

void FitCF3DWithSelectedBins(
    const string &rpath, const string &rfilename, const string &wpath,
    const string &wfilename, const string &txtFilename,
    std::vector<std::pair<double, double>> centBins,
    std::vector<std::pair<double, double>> mTBins, bool useFullModel,
    const Levy3DFitOptions &fitOptions = Levy3DFitOptions()) {
  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!rf || rf->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file " << rpath << "/" << rfilename
         << ".root" << endl;
    delete rf;
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  if (!InitializeROOTFile(wpath, wfilename)) {
    rf->Close();
    delete rf;
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  vector<Levy3DFitResult> fitResults;
  TIter nextKey(rf->GetListOfKeys());
  TKey *key = nullptr;

  while ((key = (TKey *)nextKey())) {
    string objectName = key->GetName();
    if (!EndsWith(objectName, "_CF3D")) {
      continue;
    }

    Levy3DFitResult fitResult;
    if (!ParseCF3DHistogramName(objectName, fitResult)) {
      cout << "WARNING: cannot parse histogram name: " << objectName << endl;
      continue;
    }

    if (!MatchSelectedBin(fitResult, centBins, mTBins)) {
      cout << "Skipping histogram " << objectName
           << " due to unmatched centrality or mT bin." << endl;
      continue;
    }

    cout << "Fitting selected histogram " << objectName << endl;

    auto hCF = LoadStoredCFHistogram(rf, objectName);
    if (!hCF) {
      cout << "WARNING: cannot read histogram " << objectName << endl;
      continue;
    }

    if (FitAndWriteSingleCFHistogram(hCF, fitResult, wpath, wfilename,
                                     useFullModel, fitOptions)) {
      fitResults.push_back(fitResult);
    }

    delete hCF;
  }

  auto wf = GetROOT(wpath, wfilename, "update");
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot update output ROOT file " << wpath << "/"
         << wfilename << ".root for summary graphs" << endl;
    delete wf;
  } else {
    WriteR2Graphs(wf, fitResults);
    wf->Close();
    delete wf;
  }

  WriteFitResultsSummary(wpath + "/" + txtFilename + ".txt", fitResults);

  rf->Close();
  delete rf;
  cout << (useFullModel ? "Full 3D Levy fit" : "Diagonal 3D Levy fit")
       << " results have been written to " << wpath << "/" << wfilename
       << ".root and " << wpath << "/" << txtFilename << ".txt" << endl;
  TH1::AddDirectory(oldAddDirectory);
}

void WriteFitResultsSummary(const string &txtPath,
                            const vector<Levy3DFitResult> &results) {
  ofstream out(txtPath);
  out << "# fitModel usesCoulomb usesCoreHaloLambda histName centLow centHigh "
         "mTLow mTHigh phi isPhiIntegrated norm normErr lambda lambdaErr "
         "Rout2 Rout2Err Rside2 Rside2Err Rlong2 Rlong2Err Routside2 "
         "Routside2Err "
         "Routlong2 Routlong2Err Rsidelong2 Rsidelong2Err alpha alphaErr "
         "chi2 ndf status\n";
  out << std::fixed << std::setprecision(6);

  for (const auto &result : results) {
    out << result.fitModel << " " << (result.usesCoulomb ? 1 : 0) << " "
        << (result.usesCoreHaloLambda ? 1 : 0) << " " << result.histName << " "
        << result.centLow << " " << result.centHigh << " " << result.mTLow
        << " " << result.mTHigh << " " << result.phi << " "
        << (result.isPhiIntegrated ? 1 : 0) << " " << result.norm << " "
        << result.normErr << " " << result.lambda << " " << result.lambdaErr
        << " " << result.rout2 << " " << result.rout2Err << " " << result.rside2
        << " " << result.rside2Err << " " << result.rlong2 << " "
        << result.rlong2Err << " " << result.routside2 << " "
        << result.routside2Err << " " << result.routlong2 << " "
        << result.routlong2Err << " " << result.rsidelong2 << " "
        << result.rsidelong2Err << " " << result.alpha << " " << result.alphaErr
        << " " << result.chi2 << " " << result.ndf << " " << result.status
        << "\n";
  }
}

void WriteR2Graphs(TFile *wf, const vector<Levy3DFitResult> &results) {
  map<string, vector<Levy3DFitResult>> groupedResults;
  for (const auto &result : results) {
    if (result.isPhiIntegrated) {
      continue;
    }
    groupedResults[result.groupKey].push_back(result);
  }

  auto graphDir = wf->mkdir("R2_vs_phi");
  graphDir->cd();

  for (auto &[groupKey, groupResults] : groupedResults) {
    if (groupResults.empty()) {
      continue;
    }
    sort(groupResults.begin(), groupResults.end(),
         [](const Levy3DFitResult &lhs, const Levy3DFitResult &rhs) {
           return lhs.phi < rhs.phi;
         });

    const int nPoints = static_cast<int>(groupResults.size());
    const bool hasOffDiagonal = std::any_of(
        groupResults.begin(), groupResults.end(),
        [](const Levy3DFitResult &result) { return result.hasOffDiagonal; });
    const bool usesCoreHaloLambda =
        std::any_of(groupResults.begin(), groupResults.end(),
                    [](const Levy3DFitResult &result) {
                      return result.usesCoreHaloLambda;
                    });
    auto gRout2 = new TGraphErrors(nPoints);
    auto gRside2 = new TGraphErrors(nPoints);
    auto gRlong2 = new TGraphErrors(nPoints);
    TGraphErrors *gRoutside2 =
        hasOffDiagonal ? new TGraphErrors(nPoints) : nullptr;
    TGraphErrors *gRoutlong2 =
        hasOffDiagonal ? new TGraphErrors(nPoints) : nullptr;
    TGraphErrors *gRsidelong2 =
        hasOffDiagonal ? new TGraphErrors(nPoints) : nullptr;
    auto gAlpha = new TGraphErrors(nPoints);
    TGraphErrors *gLambda =
        usesCoreHaloLambda ? new TGraphErrors(nPoints) : nullptr;

    for (int i = 0; i < nPoints; ++i) {
      const auto &result = groupResults[i];
      gRout2->SetPoint(i, result.phi, result.rout2);
      gRout2->SetPointError(i, 0., result.rout2Err);

      gRside2->SetPoint(i, result.phi, result.rside2);
      gRside2->SetPointError(i, 0., result.rside2Err);

      gRlong2->SetPoint(i, result.phi, result.rlong2);
      gRlong2->SetPointError(i, 0., result.rlong2Err);

      if (hasOffDiagonal) {
        gRoutside2->SetPoint(i, result.phi, result.routside2);
        gRoutside2->SetPointError(i, 0., result.routside2Err);

        gRoutlong2->SetPoint(i, result.phi, result.routlong2);
        gRoutlong2->SetPointError(i, 0., result.routlong2Err);

        gRsidelong2->SetPoint(i, result.phi, result.rsidelong2);
        gRsidelong2->SetPointError(i, 0., result.rsidelong2Err);
      }

      gAlpha->SetPoint(i, result.phi, result.alpha);
      gAlpha->SetPointError(i, 0., result.alphaErr);

      if (usesCoreHaloLambda) {
        gLambda->SetPoint(i, result.phi, result.lambda);
        gLambda->SetPointError(i, 0., result.lambdaErr);
      }
    }

    auto fitCosRout2 =
        new TF1((groupKey + "_Rout2_phi_fit").c_str(), "[0]+2.0*[1]*cos(2.0*x)",
                -TMath::Pi() / 2., TMath::Pi() / 2.);
    fitCosRout2->SetParNames("Rout2_0", "Rout2_2");
    fitCosRout2->SetParameters(groupResults.front().rout2, 0.0);

    auto fitCosRside2 =
        new TF1((groupKey + "_Rside2_phi_fit").c_str(),
                "[0]+2.0*[1]*cos(2.0*x)", -TMath::Pi() / 2., TMath::Pi() / 2.);
    fitCosRside2->SetParNames("Rside2_0", "Rside2_2");
    fitCosRside2->SetParameters(groupResults.front().rside2, 0.0);

    auto fitCosRlong2 =
        new TF1((groupKey + "_Rlong2_phi_fit").c_str(),
                "[0]+2.0*[1]*cos(2.0*x)", -TMath::Pi() / 2., TMath::Pi() / 2.);
    fitCosRlong2->SetParNames("Rlong2_0", "Rlong2_2");
    fitCosRlong2->SetParameters(groupResults.front().rlong2, 0.0);

    TF1 *fitSinRoutside2 = nullptr;
    TF1 *fitCosRoutlong2 = nullptr;
    TF1 *fitSinRsidelong2 = nullptr;
    if (hasOffDiagonal) {
      fitSinRoutside2 = new TF1((groupKey + "_Routside2_phi_fit").c_str(),
                                "[0]+2.0*[1]*sin(2.0*x)", -TMath::Pi() / 2.,
                                TMath::Pi() / 2.);
      fitSinRoutside2->SetParNames("Routside2_0", "Routside2_2");
      fitSinRoutside2->SetParameters(groupResults.front().routside2, 0.0);

      fitCosRoutlong2 = new TF1((groupKey + "_Routlong2_phi_fit").c_str(),
                                "[0]+2.0*[1]*cos(2.0*x)", -TMath::Pi() / 2.,
                                TMath::Pi() / 2.);
      fitCosRoutlong2->SetParNames("Routlong2_0", "Routlong2_2");
      fitCosRoutlong2->SetParameters(groupResults.front().routlong2, 0.0);

      fitSinRsidelong2 = new TF1((groupKey + "_Rsidelong2_phi_fit").c_str(),
                                 "[0]+2.0*[1]*sin(2.0*x)", -TMath::Pi() / 2.,
                                 TMath::Pi() / 2.);
      fitSinRsidelong2->SetParNames("Rsidelong2_0", "Rsidelong2_2");
      fitSinRsidelong2->SetParameters(groupResults.front().rsidelong2, 0.0);
    }

    auto fitConstAlpha = new TF1((groupKey + "_alpha_phi_fit").c_str(), "[0]",
                                 -TMath::Pi() / 2., TMath::Pi() / 2.);
    fitConstAlpha->SetParName(0, "alpha0");
    fitConstAlpha->SetParameter(0, groupResults.front().alpha);

    TF1 *fitConstLambda = nullptr;
    if (usesCoreHaloLambda) {
      fitConstLambda = new TF1((groupKey + "_lambda_phi_fit").c_str(), "[0]",
                               -TMath::Pi() / 2., TMath::Pi() / 2.);
      fitConstLambda->SetParName(0, "lambda0");
      fitConstLambda->SetParameter(0, groupResults.front().lambda);
    }

    string routName = groupKey + "_Rout2_vs_phi";
    gRout2->SetName(routName.c_str());
    gRout2->SetTitle((routName + ";#phi_{pair}-#Psi_{EP} (rad);R_{out}^{2} "
                                 "(fm^{2})")
                         .c_str());
    gRout2->SetMarkerStyle(20);
    gRout2->SetMarkerColor(kRed + 1);
    gRout2->SetLineColor(kRed + 1);
    gRout2->Fit(fitCosRout2, "QN");
    gRout2->Write();
    fitCosRout2->Write();

    string rsideName = groupKey + "_Rside2_vs_phi";
    gRside2->SetName(rsideName.c_str());
    gRside2->SetTitle((rsideName + ";#phi_{pair}-#Psi_{EP} (rad);R_{side}^{2} "
                                   "(fm^{2})")
                          .c_str());
    gRside2->SetMarkerStyle(21);
    gRside2->SetMarkerColor(kBlue + 1);
    gRside2->SetLineColor(kBlue + 1);
    gRside2->Fit(fitCosRside2, "QN");
    gRside2->Write();
    fitCosRside2->Write();

    string rlongName = groupKey + "_Rlong2_vs_phi";
    gRlong2->SetName(rlongName.c_str());
    gRlong2->SetTitle((rlongName + ";#phi_{pair}-#Psi_{EP} (rad);R_{long}^{2} "
                                   "(fm^{2})")
                          .c_str());
    gRlong2->SetMarkerStyle(22);
    gRlong2->SetMarkerColor(kGreen + 2);
    gRlong2->SetLineColor(kGreen + 2);
    gRlong2->Fit(fitCosRlong2, "QN");
    gRlong2->Write();
    fitCosRlong2->Write();

    if (hasOffDiagonal) {
      string routsideName = groupKey + "_Routside2_vs_phi";
      gRoutside2->SetName(routsideName.c_str());
      gRoutside2->SetTitle(
          (routsideName +
           ";#phi_{pair}-#Psi_{EP} (rad);R_{outside}^{2} (fm^{2})")
              .c_str());
      gRoutside2->SetMarkerStyle(34);
      gRoutside2->SetMarkerColor(kCyan + 1);
      gRoutside2->SetLineColor(kCyan + 1);
      gRoutside2->Fit(fitSinRoutside2, "QN");
      gRoutside2->Write();
      fitSinRoutside2->Write();

      string routlongName = groupKey + "_Routlong2_vs_phi";
      gRoutlong2->SetName(routlongName.c_str());
      gRoutlong2->SetTitle(
          (routlongName +
           ";#phi_{pair}-#Psi_{EP} (rad);R_{outlong}^{2} (fm^{2})")
              .c_str());
      gRoutlong2->SetMarkerStyle(29);
      gRoutlong2->SetMarkerColor(kViolet + 1);
      gRoutlong2->SetLineColor(kViolet + 1);
      gRoutlong2->Fit(fitCosRoutlong2, "QN");
      gRoutlong2->Write();
      fitCosRoutlong2->Write();

      string rsidelongName = groupKey + "_Rsidelong2_vs_phi";
      gRsidelong2->SetName(rsidelongName.c_str());
      gRsidelong2->SetTitle(
          (rsidelongName +
           ";#phi_{pair}-#Psi_{EP} (rad);R_{sidelong}^{2} (fm^{2})")
              .c_str());
      gRsidelong2->SetMarkerStyle(47);
      gRsidelong2->SetMarkerColor(kAzure + 2);
      gRsidelong2->SetLineColor(kAzure + 2);
      gRsidelong2->Fit(fitSinRsidelong2, "QN");
      gRsidelong2->Write();
      fitSinRsidelong2->Write();
    }

    string alphaName = groupKey + "_alpha_vs_phi";
    gAlpha->SetName(alphaName.c_str());
    gAlpha->SetTitle(
        (alphaName + ";#phi_{pair}-#Psi_{EP} (rad);#alpha").c_str());
    gAlpha->SetMarkerStyle(23);
    gAlpha->SetMarkerColor(kMagenta + 1);
    gAlpha->SetLineColor(kMagenta + 1);
    gAlpha->Fit(fitConstAlpha, "QN");
    gAlpha->Write();
    fitConstAlpha->Write();

    if (usesCoreHaloLambda) {
      string lambdaName = groupKey + "_lambda_vs_phi";
      gLambda->SetName(lambdaName.c_str());
      gLambda->SetTitle(
          (lambdaName + ";#phi_{pair}-#Psi_{EP} (rad);#lambda").c_str());
      gLambda->SetMarkerStyle(33);
      gLambda->SetMarkerColor(kOrange + 7);
      gLambda->SetLineColor(kOrange + 7);
      gLambda->Fit(fitConstLambda, "QN");
      gLambda->Write();
      fitConstLambda->Write();
    }

    delete gRout2;
    delete gRside2;
    delete gRlong2;
    delete gAlpha;
    delete fitCosRout2;
    delete fitCosRside2;
    delete fitCosRlong2;
    delete fitConstAlpha;
    if (usesCoreHaloLambda) {
      delete gLambda;
      delete fitConstLambda;
    }

    if (hasOffDiagonal) {
      delete gRoutside2;
      delete gRoutlong2;
      delete gRsidelong2;
      delete fitSinRoutside2;
      delete fitCosRoutlong2;
      delete fitSinRsidelong2;
    }
  }

  wf->cd();
}

void FitCF3DWithLevy(const string &rpath, const string &rfilename,
                     const string &wpath, const string &wfilename,
                     const string &txtFilename,
                     const Levy3DFitOptions &fitOptions = Levy3DFitOptions()) {
  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!InitializeROOTFile(wpath, wfilename)) {
    if (rf) {
      rf->Close();
      delete rf;
    }
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  if (!rf || rf->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file " << rpath << "/" << rfilename
         << ".root" << endl;
    delete rf;
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  vector<Levy3DFitResult> fitResults;
  TIter nextKey(rf->GetListOfKeys());
  TKey *key = nullptr;

  while ((key = (TKey *)nextKey())) {
    string objectName = key->GetName();
    if (!EndsWith(objectName, "_CF3D")) {
      continue;
    }

    auto hCF = LoadStoredCFHistogram(rf, objectName);
    if (!hCF) {
      continue;
    }

    Levy3DFitResult fitResult;
    if (!ParseCF3DHistogramName(objectName, fitResult)) {
      cout << "WARNING: cannot parse histogram name: " << objectName << endl;
      delete hCF;
      continue;
    }

    cout << "Fitting " << objectName << endl;
    if (FitAndWriteSingleCFHistogram(hCF, fitResult, wpath, wfilename, false,
                                     fitOptions)) {
      fitResults.push_back(fitResult);
    }
    delete hCF;
  }

  auto wf = GetROOT(wpath, wfilename, "update");
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot update output ROOT file " << wpath << "/"
         << wfilename << ".root for summary graphs" << endl;
    delete wf;
  } else {
    WriteR2Graphs(wf, fitResults);
    wf->Close();
    delete wf;
  }
  WriteFitResultsSummary(wpath + "/" + txtFilename + ".txt", fitResults);

  rf->Close();
  delete rf;
  cout << "Diagonal 3D Levy fit results have been written to " << wpath << "/"
       << wfilename << ".root and " << wpath << "/" << txtFilename << ".txt"
       << endl;
  TH1::AddDirectory(oldAddDirectory);
}

void FitCF3DWithFullLevy(
    const string &rpath, const string &rfilename, const string &wpath,
    const string &wfilename, const string &txtFilename,
    const Levy3DFitOptions &fitOptions = Levy3DFitOptions()) {
  const bool oldAddDirectory = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!InitializeROOTFile(wpath, wfilename)) {
    if (rf) {
      rf->Close();
      delete rf;
    }
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  if (!rf || rf->IsZombie()) {
    cout << "ERROR: cannot open input ROOT file " << rpath << "/" << rfilename
         << ".root" << endl;
    delete rf;
    TH1::AddDirectory(oldAddDirectory);
    return;
  }

  vector<Levy3DFitResult> fitResults;
  TIter nextKey(rf->GetListOfKeys());
  TKey *key = nullptr;

  while ((key = (TKey *)nextKey())) {
    string objectName = key->GetName();
    if (!EndsWith(objectName, "_CF3D")) {
      continue;
    }

    auto hCF = LoadStoredCFHistogram(rf, objectName);
    if (!hCF) {
      continue;
    }

    Levy3DFitResult fitResult;
    if (!ParseCF3DHistogramName(objectName, fitResult)) {
      cout << "WARNING: cannot parse histogram name: " << objectName << endl;
      delete hCF;
      continue;
    }

    cout << "Fitting full Levy model for " << objectName << endl;
    if (FitAndWriteSingleCFHistogram(hCF, fitResult, wpath, wfilename, true,
                                     fitOptions)) {
      fitResults.push_back(fitResult);
    }
    delete hCF;
  }

  auto wf = GetROOT(wpath, wfilename, "update");
  if (!wf || wf->IsZombie()) {
    cout << "ERROR: cannot update output ROOT file " << wpath << "/"
         << wfilename << ".root for summary graphs" << endl;
    delete wf;
  } else {
    WriteR2Graphs(wf, fitResults);
    wf->Close();
    delete wf;
  }
  WriteFitResultsSummary(wpath + "/" + txtFilename + ".txt", fitResults);

  rf->Close();
  delete rf;
  cout << "Full 3D Levy fit results have been written to " << wpath << "/"
       << wfilename << ".root and " << wpath << "/" << txtFilename << ".txt"
       << endl;
  TH1::AddDirectory(oldAddDirectory);
}

void _3d_cf_from_exp() {
  // string suffix = "_its";
  // string suffix = "_2400bins";
  string suffix = "_no_epmix";
  // string suffix = "_cent_mix";
  // string data_set = "23zzi_pass5";
  // string data_set = "23zzh_pass5_noepmix";
  // string data_set = "23zzm_pass5_medium";
  // string data_set = "23zzh_pass5_t3_morebin";
  // string data_set = "23zzh_pass5_t4_morebin";
  //   string data_set = "23zzh_pass5_t4_fulltof";
  // string data_set = "23zzh_pass5_t3_morebin";
  // string data_set = "23zzh_pass5_small_t3";
  // string data_set = "23zzh_pass5_small_dbgt4";
  // string data_set = "23zzh_pass5_small_noepmix";
  string data_set = "23zzh_pass5_medium";

  // string data_set = "25ae_pass2";
  // string data_set = "25ae_pass2_test";
  // string data_set = "25ae_pass2_noepmix";
  // string data_set = "25ae_pass2_noepmix_moremult";
  // string data_set = "25ae_pass2_epmix";
  // string data_set = "25ae_pass2_moredepth";

  // string rpath1 = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  string rpath_pbpb =
      "/Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/PbPb";
  string rpath_oo = "/Users/allenzhou/ALICE/alidata/hyperloop_res/femtoep/OO";
  string wpath_pbpb = "/Users/allenzhou/ALICE/alidata/femtoep_res/PbPb";
  string wpath_oo = "/Users/allenzhou/ALICE/alidata/femtoep_res/OO";

  string wpath, rpath;

  // enum system = {"PbPb", "OO"}; forbidden
  //   string rpath = rpath_pbpb;
  bool kisoo = 0;
  bool kispbpb = 1;

  if (kisoo && kispbpb) {
    cout << "ERROR: both kisoo and kispbpb cannot be true at the same time."
         << endl;
    return;
  }
  if (kisoo) {
    rpath = rpath_oo;
    wpath = wpath_oo;
  } else if (kispbpb) {
    rpath = rpath_pbpb;
    wpath = wpath_pbpb;
  }

  // string rpath2 = "/Users/allenzhou/ALICE/scripts/femtoep";
  // string rpath2 = rpath1;
  // string rpath2 = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  // string rfilename = "23_zzh_pass5_small_3";
  // string rfilename = "AnalysisResults_23zzf_pass5";
  // string rfilename = "AnalysisResults_" + data_set;
  string rfilename = "AnalysisResults_3d_" + data_set;
  // string rtaskname = "femto-dream-pair-task-track-track";
  string rtaskname_main = "femto-dream-pair-task-track-track";
  string rtaskname_sub = "femto-dream-pair-task-track-track" + suffix;
  string rsubtask_se = "SameEvent_3Dqn";
  string rsubtask_me = "MixedEvent_3Dqn";

  // string wpath = "/Users/allenzhou/ALICE/scripts/femtoep";
  string wfilename_main = "EP_dependence_CF_" + data_set;
  string wfilename_sub = "EP_dependence_CF_" + data_set + suffix;
  // string wfilename2 = "Check EP mixing 2";

  string wfilename, rtaskname;

  bool kis_subwagon = 0;

  if (kis_subwagon) {
    rtaskname = rtaskname_sub;
    wfilename = wfilename_sub;
  } else {
    rtaskname = rtaskname_main;
    wfilename = wfilename_main;
  }

  // std::vector<std::pair<double, double>> centBins = {
  //     {0, 30}, {30, 70}, {70, 100}};
  // std::vector<std::pair<double, double>> mTBins = {
  //     {0.5, 0.7}, {0.7, 1.0}, {1.0, 1.5}};
  std::vector<std::pair<double, double>> centBins = {
      {0, 10}, {10, 30}, {30, 50}, {50, 80}, {80, 100}};
  // std::vector<std::pair<double, double>> mTBins_kt = {{0.423, 8.001}};
  std::vector<std::pair<double, double>> mTBins_kt = {{0.244131, 0.331059},
                                                      {0.331059, 0.423792},
                                                      {0.423792, 0.51923},
                                                      {0.51923, 0.713863},
                                                      {0.244131, 0.713863}};
  // std::vector<std::pair<double, double>> mTBins_kt = {
  //     {0.423, 0.519}, {0.519, 0.812}, {0.812, 1.208},
  //     {1.208, 1.606}, {2.005, 2.504}, {2.504, 3.003},
  //     {3.004, 4.002}, {4.002, 5.002}, {5.002, 8.001}};
  std::vector<std::pair<double, double>> mTBins_mt = {
      {0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7},
      {0.7, 0.8}, {0.8, 1.0}, {1.0, 1.2}, {1.2, 1.6}, {1.6, 2}};
  //   std::vector<std::pair<double, double>> mTBins = {{0.244131, 0.331059},
  //                                                    {0.331059, 0.423792},
  //                                                    {0.423792, 0.51923},
  //                                                    {0.51923, 0.713863},
  //                                                    {0.244131, 0.713863}};
  std::vector<std::string> pairphiBins = {"In_plane", "Out_of_plane"};

  std::vector<EPBin> epBins_se = {
      {0, TMath::Pi() / 4, 3 * TMath::Pi() / 4, TMath::Pi(), "In_plane"},
      {TMath::Pi() / 4, 3 * TMath::Pi() / 4, -1, -1, "Out_of_plane"}};
  std::vector<EPBin> epBins_me = {{{0, TMath::Pi(), -1, -1, "Min bias EP"}}};

  bool doBuildCF3D = true;
  bool doFitDiag = false;
  bool doFitFull = false;
  bool fitUseCoulomb = true;
  bool fitUseCoreHaloLambda = true;

  std::vector<std::pair<double, double>> fitCentBins = {{0, 10}, {10, 30}};
  // 下面mtbin留两位数就行
  std::vector<std::pair<double, double>> fitMTBins_mt = {
      {0.20, 0.30},
      {0.30, 0.40},
      {0.40, 0.50},
  };
  // std::vector<std::pair<double, double>> fitMTBins_kt = {{0.42, 0.81}};
  // std::vector<std::pair<double, double>> fitMTBins_kt = {{0.42, 0.52},
  //                                                        {0.52, 0.81}};
  std::vector<std::pair<double, double>> fitMTBins_kt = {
      {0.42, 0.52}, {0.52, 0.71}, {0.24, 0.71}};

  std::vector<std::pair<double, double>> fitMTBins = fitMTBins_kt;
  std::vector<std::pair<double, double>> mTBins = mTBins_kt;

  Levy3DFitOptions fitOptions;
  fitOptions.useCoulomb = fitUseCoulomb;
  fitOptions.useCoreHaloLambda = fitUseCoreHaloLambda;
  fitOptions.fitQMax = 0.15;
  const string fitOptionTag = BuildFitOptionTag(fitOptions);

  string fitInputPath = wpath;
  string fitInputFilename = wfilename;
  string fitDiagOutput = "EP_dependence_CF_fit_" + fitOptionTag + "_" +
                         data_set + (kis_subwagon ? suffix : "");
  string fitDiagTxt = "EP_dependence_CF_fit_params_" + fitOptionTag + "_" +
                      data_set + (kis_subwagon ? suffix : "");
  string fitFullOutput = "EP_dependence_CF_full_fit_" + fitOptionTag + "_" +
                         data_set + (kis_subwagon ? suffix : "");
  string fitFullTxt = "EP_dependence_CF_full_fit_params_" + fitOptionTag + "_" +
                      data_set + (kis_subwagon ? suffix : "");

  //   CFCalcWith_Cent_Mt_pairphi_full(rpath, rfilename, rtaskname, rsubtask_se,
  //                                   rsubtask_me, wpath, wfilename, centBins,
  //                                   mTBins, epBins_se, epBins_me,
  //                                   pairphiBins, {0.5, 0.8}, {0., 0.8});
  if (doBuildCF3D) {
    CFCalc3D(rpath, rfilename, rtaskname, rsubtask_se, rsubtask_me, wpath,
             wfilename, centBins, mTBins);
  }

  if (doFitDiag) {
    FitCF3DWithSelectedBins(fitInputPath, fitInputFilename, wpath,
                            fitDiagOutput, fitDiagTxt, fitCentBins, fitMTBins,
                            false, fitOptions);
  }

  if (doFitFull) {
    FitCF3DWithSelectedBins(fitInputPath, fitInputFilename, wpath,
                            fitFullOutput, fitFullTxt, fitCentBins, fitMTBins,
                            true, fitOptions);
  }
}
