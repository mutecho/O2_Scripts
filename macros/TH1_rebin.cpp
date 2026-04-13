#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TMath.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include "TF1.h"
#include "TH1.h"
#include "TTree.h"

using namespace std;

bool FileExists_bool(const string &filename) {
  ifstream file(filename);
  return file.good();
}

void FileExists_warn(const string &filename) {
  if (!FileExists_bool(filename)) {
    cerr << "Error: File: " << filename << " doesn't exist!" << endl;
  }
}

TFile *GetROOT(const string Readpath, const string ReadFilename, const string Option) {
  string lowerOption = Option;
  transform(lowerOption.begin(), lowerOption.end(), lowerOption.begin(), [](unsigned char c) {
    return tolower(c);
  });  // 全转小写
  if (!(lowerOption == "read" || lowerOption == "write" || lowerOption == "recreate")) {
    cout << "Invalid Option: " << Option << ", should be one of: read, write, recreate" << endl;
  }
  string rfilename = Readpath + "/" + ReadFilename + ".root";  // 加上root后缀
  FileExists_warn(rfilename);

  return new TFile(rfilename.c_str(), Option.c_str());
}

TH1D *calc_cf_from_sme_rerange(TH1D *h_se,
                               TH1D *h_me,
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
  double factorN = h_me_c->Integral(binNorm[0], binNorm[1]) / h_se_c->Integral(binNorm[0], binNorm[1]);
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
    int newBin = i - binRange[0] + 1;  // 从1开始填
    h_cf_re->SetBinContent(newBin, content);
    h_cf_re->SetBinError(newBin, error);
  }
  return h_cf_re;
}

TH1D *TH1D_rebin(TH1D *h, int rebinFactor = 2) {
  int nBins = h->GetNbinsX();
  int nBinsRebinned = nBins / rebinFactor;

  TH1D *hRebinned = dynamic_cast<TH1D *>(h->Rebin(rebinFactor, "h_rebinned"));
  //   TH1D *hRebinned = new TH1D("", "", nBinsRebinned,
  //   h->GetXaxis()->GetXmin(),
  //                              h->GetXaxis()->GetXmax());

  //   for (int i = 0; i < nBinsRebinned; i++) {
  //     double contentSum = 0.0;
  //     double errorSumSq = 0.0;

  //     for (int j = 0; j < rebinFactor; j++) {
  //       int binIndex = i * rebinFactor + j + 1; // +1 because bins start at 1
  //       contentSum += h->GetBinContent(binIndex);
  //       errorSumSq += pow(h->GetBinError(binIndex), 2);
  //     }

  //     hRebinned->SetBinContent(i + 1, contentSum);
  //     hRebinned->SetBinError(i + 1, sqrt(errorSumSq));
  //   }

  return hRebinned;

  // Copy rebinned content back to original histogram
  //   h->Reset();
  //   h->SetBins(nBinsRebinned, hRebinned->GetXaxis()->GetXmin(),
  //              hRebinned->GetXaxis()->GetXmax());

  //   for (int i = 0; i < nBinsRebinned; i++) {
  //     h->SetBinContent(i + 1, hRebinned->GetBinContent(i + 1));
  //     h->SetBinError(i + 1, hRebinned->GetBinError(i + 1));
  //   }

  //   delete hRebinned;
}

void TH1_rebin() {
  string suffix = "_2400bins";
  string data_set = "23zzh_pass5_morebin";
  string rpath = "/Users/allenzhou/ALICE/alidata/femtoep_res";

  string rfilename = "EP_dependence_CF_23zzh_pass5_morebin";
  //   string rtaskname = "cent=0-10_mT=0.24-0.33/CF_reranged_Minbias_EP";
  // string rtaskname = "cent=0-10_mT=0.24-0.33_CF_reranged_Minbias_EP";
  string rtaskname1 = "cent=0-10_mT=0.2-0.3_SE_Minbias_EP";
  string rtaskname2 = "cent=0-10_mT=0.2-0.3_ME_Minbias_EP";

  // string wpath = "/Users/allenzhou/ALICE/scripts/femtoep";
  string wpath = "/Users/allenzhou/ALICE/alidata/femtoep_res";
  // string wfilename1 = "EP_dependence_CF_" + data_set;
  string wfilename = "test";

  auto rf = GetROOT(rpath, rfilename, "read");
  if (!rf || rf->IsZombie()) {
    cout << "Error opening file: " << rpath + "/" + rfilename + ".root" << endl;
    return;
  }

  auto wf = GetROOT(wpath, wfilename, "recreate");
  if (!wf || wf->IsZombie()) {
    cout << "Error creating file: " << wpath + "/" + wfilename + ".root" << endl;
    rf->Close();
    return;
  }

  TH1D *hs = (TH1D *)rf->Get(rtaskname1.c_str());
  if (!hs) {
    cout << "Error: Histogram " << rtaskname1 << " not found in file " << rfilename << endl;
    rf->Close();
    return;
  }
  // hs->Sumw2();
  wf->WriteObject(hs, "original_same");

  TH1D *hs_rebinned = TH1D_rebin(hs, 5);
  if (!hs_rebinned) {
    cout << "Error: Rebinning failed." << endl;
    rf->Close();
    return;
  }
  wf->WriteObject(hs_rebinned, "rebinned_same");

  TH1D *hm = (TH1D *)rf->Get(rtaskname2.c_str());
  if (!hm) {
    cout << "Error: Histogram " << rtaskname2 << " not found in file " << rfilename << endl;
    rf->Close();
    return;
  }
  // hm->Sumw2();
  wf->WriteObject(hm, "original_mixed");

  TH1D *hm_rebinned = TH1D_rebin(hm, 5);
  if (!hm_rebinned) {
    cout << "Error: Rebinning failed." << endl;
    rf->Close();
    return;
  }
  wf->WriteObject(hm_rebinned, "rebinned_mixed");

  auto *hcf_rebinned = calc_cf_from_sme_rerange(hs_rebinned, hm_rebinned, {0.5, 0.8}, {0., 0.25});
  hcf_rebinned->Sumw2();
  // hcf_rebinned->Draw();

  wf->WriteObject(hcf_rebinned, "rebinned_CF");
  return;
}