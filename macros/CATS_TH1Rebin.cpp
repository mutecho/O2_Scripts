// 本文件用于从只含一张直方图（TH1）的图中，重新分bin

#include "TRandom3.h"
#include <cstdio>
#include <cfloat>
#include <cstring>
#include "TF1.h"
#include "TProfile.h"
#include "TAxis.h"
#include "TLine.h"
#include "TString.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#define pi_s 3.14159265358979323844
// using std::cout;
// using std::endl;
using namespace std;

void TH1Rebin()
{
    double bincenter, CFvalue;
    // string filename = "PLambda-corr_1260.0-1320.0__nocoul_pp13_EPOS_list"; // 选择要读取数据的root文件名字
    // string filename = "AMPT_CATS_bs"; // 选择要读取数据的root文件名字
    string filename = "AMPT_CATS"; // 选择要读取数据的root文件名字
    // string filename = "PLambda-corr_bs_1260.0-1320.0__nocoul_pp13_EPOS_list"; // 选择要读取数据的root文件名字
    // string readpath = "/storage/fdunphome/zhouziqian/MyTools/EPOS_CATS/result/pp13_EPOS_list";
    string readpath = "/storage/fdunphome/zhouziqian/MyTools/data/CATS";
    // string writepath = "/storage/fdunphome/zhouziqian/MyTools/data/CATS/AMPT_pp13_3q1.5";
    TString rfilename = readpath + "/" + filename + ".root"; // 加上root后缀
    unique_ptr<TFile> rootf = make_unique<TFile>(rfilename, "update");
    // TFile fnew(rfilename);
    // TH1D *CF = (TH1D*)fnew.Get("CF_0");
    TProfile *CF = (TProfile *)rootf->Get("correlation_EPOS__nocoul");
    // TProfile *CF = (TProfile *)rootf->Get("correlation_EventGen__nocoul");
    TH1 *CF_new = CF->Rebin(20, "rebinned cf");
    CF_new->Write();
    cout << "文件：" << filename << ".root重写完毕" << endl;
    rootf->Close();
}