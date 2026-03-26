#include <TFile.h>
#include <THnSparse.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <iostream>

void femtoep()
{
    // 1. 打开文件
    TFile *f = TFile::Open("/Users/allenzhou/ALICE/scripts/femtoep/AnalysisResults_femtoep_sm.root");
    auto h4 = (THnSparseF*)f->Get("femto-dream-pair-task-track-track/SameEvent_EP/relPairkstarmTMultMultPercentileQn");

    // 2. 获取每个维度的轴指针
    auto ax_kstar = h4->GetAxis(0);
    auto ax_mT    = h4->GetAxis(1);
    auto ax_cent  = h4->GetAxis(2);
    auto ax_EP    = h4->GetAxis(3);

    // 定义中心度和 mT 区间
    std::vector<std::pair<double,double>> centBins = {{0,10}, {10,30}, {30,50}};
    std::vector<std::pair<double,double>> mTBins   = {{0.5,0.7}, {0.7,1.0}, {1.0,1.5}};

    // 定义 EP 对称区间
    struct EPBin {
        double low1, high1;
        double low2, high2;
        const char* name;
    };
    std::vector<EPBin> epBins = {
        {0, TMath::Pi()/4, 3*TMath::Pi()/4, TMath::Pi(), "In-plane"},
        {TMath::Pi()/4, 3*TMath::Pi()/4, -1, -1, "Out-of-plane"}
    };

    int colorList[] = {kRed+1, kBlue+1, kGreen+2, kOrange+7};

    // 主循环：每个中心度和 mT 各一张图
    for (auto [centLow, centHigh] : centBins) {
        for (auto [mTLow, mTHigh] : mTBins) {

            // 创建新的画布
            TString cname = Form("c_cent%.0f_%.0f_mT%.1f_%.1f", centLow, centHigh, mTLow, mTHigh);
            TCanvas *c = new TCanvas(cname, cname, 800, 600);
            c->SetMargin(0.13, 0.05, 0.12, 0.07);
            c->SetGrid();

            TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetTextSize(0.04);

            std::cout << Form("Centrality %.0f–%.0f, mT %.2f–%.2f", centLow, centHigh, mTLow, mTHigh) << std::endl;

            // 限制中心度与 mT
            ax_cent->SetRangeUser(centLow, centHigh);
            ax_mT->SetRangeUser(mTLow, mTHigh);

            int iColor = 0;

            for (auto &ep : epBins) {
                // 如果是双区间（in-plane），要合并两部分
                TH1D *h_kstar_total = nullptr;

                if (ep.low2 > 0) { // 有镜像部分
                    ax_EP->SetRangeUser(ep.low1, ep.high1);
                    auto h1 = (TH1D*)h4->Projection(0);
                    ax_EP->SetRangeUser(ep.low2, ep.high2);
                    auto h2 = (TH1D*)h4->Projection(0);
                    h1->Add(h2);
                    h_kstar_total = h1;
                } else {
                    ax_EP->SetRangeUser(ep.low1, ep.high1);
                    h_kstar_total = (TH1D*)h4->Projection(0);
                }

                // 设置样式
                h_kstar_total->SetLineColor(colorList[iColor++ % 4]);
                h_kstar_total->SetLineWidth(2);
                h_kstar_total->SetTitle(Form("Centrality %.0f–%.0f%%, mT %.2f–%.2f", centLow, centHigh, mTLow, mTHigh));
                h_kstar_total->GetXaxis()->SetTitle("k* (GeV/c)");
                h_kstar_total->GetYaxis()->SetTitle("Counts");

                // 归一化
                if (h_kstar_total->Integral() > 0)
                    h_kstar_total->Scale(1.0 / h_kstar_total->Integral());

                if (iColor == 1)
                    h_kstar_total->Draw("hist");
                else
                    h_kstar_total->Draw("hist same");

                leg->AddEntry(h_kstar_total, ep.name, "l");
            }

            leg->Draw();
            c->Update();
        }
    }
}