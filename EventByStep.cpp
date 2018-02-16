{
	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);//To display the mean and RMS: SetOptStat("mr"), nemruoi, ;
	gStyle->SetOptDate(0);//display date position
	/*gStyle->SetCanvasDefH(600);//Height of canvas
	  gStyle->SetCanvasDefW(600);//Width of canvas
	  gStyle->SetCanvasDefX(0);//Position on screen
	  gStyle->SetCanvasDefY(0);*/
	gStyle->SetPalette(0);

	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetPadTopMargin(0.07);
	gStyle->SetPadBottomMargin(0.1);
	gStyle->SetPadBorderMode(0);

	gStyle->SetLabelFont(42, "XY");
	gStyle->SetLabelOffset(0.007, "XY");
	gStyle->SetLabelSize(0.04, "Y");
	gStyle->SetLabelSize(0.06, "X");

	gStyle->SetTitleFont(42, "XY");
	gStyle->SetTitleSize(0.06, "X");
	gStyle->SetTitleSize(0.06, "Y");
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetTitleYOffset(1.15);

	gStyle->SetAxisColor(1, "XYZ");
	gStyle->SetTickLength(0.03, "XYZ");
	gStyle->SetNdivisions(510, "XYZ");
	gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	gStyle->SetPadTickY(1);

	TLatex lt1;
	lt1.SetTextAlign(12);
	lt1.SetNDC();
	lt1.SetTextFont(132);
	lt1.SetTextAngle(0);
	lt1.SetTextSize(0.045);

	TLatex lt2;
	lt2.SetTextAlign(12);
	lt2.SetNDC();
	lt2.SetTextFont(61);
	lt2.SetTextAngle(0);
	lt2.SetTextSize(0.058);

	TLatex lt3;
	lt3.SetTextAlign(12);
	lt3.SetNDC();
	lt3.SetTextAngle(0);
	lt3.SetTextFont(52);
	lt3.SetTextSize(0.045);

	TLatex lt4;
	lt4.SetTextAlign(32);
	lt4.SetNDC();
	lt4.SetTextAngle(0);
	lt4.SetTextFont(42);
	lt4.SetTextSize(0.05);

	THStack *hs;
	TLegend *l_;

	float x_1 = 0.12; //right top side x_1 = 0.73
	const float y_1 = 0.97; //right top side y_1 = 0.84
	float x_2 = x_1+0.095; //right top side y_2 = y_1-0.07
	float y_2 = y_1-0.005;
	//------------------------------------Coordinate of first LatexBox---------------------------------------
	float xx_1 = 0.18;
	float yy_1 = 0.83;
	//--------------------------------------Set Maximum histo_TTTT[NVar][NStep][nCh]---------------------------------
	//float ymin_1 = 0;
	//-----------------------------------------ExtraText---------------------------------------
	float tx = 0.98;
	float ty = 0.97;
	//-------------------------Legend coordinate--------------------
	float lx1 = 0.65;
	float ly1 = 0.59;
	float lx2 = 0.95;
	float ly2 = 0.86;

	c1 = new TCanvas;
	//c1->SetLogy();
	l_ = new TLegend(lx1,ly1,lx2,ly2);
	l_->SetFillColor(0);
	l_->SetLineColor(0);
	l_->SetLineStyle(kSolid);
	l_->SetLineWidth(1);
	l_->SetFillStyle(1001);
	l_->SetTextFont(42);
	l_->SetTextSize(0.035);

	const int n = 6;
	float x[n] = {1,2,3,4,5,6};
	float ev_sig[n] = {215279,215279,157703,157703,125524,800};
	float ev_others[n] = {35959,35959,26336,26336,20459,700};
	float ev_wjet[n] = {11421,11421,1130,1130,0,700};
	float ev_STop[n] = {24588,24588,9575,9575,6711,300};
	float ev_Diboson[n] = {28938,28938,2186,2186,200,30};
	float ev_DY[n] = {80722,80722,5950,5950,623,90};
	float ev_data[n] = {366556,366556,186131,186131,139439,90};

	TH1F *hist_sig = new TH1F("hist_sig","", n-1, x);
	TH1F *hist_others = new TH1F("hist_others","",n-1, x);
	TH1F *hist_wjet = new TH1F("hist_wjet","",n-1, x);
	TH1F *hist_STop = new TH1F("hist_STop","",n-1, x);
	TH1F *hist_Diboson = new TH1F("hist_Diboson","",n-1, x);
	TH1F *hist_DY = new TH1F("hist_DY","",n-1, x);
	TH1F *hist_data = new TH1F("hist_data","",n-1, x);

	//double er_data[n] = {0,};
	double er_data = 0;
	for(int i = 0; i < n; i++){
		hist_sig->SetBinContent(i+1, x[i]);
		hist_sig->SetBinContent(i+1, ev_sig[i]);
		hist_others->SetBinContent(i+1, x[i]);
		hist_others->SetBinContent(i+1, ev_others[i]);
		hist_wjet->SetBinContent(i+1, x[i]);
		hist_wjet->SetBinContent(i+1, ev_wjet[i]);
		hist_STop->SetBinContent(i+1, x[i]);
		hist_STop->SetBinContent(i+1, ev_STop[i]);
		hist_Diboson->SetBinContent(i+1, x[i]);
		hist_Diboson->SetBinContent(i+1, ev_Diboson[i]);
		hist_DY->SetBinContent(i+1, x[i]);
		hist_DY->SetBinContent(i+1, ev_DY[i]);
		hist_data->SetBinContent(i+1, x[i]);
		hist_data->SetBinContent(i+1, ev_data[i]);
		//ev[i] = hist->GetBinContent(i+1);
		er_data = hist_data->GetBinContent(i+1);
		//hist->SetBinError(i+1,sqrt(ev));
		//hist->SetBinError(i+1,sqrt(y[i]));
		//hist->GetYaxis()->SetRangeUser(10.,100.);
	}

	/*for(int i = 0; i < n; i++){
	//hist_others->Fill(x[i],y_[i]);
	//hist_others->SetBinError(i+1,sqrt(y[i]));
	}*/

	hist_sig->SetLineColor(2);
	hist_sig->SetFillColor(2);
	hist_others->SetLineColor(906);
	hist_others->SetFillColor(906);
	hist_wjet->SetLineColor(3);
	hist_wjet->SetFillColor(3);
	hist_STop->SetLineColor(42);
	hist_STop->SetFillColor(42);
	hist_Diboson->SetLineColor(7);
	hist_Diboson->SetFillColor(7);
	hist_DY->SetLineColor(4);
	hist_DY->SetFillColor(4);
	hist_data->SetLineColor(4);
	hist_data->SetFillColor(4);

	hist_data->SetBinError(n,sqrt(er_data));
	hist_data->SetMarkerSize(1.1);
	hist_data->SetMarkerStyle(20);

	hist_sig->GetXaxis()->SetBinLabel(1,"step1");
	hist_sig->GetXaxis()->SetBinLabel(2,"step2");
	hist_sig->GetXaxis()->SetBinLabel(3,"step3");
	hist_sig->GetXaxis()->SetBinLabel(4,"step4");
	hist_sig->GetXaxis()->SetBinLabel(5,"step5");
	//hist->Scale(100/1);
	//hist->SetMinimum(0);

	l_->AddEntry(hist_sig,"sig","f");
	l_->AddEntry(hist_others,"others","f");
	l_->AddEntry(hist_wjet,"wjet","f");
	l_->AddEntry(hist_STop,"Stop","f");
	l_->AddEntry(hist_Diboson,"Diboson","f");
	l_->AddEntry(hist_DY,"DY","f");
	l_->AddEntry(hist_data,"data","f");

	hs = new THStack();

	hs->Add(hist_sig);
	hs->Add(hist_others);
	hs->Add(hist_wjet);
	hs->Add(hist_STop);
	hs->Add(hist_Diboson);
	hs->Add(hist_DY);
	hs->Draw("hist");
	hist_data->Draw("same");

	lt1.DrawLatex(xx_1,yy_1,"e^{#pm}#mu^{#mp}");
	lt2.DrawLatex(x_1,y_1,"CMS");
	lt3.DrawLatex(x_2,y_2,"Preliminary");
	lt4.DrawLatex(tx,ty,"35.9 fb^{-1}, #sqrt{s} = 13 TeV");
	l_->Draw();
	c1->SaveAs("/afs/cern.ch/work/y/yjeong/catMacro/plots/StepByEvents.png");
}
