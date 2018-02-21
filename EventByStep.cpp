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
	gStyle->SetTitleXOffset(0.9);
	gStyle->SetTitleYOffset(0.9);

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
	c1->RedrawAxis();

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

	float ev_sig[n] = {221758+36466,221758+36466,157706+26337,160616+27243,125519+20458};
	float ev_wjet[n] = {11692,11692,1130,1154,0};
	float ev_STop[n] = {25389,25389,9574,10006,6711};
	float ev_Diboson[n] = {29945,29945,2186,2327,200};
	float ev_DY[n] = {80763,80763,5949,6100,623};
	float ev_data[n] = {374170,374170,189973,190833,139426};

	TH1F *hist_sig = new TH1F("hist_sig","", n-1, x);
	TH1F *hist_wjet = new TH1F("hist_wjet","",n-1, x);
	TH1F *hist_STop = new TH1F("hist_STop","",n-1, x);
	TH1F *hist_Diboson = new TH1F("hist_Diboson","",n-1, x);
	TH1F *hist_DY = new TH1F("hist_DY","",n-1, x);
	TH1F *hist_data = new TH1F("hist_data","",n-1, x);

	double er_data = 0;
	for(int i = 0; i < n; i++){
		hist_sig->SetBinContent(i+1, x[i]);
		hist_sig->SetBinContent(i+1, ev_sig[i]);
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
	//hist_others->SetLineColor(906);
	//hist_others->SetFillColor(906);
	hist_wjet->SetLineColor(3);
	hist_wjet->SetFillColor(3);
	hist_STop->SetLineColor(42);
	hist_STop->SetFillColor(42);
	hist_Diboson->SetLineColor(7);
	hist_Diboson->SetFillColor(7);
	hist_DY->SetLineColor(4);
	hist_DY->SetFillColor(4);

	hist_data->SetBinError(n,sqrt(er_data));
	hist_data->SetMarkerSize(1.1);
	hist_data->SetMarkerStyle(20);

	hist_data->GetXaxis()->SetBinLabel(1,"Dilepton");
	hist_data->GetXaxis()->SetBinLabel(2,"Z veto");
	hist_data->GetXaxis()->SetBinLabel(3,">1 jets");
	hist_data->GetXaxis()->SetBinLabel(4,"MET");
	hist_data->GetXaxis()->SetBinLabel(5,">0 b jets");
	//hist->Scale(100/1);
	//hist->SetMinimum(0);

	l_->AddEntry(hist_sig,"t#bar{t}","f");
	//l_->AddEntry(hist_others,"tt_others","f");
	l_->AddEntry(hist_wjet,"W+jet","f");
	l_->AddEntry(hist_STop,"Single Top","f");
	l_->AddEntry(hist_Diboson,"Diboson","f");
	l_->AddEntry(hist_DY,"Z/#gamma^{*}#rightarrow#font[12]{l#lower[-0.4]{+}l#lower[-0.4]{#font[122]{\55}}}","f");
	l_->AddEntry(hist_data,"data","lp");

	hs = new THStack();

	hs->Add(hist_DY);
	hs->Add(hist_Diboson);
	hs->Add(hist_STop);
	hs->Add(hist_wjet);
	//hs->Add(hist_others);
	hs->Add(hist_sig);

	double ymax = 0;
	ymax = hist_data->GetMaximum();
	hist_data->SetMaximum(ymax*1.3);
	hist_data->SetMinimum(0);
	hist_data->GetYaxis()->SetTitle("Number of Events");
	//hist_data->GetXaxis()->SetTitle("Number of Steps");

	hist_data->Draw();
	hs->Draw("histsame");
	hist_data->Draw("esame");

	lt1.DrawLatex(xx_1,yy_1,"e^{#pm}#mu^{#mp}");
	lt2.DrawLatex(x_1,y_1,"CMS");
	lt3.DrawLatex(x_2,y_2,"Preliminary");
	lt4.DrawLatex(tx,ty,"35.9 fb^{-1}, #sqrt{s} = 13 TeV");
	l_->Draw();
	c1->SaveAs("/afs/cern.ch/work/y/yjeong/catMacro/plots/EventByStep.png");
}
