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
	gStyle->SetLabelOffset(0.006, "XY");
	gStyle->SetLabelSize(0.05, "XY");

	gStyle->SetTitleFont(42, "XY");
	gStyle->SetTitleXOffset(0.9);

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
	lt1.SetTextSize(0.075);

	TLatex lt2;
	lt2.SetTextAlign(12);
	lt2.SetNDC();
	lt2.SetTextFont(61);
	lt2.SetTextAngle(0);
	lt2.SetTextSize(0.066);

	TLatex lt3;
	lt3.SetTextAlign(12);
	lt3.SetNDC();
	lt3.SetTextAngle(0);
	lt3.SetTextFont(52);
	lt3.SetTextSize(0.053);

	TLatex lt4;
	lt4.SetTextAlign(32);
	lt4.SetNDC();
	lt4.SetTextAngle(0);
	lt4.SetTextFont(42);
	lt4.SetTextSize(0.065);

	float x_1 = 0.12; //right top side x_1 = 0.73
	const float y_1 = 0.97; //right top side y_1 = 0.84
	float x_2 = x_1+0.075; //right top side y_2 = y_1-0.07
	float y_2 = y_1-0.005;
	//------------------------------------Coordinate of first LatexBox---------------------------------------
	float xx_1 = 0.18;
	float yy_1 = 0.83;
	//--------------------------------------Set Maximum histo_TTTT---------------------------------
	//float ymin_1 = 0;
	//-----------------------------------------ExtraText---------------------------------------
	float tx = 0.98;
	float ty = 0.97;
	//-------------------------Legend coordinate--------------------
	float lx1 = 0.75;
	float ly1 = 0.54;
	float lx2 = 0.95;
	float ly2 = 0.86;

	TString Save_dir;
	Save_dir = "/afs/cern.ch/work/y/yjeong/catMacro/plots/";

	int Var=4;

	TString Variable;
	if(Var==0)Variable = "";//1
	if(Var==1)Variable = "M(ll)";//2
	if(Var==2)Variable = "MET";//3
	if(Var==3)Variable = "NJet";//4
	if(Var==4)Variable = "NBJet";//5

	THStack *hs;
	TLegend *l_;

	TH1F *hist_sig, *hist_wjet, *hist_STop, *hist_Diboson, *hist_DY, *hist_data;

	TH1F *hist_MC, *hist_Ratio;

	c1 = new TCanvas;
	//c1->SetLogy();
	c1->RedrawAxis();

	TPad *plotpad, *ratiopad;

	plotpad = new TPad("","",0.02,0.3,0.98,0.98);
	ratiopad = new TPad("","",0.02,0.1,0.98,0.3);

	plotpad->Draw();
	ratiopad->Draw();
	plotpad->RedrawAxis();

	plotpad->cd();
	//plotpad->RedrawAxis();
	gPad->SetBottomMargin(0);

	l_ = new TLegend(lx1,ly1,lx2,ly2);
	l_->SetFillColor(0);
	l_->SetLineColor(0);
	l_->SetLineStyle(kSolid);
	l_->SetLineWidth(1);
	l_->SetFillStyle(1001);
	l_->SetTextFont(42);
	l_->SetTextSize(0.055);

	const int n = 6;
	const int nVariable = 5;
	float x[n] = {1,2,3,4,5,6};

	float ev_sig[nVariable][5] = {
		{221758+36466,221758+36466,157706+26337,160616+27243,125519+20458},//total
		{221758+36466,221758+36466,162305+26705,162305+26705,129035+20733},//dilep.M
		{218657+37030,218657+37030,160616+27243,160616+27243,127721+21118},//MET
		{215279+35959,215279+35959,157706+26337,157706+26337,125526+20459},//NJet
		{215271+35958,215271+35958,157697+26335,157697+26335,125519+20458}//NBJet
	};
	float ev_wjet[nVariable][5] = {
		{11692,11692,1130,1154,0},
		{11692,11692,1237,1237,0},
		{11577,11577,1154,1154,0},
		{11420,11420,1130,1130,0},
		{11420,11420,1130,1130,0}
	};
	float ev_STop[nVariable][5] = {
		{25389,25389,9574,10006,6711},
		{25389,25389,9973,9973,6984},
		{25123,25123,10006,10006,7022},
		{24587,24587,9574,9574,6711},
		{24587,24587,9574,9574,6711}
	};
	float ev_Diboson[nVariable][5] = {
		{29945,29945,2186,2327,200},
		{29945,29945,2327,2327,213},
		{29221,29221,2327,2327,218},
		{28938,28938,2186,2186,200},
		{28938,28938,2186,2186,200}
	};
	float ev_DY[nVariable][5] = {
		{80763,80763,5949,6100,623},
		{80763,80763,5966,5966,626},
		{80963,80963,6100,6100,647},
		{80720,80720,5949,5949,623},
		{80720,80720,5949,5949,623}
	};
	float ev_data[nVariable][5] = {
		{374170,374170,189973,190833,139426},
		{374170,374170,190833,190833,142834},
		{371168,371168,189973,189973,142046},
		{366529,366529,186119,186119,139433},
		{366520,366520,186110,186110,139426}
	};//-*/

	hist_sig = new TH1F("hist_sig","", n-1, x);
	hist_wjet = new TH1F("hist_wjet","",n-1, x);
	hist_STop = new TH1F("hist_STop","",n-1, x);
	hist_Diboson = new TH1F("hist_Diboson","",n-1, x);
	hist_DY = new TH1F("hist_DY","",n-1, x);
	hist_data = new TH1F("hist_data","",n-1, x);

	double er_data = 0;
	for(int i = 0; i < n; i++){
		hist_sig->SetBinContent(i+1, x[i]);
		hist_sig->SetBinContent(i+1, ev_sig[Var][i]);
		hist_wjet->SetBinContent(i+1, x[i]);
		hist_wjet->SetBinContent(i+1, ev_wjet[Var][i]);
		hist_STop->SetBinContent(i+1, x[i]);
		hist_STop->SetBinContent(i+1, ev_STop[Var][i]);
		hist_Diboson->SetBinContent(i+1, x[i]);
		hist_Diboson->SetBinContent(i+1, ev_Diboson[Var][i]);
		hist_DY->SetBinContent(i+1, x[i]);
		hist_DY->SetBinContent(i+1, ev_DY[Var][i]);
		hist_data->SetBinContent(i+1, x[i]);
		hist_data->SetBinContent(i+1, ev_data[Var][i]);
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

	l_->AddEntry(hist_sig,"t#bar{t}","f");
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
	hs->Add(hist_sig);

	plotpad->cd();

	double ymax = 0;
	ymax = hist_data->GetMaximum();
	hist_data->SetMaximum(ymax*1.3);
	hist_data->SetMinimum(0.1);
	hist_data->GetYaxis()->SetTitle("Number of Events");
	hist_data->GetYaxis()->SetTitleOffset(0.7);
	hist_data->GetYaxis()->SetTitleSize(0.09);
	if(Var==0)hist_data->GetXaxis()->SetTitle("");
	if(Var>0) hist_data->GetXaxis()->SetTitle(Variable);

	hist_data->Draw();
	hs->Draw("histsame");
	hist_data->Draw("esame");

	lt1.DrawLatex(xx_1,yy_1,"e^{#pm}#mu^{#mp}");
	lt2.DrawLatex(x_1,y_1,"CMS");
	lt3.DrawLatex(x_2,y_2,"Preliminary");
	lt4.DrawLatex(tx,ty,"35.9 fb^{-1}, #sqrt{s} = 13 TeV");
	l_->Draw();

	hist_MC = new TH1F("hist_MC","",n-1, x);
	hist_MC->Add(hist_sig);
	hist_MC->Add(hist_wjet);
	hist_MC->Add(hist_STop);
	hist_MC->Add(hist_Diboson);
	hist_MC->Add(hist_DY);

	hist_Ratio = new TH1F("hist_Ratio","",n-1,x);
	hist_Ratio->Divide(hist_MC,hist_data,1,1,"b");

	if(Var==0){
		hist_Ratio->GetXaxis()->SetBinLabel(1,"Dilepton");
		hist_Ratio->GetXaxis()->SetBinLabel(2,"Z veto");
		hist_Ratio->GetXaxis()->SetBinLabel(3,">1 jets");
		hist_Ratio->GetXaxis()->SetBinLabel(4,"MET");
		hist_Ratio->GetXaxis()->SetBinLabel(5,">0 b jets");
	}

	if(Var>0){
		hist_Ratio->GetXaxis()->SetBinLabel(1,"step1");
		hist_Ratio->GetXaxis()->SetBinLabel(2,"step2");
		hist_Ratio->GetXaxis()->SetBinLabel(3,"step3");
		hist_Ratio->GetXaxis()->SetBinLabel(4,"step4");
		hist_Ratio->GetXaxis()->SetBinLabel(5,"step5");
	}

	ratiopad->cd();
	gPad->SetTopMargin(0);
	gPad->SetBottomMargin(0);
	ratiopad->SetGridy();
	hist_Ratio->SetMarkerStyle(20);
	hist_Ratio->SetMarkerSize(1.2);
	hist_Ratio->GetXaxis()->SetTitle(Variable);
	hist_Ratio->GetYaxis()->SetTitle("Data / MC");
	hist_Ratio->GetYaxis()->SetTitleSize(0.2);
	hist_Ratio->GetYaxis()->SetTitleOffset(0.23);
	hist_Ratio->GetYaxis()->SetLabelSize(0.14);
	hist_Ratio->GetYaxis()->CenterTitle();
	hist_Ratio->GetYaxis()->SetNdivisions(6);
	hist_Ratio->GetXaxis()->SetLabelSize(0.28);
	hist_Ratio->GetXaxis()->SetTitleSize(0.25);
	hist_Ratio->SetAxisRange(0.5,1.5,"y");
	hist_Ratio->Draw("e");

	c1->SaveAs(Save_dir+Variable+"_EventByStep.png");
}
