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
	gStyle->SetLabelSize(0.05, "YX");

	gStyle->SetTitleFont(42, "XY");
	gStyle->SetTitleXOffset(1);

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
	lt1.SetTextSize(0.055);

	TLatex lt2;
	lt2.SetTextAlign(12);
	lt2.SetNDC();
	lt2.SetTextFont(61);
	lt2.SetTextAngle(0);
	lt2.SetTextSize(0.073);

	TLatex lt3;
	lt3.SetTextAlign(12);
	lt3.SetNDC();
	lt3.SetTextAngle(0);
	lt3.SetTextFont(52);
	lt3.SetTextSize(0.06);

	TLatex lt4;
	lt4.SetTextAlign(32);
	lt4.SetNDC();
	lt4.SetTextAngle(0);
	lt4.SetTextFont(42);
	lt4.SetTextSize(0.06);

	TString Save_dir;
	Save_dir = "/afs/cern.ch/work/y/yjeong/catMacro/plots/";

	float x_1 = 0.12; //right top side x_1 = 0.73
	const float y_1 = 0.97; //right top side y_1 = 0.84
	float x_2 = x_1+0.08; //right top side y_2 = y_1-0.07
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

	const int nVariable = 5;
	const int nXbin = 6;
	const int nTree = 5;
	float x[nXbin] = {1,2,3,4,5,6};

	THStack *hs[nVariable];
	TLegend *l_[nVariable];
	TCanvas *c1[nVariable];
	TPad *plotpad_[nVariable];
	TPad *ratiopad_[nVariable];

	TH1F *hist_MC[nVariable];
	TH1F *hist_data[nVariable];
	TH1F *hist_[nVariable][nTree];
	TH1F *hist_Ratio[nVariable];

	TString Variable[nVariable] = {"","M(ll)","MET","NJet","NBJet"};

	cout<<"4"<<endl;
	float ev_sig[nVariable][nXbin] = {
		{221758+36466,221758+36466,157706+26337,160616+27243,125519+20458},//total
		{221758+36466,221758+36466,162305+26705,162305+26705,129035+20733},//dilep.M
		{218657+37030,218657+37030,160616+27243,160616+27243,127721+21118},//MET
		{215279+35959,215279+35959,157706+26337,157706+26337,125526+20459},//NJet
		{215271+35958,215271+35958,157697+26335,157697+26335,125519+20458}//NBJet
	};
	float ev_wjet[nVariable][nXbin] = {
		{11692,11692,1130,1154,0},
		{11692,11692,1237,1237,0},
		{11577,11577,1154,1154,0},
		{11420,11420,1130,1130,0},
		{11420,11420,1130,1130,0}
	};
	float ev_STop[nVariable][nXbin] = {
		{25389,25389,9574,10006,6711},
		{25389,25389,9973,9973,6984},
		{25123,25123,10006,10006,7022},
		{24587,24587,9574,9574,6711},
		{24587,24587,9574,9574,6711}
	};
	float ev_Diboson[nVariable][nXbin] = {
		{29945,29945,2186,2327,200},
		{29945,29945,2327,2327,213},
		{29221,29221,2327,2327,218},
		{28938,28938,2186,2186,200},
		{28938,28938,2186,2186,200}
	};
	float ev_DY[nVariable][nXbin] = {
		{80763,80763,5949,6100,623},
		{80763,80763,5966,5966,626},
		{80963,80963,6100,6100,647},
		{80720,80720,5949,5949,623},
		{80720,80720,5949,5949,623}
	};
	float ev_data[nVariable][nXbin] = {
		{374170,374170,189973,190833,139426},
		{374170,374170,190833,190833,142834},
		{371168,371168,189973,189973,142046},
		{366529,366529,186119,186119,139433},
		{366520,366520,186110,186110,139426}
	};//-*/

	double er_data = 0;
	for(int nData = 0; nData < nVariable; nData++){
		c1[nData] = new TCanvas;
		c1[nData]->RedrawAxis();

		l_[nData] = new TLegend(lx1,ly1,lx2,ly2);
		l_[nData]->SetFillColor(0);
		l_[nData]->SetLineColor(0);
		l_[nData]->SetLineStyle(kSolid);
		l_[nData]->SetLineWidth(1);
		l_[nData]->SetFillStyle(1001);
		l_[nData]->SetTextFont(42);
		l_[nData]->SetTextSize(0.045);

		plotpad_[nData] = new TPad(Form("plottitle_%d",nData),Form(""),0.02,0.3,0.98,0.98);
		ratiopad_[nData] = new TPad(Form("ratiotitle_%d",nData),Form(""),0.02,0.1,0.98,0.3);

		plotpad_[nData]->Draw();
		ratiopad_[nData]->Draw();
		plotpad_[nData]->cd();
		plotpad_[nData]->RedrawAxis();
		gPad->SetBottomMargin(0);

		cout<<"5"<<endl;
		for(int nT=0; nT < nTree; nT++){
			hist_[nT][nData] = new TH1F(Form("hist_%d%d",nT,nData),Form(""),nXbin-1,x);

			for(int i = 0; i < nXbin; i++){
				if(nT == 0){
					hist_[nT][nData]->SetBinContent(i+1, x[i]);
					hist_[nT][nData]->SetBinContent(i+1, ev_sig[nData][i]);
					hist_[nT][nData]->SetLineColor(2);
					hist_[nT][nData]->SetFillColor(2);
				}
				if(nT == 1){
					hist_[nT][nData]->SetBinContent(i+1, x[i]);
					hist_[nT][nData]->SetBinContent(i+1, ev_wjet[nData][i]);
					hist_[nT][nData]->SetLineColor(3);
					hist_[nT][nData]->SetFillColor(3);
				}
				if(nT == 2){
					hist_[nT][nData]->SetBinContent(i+1, x[i]);
					hist_[nT][nData]->SetBinContent(i+1, ev_STop[nData][i]);
					hist_[nT][nData]->SetLineColor(42);
					hist_[nT][nData]->SetFillColor(42);
				}
				if(nT == 3){
					hist_[nT][nData]->SetBinContent(i+1, x[i]);
					hist_[nT][nData]->SetBinContent(i+1, ev_Diboson[nData][i]);
					hist_[nT][nData]->SetLineColor(7);
					hist_[nT][nData]->SetFillColor(7);
				}
				if(nT == 4){
					hist_[nT][nData]->SetBinContent(i+1, x[i]);
					hist_[nT][nData]->SetBinContent(i+1, ev_DY[nData][i]);
					hist_[nT][nData]->SetLineColor(4);
					hist_[nT][nData]->SetFillColor(4);
				}
				if(nT == 5){
					hist_data[nData]->SetBinContent(i+1, x[i]);
					hist_data[nData]->SetBinContent(i+1, ev_data[nData][i]);
					hist_data[nData]->SetBinError(i+1,sqrt(er_data));
					hist_data[nData]->SetMarkerSize(1.1);
					hist_data[nData]->SetMarkerStyle(20);
				}
			}
		}
	}
	cout<<"6"<<endl;

	for(int nData = 0; nData < nVariable; nData++){
		hist_Ratio[nData] = new TH1F(Form("hist_Ratio%d",nData),Form(""),nXbin-1,x);
		hist_MC[nData] = new TH1F(Form("hist_MC%d",nData),Form(""),nXbin-1,x);
	}

	//hist->Scale(100/1);
	//hist->SetMinimum(0);
	for(int nData = 0; nData < nVariable; nData++){
		for(int nT = 0; nT < nTree; nT++){
			if(nT==0)l_[nData]->AddEntry(hist_[nT][nData],"t#bar{t} Powheg","f");
			if(nT==1)l_[nData]->AddEntry(hist_[nT][nData],"W+jet","f");
			if(nT==2)l_[nData]->AddEntry(hist_[nT][nData],"Single Top","f");
			if(nT==3)l_[nData]->AddEntry(hist_[nT][nData],"Diboson","f");
			if(nT==4)l_[nData]->AddEntry(hist_[nT][nData],"Z/#gamma^{*}#rightarrow#font[12]{l#lower[-0.4]{+}l#lower[-0.4]{#font[122]{\55}}}","f");
		}
		l_[nData]->AddEntry(hist_data[nData],"data","lp");
	}
	for(int nData = 0; nData < nVariable; nData++){
		hs[nData] = new THStack(Form("hs%d",nData),Form(""));
		for(int nT = 0; nT < nTree-1; nT++){
			hs[nData]->Add(hist_[nData][nT]);
		}
	}
	cout << "66" << endl;
	double ymax = 0;
	for(int nData = 0; nData < nVariable; nData++){
		ymax = hist_data[nData]->GetMaximum();
		hist_data[nData]->SetMaximum(ymax*1.3);
		hist_data[nData]->SetMinimum(0.1);
		hist_data[nData]->GetYaxis()->SetTitle("Number of Events");
		hist_data[nData]->GetYaxis()->SetTitleOffset(0.7);
		hist_data[nData]->GetYaxis()->SetTitleSize(0.08);
		hist_data[nData]->GetXaxis()->SetTitle(Variable[nData]);
		hist_data[nData]->Draw();
	}
	cout<<"7"<<endl;

	for(int nData = 0; nData < nVariable; nData++){
		hs[nData]->Draw("histsame");
		lt1.DrawLatex(xx_1,yy_1,"e^{#pm}#mu^{#mp}");
		lt2.DrawLatex(x_1,y_1,"CMS");
		lt3.DrawLatex(x_2,y_2,"Preliminary");
		lt4.DrawLatex(tx,ty,"35.9 fb^{-1}, #sqrt{s} = 13 TeV");
	}
	for(int nData = 0; nData < nVariable; nData++){
		hist_data[nData]->Draw("esame");
		l_[nData]->Draw();
	}
	cout<<"8"<<endl;
	for(int nData = 0; nData < nVariable; nData++){
		for(int nT = 0; nT < nTree; nT++){
			hist_MC[nData]->Add(hist_[nT][nData]);
		}
	}
	for(int nData = 0; nData < nVariable; nData++){
		hist_Ratio[nData]->Divide(hist_MC[nData],hist_data[Data],1,1,"b");
	}
	cout<<"9"<<endl;

	for(int nData = 0; nData < nVariable; nData++){
		if(nData == 0 ){
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(1,"Dilepton");
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(2,"Z veto");
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(3,">1 jets");
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(4,"MET");
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(5,">0 b jets");//-*/
		}
		if(nData >= 1){
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(1,"step1");
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(2,"step2");
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(3,"step3");
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(4,"step4");
			hist_Ratio[nData]->GetXaxis()->SetBinLabel(5,"step5");//-*/
		}
	}

	for(int nData = 0; nData < nVariable; nData++){
		ratiopad_[nData]->cd();
		ratiopad_[nData]->SetGridy();
		gPad->SetTopMargin(0);
		gPad->SetBottomMargin(0);

		hist_Ratio[nData]->GetXaxis()->SetTitle(Variable[nData]);
		hist_Ratio[nData]->SetMarkerStyle(20);
		hist_Ratio[nData]->SetMarkerSize(1.2);
		hist_Ratio[nData]->GetYaxis()->SetTitle("Data / MC");
		hist_Ratio[nData]->GetYaxis()->SetTitleSize(0.17);
		hist_Ratio[nData]->GetYaxis()->SetTitleOffset(0.3);
		hist_Ratio[nData]->GetXaxis()->SetTitleSize(0.23);
		hist_Ratio[nData]->GetYaxis()->SetLabelSize(0.13);
		hist_Ratio[nData]->GetXaxis()->SetLabelSize(0.25);
		hist_Ratio[nData]->GetYaxis()->CenterTitle();
		hist_Ratio[nData]->GetYaxis()->SetNdivisions(6);
		hist_Ratio[nData]->SetAxisRange(0.5,1.5,"y");
		hist_Ratio[nData]->Draw("e");

		c1[nData]->cd();
		c1[nData]->SaveAs(Save_dir+"EventByStep_"+Variable[nData]+".png");
	}
}
