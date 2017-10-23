{
	gROOT->SetStyle("Plain");//"Pub","Plain"
	gStyle->SetOptStat(0);//To display the mean and RMS: SetOptStat("mr"), nemruoi, ;
	gStyle->SetOptDate(0);//display date position
	/*gStyle->SetCanvasDefH(600);//Height of canvas
	  gStyle->SetCanvasDefW(600);//Width of canvas
	  gStyle->SetCanvasDefX(0);//POsition on screen
	  gStyle->SetCanvasDefY(0);*/

	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetPadTopMargin(0.07);
	gStyle->SetPadBottomMargin(0.1);
	gStyle->SetPadBorderMode(0);

	gStyle->SetLabelColor(1, "XYZ");
	gStyle->SetLabelFont(42, "XYZ");
	gStyle->SetLabelOffset(0.007, "XYZ");
	gStyle->SetLabelSize(0.05, "XYZ");

	gStyle->SetTitleColor(1, "XYZ");
	gStyle->SetTitleFont(42, "XYZ");
	gStyle->SetTitleSize(0.07, "X");
	gStyle->SetTitleSize(0.06, "Y");
	gStyle->SetTitleXOffset(1.1);
	gStyle->SetTitleYOffset(0.9);

	gStyle->SetAxisColor(1, "XYZ");
	gStyle->SetTickLength(0.03, "XYZ");
	gStyle->SetNdivisions(510, "XYZ");
	gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	gStyle->SetPadTickY(1);

	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(1);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);

	//---------------------------------------------------

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
	//----------------------------------------------------


	//-----------------------------------Coordinate of CMS Simulation------------------------------------------------------------
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

	const int StepNum = 1;//Step Num total:5
	const int nVariable = 2;//number of Variable 
	const int nChannel = 4;//total: 4 ---> Dilepton, MuEl, ElEl, MuMu.
	//int NJet[] = {4,5,6,7,8,9,10};
	//int NJet[] = {6};
	const int nMonteCal = 10;
	const int nRealData = 3;

	TH1F *histo_MonteCal[StepNum][nVariable][nChannel][nMonteCal];
	TH1F *histo_nRealData[StepNum][nVariable][nChannel][nRealData];

	TH1F *histo_RealData[StepNum][nVariable][nChannel];

	TH1F *histo_SingleTop[StepNum][nVariable][nChannel];
	TH1F *histo_Diboson[StepNum][nVariable][nChannel];
	TH1F *histo_Zr[StepNum][nVariable][nChannel];

	THStack *hs[StepNum][nVariable][nChannel];

	TH1F *histo_MC[StepNum][nVariable][nChannel];
	TH1F *histo_Data[StepNum][nVariable][nChannel];
	TH1F *histo_Ratio[StepNum][nVariable][nChannel];

	//-----------------------------------------------------------

	TCanvas *canv_[StepNum][nVariable][nChannel];
	TPad *plotpad_[StepNum][nVariable][nChannel];
	TPad *ratiopad_[StepNum][nVariable][nChannel];
	TLegend *l_[StepNum][nVariable][nChannel];

	TString PATH_samples;
	//PATH_samples = "/xrootd/store/user/yjeong/4TopFullHadronic/";//KISTI
	PATH_samples = "/xrootd/store/user/yjeong/TTBarDileptonAnalyzer/TtbarDileptonAnalyzer_";//KISTI
	//PATH_samples = "/xrootd/store/user/yjeong/TtBarDileptonAnalyzer/TtBarDileptonAnalyzer_";//KISTI
	//PATH_samples = "/cms/scratch/yjeong/";//KISTI
	TString Save_dir;
	Save_dir = "/cms/scratch/yjeong/catMacro/plots/";

	TString Variable[] = {"dilep.M()","nbjet"};//==================================variable

	TString Step_Cut[] = {"step>=1","step>=2","step>=3","step>=4","step>=5","step>=6"};

	TString TCut_base;
	TCut_base = "&&tri!=0&&filtered==1&&is3lep==2";
	TString weight_cut;
	weight_cut = "*genweight";//check, reduced wjet,z-gamma.
	//weight_cut = "";
	//weight_cut = "*(genweight*puweight*mueffweight*eleffweight*tri)";
	//weight_cut = "*genweight*mueffweight";
	//weight_cut = "*genweight*tri";
	//weight_cut = "*mueffweight*eleffweight";
	//weight_cut = "*tri";//worst
	//weight_cut = "*puweight";//worst
	//weight_cut = "*eleffweight";//not bad
	//weight_cut = "*mueffweight";//bad
	//weight_cut = "*genweight*eleffweight";

	TString Advanced_cut[] = {"","","","","&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20","&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20"};

	TString tt_others[] = {"&&!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==channel && gen_partonMode!=0)","&&!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)","&&!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)","&&!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)"};//channel = 0, 1, 2, 3 -> Dileoton, MuEl, ElEl, MuMu
	TString tt_signal[] = {"&&(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==channel && gen_partonMode!=0)","&& (gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)","&&(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)","&&(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)"};//channel = 0, 1, 2, 3 -> Dileoton, MuEl, ElEl, MuMu

	TString Step_txt[] = {"step1","step2","step3","step4","step5","step6"};

	TString Ytitle[] = {"Events / 5 GeV","Events"};//=====================================variable
	TString Xtitle[] = {"M(ll) [GeV]","b Jet Multiplicity"};//========================================variable

	TString Channel_Cut[] = {"&&channel","&&channel==1","&&channel==2","&&channel==3"};//Dilepton,MuEl,ElEl,MuMu;
	TString Channel_txt[] = {"Dilepton","MuEl","ElEl","MuMu"};

	////////////////////////////////Get Samples/////////////////////////////////

	const int Sample_Num = 13;//=======================================check
	TString Sample_name[Sample_Num] = {"TT_powheg","TT_powheg","WJets","SingleTbar_tW","SingleTop_tW","ZZ","WW","WZ","DYJets","DYJets_10to50","MuonEG_Run2016","DoubleEG_Run2016","DoubleMuon_Run2016"};//===============================check

	TString Legend_Name[] = {"t#bar{t}-signal (visible)","t#bar{t}-others","W+Jets","Single Top","Single Top","Diboson","Diboson","Diboson","Z/#gamma^{*}#rightarrow#font[12]{l#lower[-0.4]{+}l#lower[-0.4]{#font[122]{\55}}}"};//===============================check

	TFile *tfile[Sample_Num];

	for(int i = 0; i < Sample_Num; i++){
		tfile[i] = new TFile(PATH_samples+Sample_name[i]+".root");
	}

	TTree *tree[Sample_Num];
	TH1D *hnevents[Sample_Num];
	double totevents[Sample_Num];
	for(int i = 0; i < Sample_Num; i++){
		tree[i] = (TTree*)tfile[i]->Get("cattree/nom");
		hnevents[i] = (TH1D*)tfile[i]->Get("cattree/nevents");
		totevents[i] = hnevents[i]->Integral();
		if(i!=1&&i<=10)cout<<Sample_name[i]<<": "<<totevents[i]<<endl;//except tt-others
		if(i>10)cout<<Sample_name[i]<<": "<<++totevents[i]<<endl;
	}
	/////////////////////////////////////////////////////////////////////////////

	for(int nCh = 0; nCh < nChannel; nCh++){
		for(int NVar = 0; NVar < nVariable; NVar++){
			for(int NStep = 0; NStep < StepNum; NStep++){
				float nbin[] = {60,6};//===================================variable
				float xmin[] = {20,0};//====================================variable
				float xmax[] = {320,6};//====================================variable
				float size = 0.8;
				int ttsignal_c = 2;
				int ttothers_c = 906;
				int wjets_c = 3;
				int STop_c = 42;
				int Diboson_c = 7;
				int Z_pshy_c = 4;
				int data_c = 1;

				canv_[NVar][NStep][nCh] = new TCanvas(Form("Canv_%d_%d_%d",NVar,NStep,nCh),Form(""),800,800);
				//if(NVar>0)canv_[NVar][NStep][nCh]->SetLogy();
				canv_[NVar][NStep][nCh]->RedrawAxis();
				//canv_[NVar][NStep][nCh]->GetFrame()->Draw();

				l_[NVar][NStep][nCh] = new TLegend(lx1,ly1,lx2,ly2);
				l_[NVar][NStep][nCh]->SetFillColor(0);
				l_[NVar][NStep][nCh]->SetLineColor(0);
				l_[NVar][NStep][nCh]->SetLineStyle(kSolid);
				l_[NVar][NStep][nCh]->SetLineWidth(1);
				l_[NVar][NStep][nCh]->SetFillStyle(1001);
				l_[NVar][NStep][nCh]->SetTextFont(42);
				l_[NVar][NStep][nCh]->SetTextSize(0.035);

				plotpad_[NVar][NStep][nCh] = new TPad(Form("title_%d_%d_%d",NVar,NStep,nCh),Form(""),0.02,0.3,0.98,0.98);//x1,y1,x2,y2
				ratiopad_[NVar][NStep][nCh] = new TPad(Form("ratiotitle_%d_%d_%d",NVar,NStep,nCh),Form(""),0.02,0.1,0.98,0.3);

				plotpad_[NVar][NStep][nCh]->Draw();
				ratiopad_[NVar][NStep][nCh]->Draw();
				plotpad_[NVar][NStep][nCh]->cd();

				plotpad_[NVar][NStep][nCh]->SetLogy();
				plotpad_[NVar][NStep][nCh]->RedrawAxis();

				gPad->SetBottomMargin(0);

				histo_SingleTop[NVar][NStep][nCh] = new TH1F(Form("histo_SingleTop_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				histo_Diboson[NVar][NStep][nCh] = new TH1F(Form("histo_Diboson_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				histo_Zr[NVar][NStep][nCh] = new TH1F(Form("histo_Zr_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);

				/////////////////////////////////////////////MonteCals////////////////////////////////////////////////////////
				for(int nMC = 0; nMC < nMonteCal; nMC++){
					if(nMC==0){//tt-signal (visible)
						histo_MonteCal[NVar][NStep][nCh][nMC] = new TH1F(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
						tree[nMC]->Project(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Variable[NVar],Step_Cut[NStep]+Channel_Cut[nCh]+TCut_base+weight_cut+tt_signal[nCh]+Advanced_cut[NStep]);
					}
					//cout<<"1"<<endl;
					if(nMC==1){//tt-others
						histo_MonteCal[NVar][NStep][nCh][nMC] = new TH1F(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
						if(nCh==0)tree[nMC]->Project(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Variable[NVar],Step_Cut[NStep]+Channel_Cut[nCh]+TCut_base+weight_cut+tt_others[nCh]+Advanced_cut[NStep]);
						if(nCh!=0)tree[nMC]->Project(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Variable[NVar],Step_Cut[NStep]+Channel_Cut[nCh]+TCut_base+weight_cut+tt_others[nCh]+Advanced_cut[NStep]+Form("&& gen_partonMode==%d",nCh));
					}
					//cout<<"2"<<endl;
					if(nMC>1){//etc (except tt-powheg)
						histo_MonteCal[NVar][NStep][nCh][nMC] = new TH1F(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
						tree[nMC]->Project(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Variable[NVar],Step_Cut[NStep]+Channel_Cut[nCh]+TCut_base+weight_cut+Advanced_cut[NStep]);
					}
					//cout<<"3"<<endl;

					//histo_MonteCal[NVar][NStep][nCh][nMC]->SetLineWidth(2);
					if(nMC == 0){//tt-signal(visible)
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetLineColor(ttsignal_c);
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetFillColor(ttsignal_c);
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetMarkerColor(ttsignal_c);
					}
					if(nMC == 1){//tt-others
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetLineColor(ttothers_c);
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetFillColor(ttothers_c);
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetMarkerColor(ttothers_c);
					}
					if(nMC == 2){//w+jets
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetLineColor(wjets_c);
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetFillColor(wjets_c);
						histo_MonteCal[NVar][NStep][nCh][nMC]->SetMarkerColor(wjets_c);
					}
				}

				//////////////////////////////////////////////RealData/////////////////////////////////////////////////////

				for(int nReal = 0; nReal < nRealData; nReal++){
					histo_nRealData[NVar][NStep][nCh][nReal] = new TH1F(Form("histo_nRealData_%d_%d_%d_%d",NVar,NStep,nCh,nReal),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
					tree[nReal+10]->Project(Form("histo_nRealData_%d_%d_%d_%d",NVar,NStep,nCh,nReal),Variable[NVar],Step_Cut[NStep]+Channel_Cut[nCh]+TCut_base);
				}
				///////////////////////////////////////////// x-section candidate ///////////////////////////////////////////

				//double BR = 0.6741;//theoritical value W->Hadron
				//4top->all hadrons = BR^4.
				//ttbar->all hadrons = BR^2.
				const double lumi = 35.9*1000;//pb-1
				//const double lumi = 2.22*1000;//pb-1
				cout<<""<<endl;
				cout<<"---------------------------------------"<<Channel_txt[nCh]<<", "<<Variable[NVar]<<", "<<Step_txt[NStep]<<"-------------------------------------"<<endl;

				cout<<"lumi : "<<lumi<<" pb-1"<<endl;
				cout<<""<<endl;
				cout<<""<<endl;

				double MonteCal_xsec[] = {831.76, 831.76, 61526.7, 35.85, 35.85, 16.523, 118.7, 47.13, 6025.2, 18610};//======================================check

				/////////////////////////////////////////////// MonteCals ///////////////////////////////////////////////////

				for(int nMC = 0; nMC < nMonteCal; nMC++){
					histo_MonteCal[NVar][NStep][nCh][nMC]->Scale(MonteCal_xsec[nMC]*lumi/totevents[nMC]);
				}

				double MonteCal_ev = 0;
				double SingleTop_ev = 0;
				double Diboson_ev = 0;
				double Zgamma_ev = 0;

				for(int nMC = 3; nMC < nMonteCal; nMC++){//singleTop, Diboson, Z-gamma
					if(nMC >= 3 && nMC <= 4){
						histo_SingleTop[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);
					}
					if(nMC >= 5 && nMC <= 7){
						histo_Diboson[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);
					}
					if( nMC == 8){
						histo_Zr[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);
					}
				}

				histo_SingleTop[NVar][NStep][nCh]->SetLineColor(STop_c);
				histo_SingleTop[NVar][NStep][nCh]->SetFillColor(STop_c);
				histo_SingleTop[NVar][NStep][nCh]->SetMarkerColor(STop_c);
				//histo_SingleTop[NVar][NStep][nCh]->SetLineWidth(2);
				histo_Diboson[NVar][NStep][nCh]->SetLineColor(Diboson_c);
				histo_Diboson[NVar][NStep][nCh]->SetFillColor(Diboson_c);
				histo_Diboson[NVar][NStep][nCh]->SetMarkerColor(Diboson_c);
				//histo_Diboson[NVar][NStep][nCh]->SetLineWidth(2);
				histo_Zr[NVar][NStep][nCh]->SetLineColor(Z_pshy_c);
				histo_Zr[NVar][NStep][nCh]->SetFillColor(Z_pshy_c);
				histo_Zr[NVar][NStep][nCh]->SetMarkerColor(Z_pshy_c);
				//histo_Zr[NVar][NStep][nCh]->SetLineWidth(2);

				for(int nMC = 0; nMC < nMonteCal; nMC++){
					if(nMC==8)l_[NVar][NStep][nCh]->AddEntry(histo_Zr[NVar][NStep][nCh],Legend_Name[nMC], "lp");
					if(nMC==6)l_[NVar][NStep][nCh]->AddEntry(histo_Diboson[NVar][NStep][nCh],Legend_Name[nMC], "lp");
					if(nMC==4)l_[NVar][NStep][nCh]->AddEntry(histo_SingleTop[NVar][NStep][nCh],Legend_Name[nMC], "lp");
					if(nMC==2)l_[NVar][NStep][nCh]->AddEntry(histo_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "lp");
					if(nMC==1)l_[NVar][NStep][nCh]->AddEntry(histo_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "lp");
					if(nMC==0)l_[NVar][NStep][nCh]->AddEntry(histo_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "lp");
				}
				//--------------------------------------------Print-----------------------------------------------

				int Int_MonteCal[] = {0,};
				int Int_SingleTop = 0;
				int Int_Diboson = 0;
				int Int_Zgamma = 0;
				int total = 0;
				int total_1 = 0;
				int total_2 = 0;
				int bkg = 0 ;

				for(int nMC = 0; nMC < nMonteCal; nMC++){
					if(nMC>=0 && nMC<=2){
						MonteCal_ev = histo_MonteCal[NVar][NStep][nCh][nMC]->GetBinContent(nbin[NVar]+1);
						Int_MonteCal[nMC] = histo_MonteCal[NVar][NStep][nCh][nMC]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" yield : "<<Int_MonteCal[nMC]<<", err : "<<sqrt(MonteCal_ev)<<endl;
					}
					if(nMC==4){
						SingleTop_ev = histo_SingleTop[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
						Int_SingleTop = histo_SingleTop[NVar][NStep][nCh]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" yield : "<<Int_SingleTop<<", err : "<<sqrt(SingleTop_ev)<<endl;
					}
					if(nMC==6){
						Diboson_ev = histo_Diboson[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
						Int_Diboson = histo_Diboson[NVar][NStep][nCh]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" yield : "<<Int_Diboson<<", err : "<<sqrt(Diboson_ev)<<endl;
					}
					if(nMC==8){
						Zgamma_ev = histo_Zr[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
						Int_Zgamma = histo_Zr[NVar][NStep][nCh]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" yield : "<<Int_Zgamma<<", err : "<<sqrt(Zgamma_ev) <<endl;
					}
				}

				////////////////////////////////////////////////// RealData ///////////////////////////////////////////////////

				histo_RealData[NVar][NStep][nCh] = new TH1F(Form("histo_RealData_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				histo_RealData[NVar][NStep][nCh]->SetLineColor(data_c);
				histo_RealData[NVar][NStep][nCh]->SetLineWidth(1);
				//histo_RealData[NVar][NStep][nCh]->SetLineStyle(2);

				//int data_err = histo_RealData[NVar][NStep][nCh]->GetEntries();
				//histo_RealData[NVar][NStep][nCh]->SetBinError(nbin[NVar]+1,sqrt(data_err));

				for(int nReal = 0; nReal < nRealData; nReal++){
					histo_RealData[NVar][NStep][nCh]->Add(histo_nRealData[NVar][NStep][nCh][nReal]);
					histo_RealData[NVar][NStep][nCh]->Scale(1,"width");//binNormalize
				}

				hs[NVar][NStep][nCh] = new THStack(Form("hs_%d_%d_%d",NVar,NStep,nCh),Form(""));
				for(int nMC = 0; nMC < nMonteCal; nMC++){
					if(nMC >= 0 && nMC <= 2)hs[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);//MC
				}
				hs[NVar][NStep][nCh]->Add(histo_SingleTop[NVar][NStep][nCh]);
				hs[NVar][NStep][nCh]->Add(histo_Diboson[NVar][NStep][nCh]);
				hs[NVar][NStep][nCh]->Add(histo_Zr[NVar][NStep][nCh]);

				hs[NVar][NStep][nCh]->Add(histo_RealData[NVar][NStep][nCh]);//data
				double revents = 0;
				revents += histo_RealData[NVar][NStep][nCh]->GetEntries();

				double ev = 0;
				ev = histo_RealData[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
				histo_RealData[NVar][NStep][nCh]->SetBinError(nbin[NVar],sqrt(ev));
				cout<<"data : "<<revents<<endl;

				for(int nMC = 0; nMC < nMonteCal-7; nMC++){
					total_1 += Int_MonteCal[nMC];
					total_2 = Int_SingleTop+Int_Diboson+Int_Zgamma;
					total = total_1+total_2;
					bkg = total-Int_MonteCal[0];
				}
				cout<<""<<endl;
				cout<<"bkg : "<<bkg<<endl;
				cout<<"total : "<<total<<endl;
				cout<<""<<endl;

				histo_RealData[NVar][NStep][nCh]->SetLineColor(1);
				histo_RealData[NVar][NStep][nCh]->SetLineWidth(1);
				histo_RealData[NVar][NStep][nCh]->SetFillStyle(3001);
				histo_RealData[NVar][NStep][nCh]->SetFillColor(14);
				histo_RealData[NVar][NStep][nCh]->SetMarkerStyle(20);
				histo_RealData[NVar][NStep][nCh]->SetMarkerSize(1.2);

				l_[NVar][NStep][nCh]->AddEntry(histo_RealData[NVar][NStep][nCh],"Data ", "lp");

				double ymax = 0;
				ymax = hs[NVar][NStep][nCh]->GetMaximum();
				hs[NVar][NStep][nCh]->SetMaximum(ymax*100);
				hs[NVar][NStep][nCh]->SetMinimum(50);
				hs[NVar][NStep][nCh]->Draw();
				hs[NVar][NStep][nCh]->GetYaxis()->SetTitle(Ytitle[NVar]);
				canv_[NVar][NStep][nCh]->Modified();

				lt1.DrawLatex(xx_1,yy_1,Channel_txt[nCh]+"_"+Step_txt[NStep]);
				lt2.DrawLatex(x_1,y_1,"CMS");
				lt3.DrawLatex(x_2,y_2,"Preliminary");
				lt4.DrawLatex(tx,ty,"35.9 fb^{-1}, #sqrt{s} = 13 TeV");
				l_[NVar][NStep][nCh]->Draw();

				histo_MC[NVar][NStep][nCh] = new TH1F(Form("histo_MC_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				for(int nMC = 0; nMC < nMonteCal; nMC++){
					histo_MC[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);
				}

				histo_Ratio[NVar][NStep][nCh] = new TH1F(Form("histo_Ratio_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);

				histo_Ratio[NVar][NStep][nCh]->Divide(histo_RealData[NVar][NStep][nCh],histo_MC[NVar][NStep][nCh],1,1,"b");

				ratiopad_[NVar][NStep][nCh]->cd();
				gPad->SetTopMargin(0);
				gPad->SetBottomMargin(0);
				ratiopad_[NVar][NStep][nCh]->SetGridy();
				histo_Ratio[NVar][NStep][nCh]->SetMarkerStyle(20);
				histo_Ratio[NVar][NStep][nCh]->SetMarkerSize(1.2);
				histo_Ratio[NVar][NStep][nCh]->GetXaxis()->SetTitle(Xtitle[NVar]);
				histo_Ratio[NVar][NStep][nCh]->GetYaxis()->SetTitle("Data / MC");
				histo_Ratio[NVar][NStep][nCh]->GetYaxis()->SetLabelSize(0.11);
				histo_Ratio[NVar][NStep][nCh]->GetXaxis()->SetLabelSize(0.13);
				histo_Ratio[NVar][NStep][nCh]->GetXaxis()->SetTitleSize(0.16);
				//histo_Ratio[NVar][NStep][nCh]->GetYaxis()->SetTitleSize(0.16);
				//histo_Ratio[NVar][NStep][nCh]->SetAxisRange(0.8,1.2,"y");
				histo_Ratio[NVar][NStep][nCh]->Draw("e");

				/*auto rp = new TRatioPlot(histo_MC[NVar][NStep][nCh],histo_RealData[NVar][NStep][nCh]);
				  rp->Draw();
				  canv_[NVar][NStep][nCh]->Update();*/
				canv_[NVar][NStep][nCh]->cd();
				canv_[NVar][NStep][nCh]->SaveAs(Save_dir+"TtbarDileptonAnalyzer_"+Variable[NVar]+"_"+Channel_txt[nCh]+"_"+Step_txt[NStep]+".png");
			}
		}
	}
	cout<<"13TeV"<<endl;
}
