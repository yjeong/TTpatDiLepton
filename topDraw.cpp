{
	gROOT->SetStyle("Plain");//"Pub","Plain"
	//gROOT->ProcessLine("#include <vector>");
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

	/*gStyle->SetFrameBorderMode(0);
	  gStyle->SetFrameBorderSize(1);
	  gStyle->SetCanvasBorderMode(0);
	  gStyle->SetGridStyle(3);
	  gStyle->SetGridWidth(1);*/

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

	const int StepNum = 5;//Step Num total:5
	const int nVariable = 1;//number of Variable 
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
	TH1F *histo_Ratio[StepNum][nVariable][nChannel];

	//----------------------------PUEventReweighting------------------------------

	TH1F *histo_nReweight_MonteCal[StepNum][nVariable][nChannel][nMonteCal];
	TH1F *histo_nReweight_MonteCal_gen[StepNum][nVariable][nChannel][nMonteCal];
	TH1F *histo_nReweight_SingleTop[StepNum][nVariable][nChannel];
	TH1F *histo_nReweight_Diboson[StepNum][nVariable][nChannel];
	TH1F *histo_nReweight_Zr[StepNum][nVariable][nChannel];
	TH1F *histo_nReweight_MC[StepNum][nVariable][nChannel];

	TH1F *histo_nReweight_Data[StepNum][nVariable][nChannel];

	//-----------------------------------------------------------

	TCanvas *canv_[StepNum][nVariable][nChannel];
	TPad *plotpad_[StepNum][nVariable][nChannel];
	TPad *ratiopad_[StepNum][nVariable][nChannel];
	TLegend *l_[StepNum][nVariable][nChannel];

	TString PATH_samples;
	//PATH_samples = "/xrootd/store/user/yjeong/4TopFullHadronic/";//KISTI
	PATH_samples = "/xrootd/store/user/yjeong/TTBarDileptonAnalyzer/TtbarDileptonAnalyzer_";//KISTI
	//PATH_samples = "/xrootd/store/user/yjeong/TtBarDileptonAnalyzer/TtBarDileptonAnalyzer_";//KISTI
	//PATH_samples = "/xrootd/store/user/dhkim/v806_data_sep_v2/TtbarDileptonAnalyzer_";//KISTI
	//PATH_samples = "/cms/scratch/yjeong/";//KISTI

	TString Save_dir;
	Save_dir = "/cms/scratch/yjeong/catMacro/plots/";

	//TString Variable[nVariable] = {"nvertex","dilep.M()","njet","nbjet","pseudottbar.M()"};//==================================variable
	TString Variable[nVariable] = {"nvertex"};//==================================variable

	/*TString Var_int[] = {"nvertex"};
	  TString Var_float[] = {"met"};
	  int Var_int_size = sizeof(Var_int)/sizeof(Var_int[0]);
	  int Var_float_size = sizeof(Var_float)/sizeof(Var_float[0]);
	  int var_int[100][100];
	  float var_float[100][100];*/

	double single_cut_var[nVariable][nMonteCal]={0,};
	int step, is3lep;
	bool filtered;
	int nvertex, njet, nbjet, event;
	float met, tri, genweight, puweight, mueffweight, eleffweight, btagweight, topPtWeight, weight;

	int gen_partonChannel, gen_partonMode, gen_pseudoChannel, channel;

	TLorentzVector* dilep = NULL;
	TLorentzVector* pseudottbar = NULL;
	TLorentzVector* pseudojet1 = NULL;
	TLorentzVector* pseudojet2 = NULL;
	TLorentzVector* lep1 = NULL;
	TLorentzVector* lep2 = NULL;

	TString Step_Cut[StepNum] = {"step>=1","step>=2","step>=3","step>=4","step>=5"};

	TString TCut_base;
	TString weight_cut;
	TCut_base = "&&tri!=0&&filtered==1&&is3lep==2";
	weight_cut = "*genweight";//check, reduced wjet,z-gamma.

	TString Advanced_cut[StepNum] = {"","","","","&&pseudojet1.Pt()>30&&pseudojet2.Pt()>30&&lep1.Pt()>20&&lep2.Pt()>20"};

	TString tt_signal[nChannel] = {"&&(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==channel && gen_partonMode!=0)","&&(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)","&&(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)","&&(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)"};//channel = 0, 1, 2, 3 -> Dileoton, MuEl, ElEl, MuMu
	TString tt_others[nChannel] = {"&&!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==channel && gen_partonMode!=0)","&&!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)","&&!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)","&&!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)"};//channel = 0, 1, 2, 3 -> Dileoton, MuEl, ElEl, MuMu

	TString Step_txt[StepNum] = {"step1","step2","step3","step4","step5"};

	//TString Ytitle[nVariable] = {"Number of Events","Events / 5 GeV","Events","Events","Events / 90 GeV"};//=====================================variable
	//TString Xtitle[nVariable] = {"Number of good vertices","M(ll) [GeV]","Jet Multiplicity","b Jet Multiplicity","M^{t#tbar{t}}"};//========================================variable
	TString Ytitle[nVariable] = {"Number of Events"};//=====================================variable
	TString Xtitle[nVariable] = {"Number of good vertices"};//========================================variable

	TString Channel_Cut[nChannel] = {"&&channel==1","&&channel==2","&&channel==3","&&(channel==1 || channel == 2 || channel == 3)"};//Dilepton,MuEl,ElEl,MuMu;
	TString Channel_txt[nChannel] = {"MuEl","ElEl","MuMu","Dilepton"};

	////////////////////////////////Get Samples/////////////////////////////////

	const int Sample_Num = 13;//=======================================check
	TString Sample_name[Sample_Num] = {"TT_powheg","TT_powheg","WJets","SingleTbar_tW","SingleTop_tW","ZZ","WW","WZ","DYJets","DYJets_10to50","MuonEG_Run2016","DoubleEG_Run2016","DoubleMuon_Run2016"};//===============================check

	TString Legend_Name[Sample_Num] = {"t#bar{t}-signal (visible)","t#bar{t}-others","W+Jets","Single Top","Single Top","Diboson","Diboson","Diboson","Z/#gamma^{*}#rightarrow#font[12]{l#lower[-0.4]{+}l#lower[-0.4]{#font[122]{\55}}}"};//===============================check

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
		//for(int l1 = 0; l1 < Var_int_size; l1++) tree[i]->SetBranchAddress(Var_int[l1],var_int[l1]);
		//for(int l1 = 0; l1 < Var_float_size; l1++) tree[i]->SetBranchAddress(Var_float[l1],var_float[l1]);
		tree[i]->SetBranchAddress("dilep",&dilep);
		tree[i]->SetBranchAddress("pseudottbar",&pseudottbar);
		tree[i]->SetBranchAddress("pseudojet1",&pseudojet1);
		tree[i]->SetBranchAddress("pseudojet2",&pseudojet2);
		tree[i]->SetBranchAddress("lep1",&lep1);
		tree[i]->SetBranchAddress("lep2",&lep2);
		tree[i]->SetBranchAddress("event",&event);
		tree[i]->SetBranchAddress("nvertex",&nvertex);
		tree[i]->SetBranchAddress("njet",&njet);
		tree[i]->SetBranchAddress("met",&met);
		tree[i]->SetBranchAddress("nbjet",&nbjet);
		tree[i]->SetBranchAddress("step",&step);
		tree[i]->SetBranchAddress("tri",&tri);
		tree[i]->SetBranchAddress("filtered",&filtered);
		tree[i]->SetBranchAddress("is3lep",&is3lep);
		tree[i]->SetBranchAddress("genweight",&genweight);
		tree[i]->SetBranchAddress("puweight",&puweight);
		tree[i]->SetBranchAddress("eleffweight",&eleffweight);
		tree[i]->SetBranchAddress("mueffweight",&mueffweight);
		tree[i]->SetBranchAddress("btagweight",&btagweight);
		tree[i]->SetBranchAddress("topPtWeight",&topPtWeight);
		tree[i]->SetBranchAddress("weight",&weight);
		tree[i]->SetBranchAddress("channel",&channel);
		tree[i]->SetBranchAddress("gen_partonChannel",&gen_partonChannel);
		tree[i]->SetBranchAddress("gen_pseudoChannel",&gen_pseudoChannel);
		tree[i]->SetBranchAddress("gen_partonMode",&gen_partonMode);
	}

	/////////////////////////////////////////////////////////////////////////////

	/*int nbin[nVariable] = {70,60,10,6,10};//===================================variable
	  int xmin[nVariable] = {0,20,0,0,300};//====================================variable
	  int xmax[nVariable] = {70,320,10,6,1200};//====================================variable
	  int ymin[nVariable] = {300,300,100,100,1100};//====================================variable
	  */
	int nbin[nVariable] = {50};//===================================variable
	int xmin[nVariable] = {0};//====================================variable
	int xmax[nVariable] = {50};//====================================variable
	int ymin[nVariable] = {10};//====================================variable

	double MonteCal_xsec[nMonteCal] = {831.76, 831.76, 61526.7, 35.85, 35.85, 16.523, 118.7, 47.13, 6025.2, 18610};//======================================check
	double Reweight_nvertex_ch0_s1[70] = {1,1,1.11669,1.14979,1.53709,1.51114,1.48553,1.46179,1.46213,1.40139,1.40683,1.39015,1.36014,1.33894,1.25762,1.25963,1.19594,1.12766,1.1163,1.05614,0.983319,0.958495,0.90247,0.886433,0.834683,0.750569,0.729631,0.69422,0.646552,0.562234,0.521072,0.446841,0.418525,0.334112,0.367731,0.309285,0.27657,0.208889,0.183123,0.204334,0.176897,0.176332,0.14754,0.153501,0.0957998,0.0699733,0.140008,0.0724668,0.119624,0.131595,0.0720464,0.102486,0.148142,0.0325122,0,0.0907374,0.119604,0,0,0.101252,0.378878,0.136711,0.180665,0,0,0,0.338323,0,0,0};

	double Reweight_nvertex_ch0_s2[70] = {1,1,1.11669,1.14979,1.53709,1.51114,1.48553,1.46179,1.46213,1.40139,1.40683,1.39015,1.36014,1.33894,1.25762,1.25963,1.19594,1.12766,1.1163,1.05614,0.983319,0.958495,0.90247,0.886433,0.834683,0.750569,0.729631,0.69422,0.646552,0.562234,0.521072,0.446841,0.418525,0.334112,0.367731,0.309285,0.27657,0.208889,0.183123,0.204334,0.176897,0.176332,0.14754,0.153501,0.0957998,0.0699733,0.140008,0.0724668,0.119624,0.131595,0.0720464,0.102486,0.148142,0.0325122,0,0.0907374,0.119604,0,0,0.101252,0.378878,0.136711,0.180665,0,0,0,0.338323,0,0,0};

	double Reweight_nvertex_ch0_s3[70] = {1,1,0.703465,1.40019,1.49645,1.47803,1.5314,1.55007,1.592,1.53754,1.44048,1.38773,1.37172,1.35659,1.29884,1.24348,1.19719,1.15386,1.07685,1.03385,0.964184,0.946938,0.910429,0.879424,0.796209,0.770304,0.749509,0.675919,0.638597,0.566006,0.522009,0.43643,0.418617,0.360905,0.357568,0.31299,0.27697,0.242652,0.196638,0.213475,0.1807,0.193461,0.159858,0.129387,0.128425,0.100497,0.113911,0.0718695,0.0825647,0.0562259,0.0498585,0.0895106,0.122015,0.139937,0,0.0951575,0.0961218,0,0,0.167064,0.138257,0.158111,0.192406,0,0,0,0.369809,0,0,0};

	double Reweight_nvertex_ch0_s4[70] = {1,1,0.703465,1.40019,1.49645,1.47803,1.5314,1.55007,1.592,1.53754,1.44048,1.38773,1.37172,1.35659,1.29884,1.24348,1.19719,1.15386,1.07685,1.03385,0.964184,0.946938,0.910429,0.879424,0.796209,0.770304,0.749509,0.675919,0.638597,0.566006,0.522009,0.43643,0.418617,0.360905,0.357568,0.31299,0.27697,0.242652,0.196638,0.213475,0.1807,0.193461,0.159858,0.129387,0.128425,0.100497,0.113911,0.0718695,0.0825647,0.0562259,0.0498585,0.0895106,0.122015,0.139937,0,0.0951575,0.0961218,0,0,0.167064,0.138257,0.158111,0.192406,0,0,0,0.369809,0,0,0};

	double Reweight_nvertex_ch0_s5[70] = {1,1,0.558507,1.39006,1.37191,1.47752,1.48753,1.61026,1.55835,1.53473,1.44318,1.38293,1.35698,1.35568,1.28899,1.25199,1.18833,1.1369,1.07694,1.02453,0.989556,0.932239,0.906399,0.874475,0.797277,0.776935,0.735366,0.676318,0.623836,0.58636,0.522878,0.47006,0.428486,0.379308,0.369805,0.319417,0.27797,0.238575,0.194423,0.210718,0.183796,0.201018,0.154176,0.145329,0.1408,0.0872786,0.107129,0.0756727,0.0951131,0.0545657,0.0804836,0.100963,0.109478,0.156926,0,0.107898,0.103514,0,0,0.100891,0,0.166641,0.212635,0,0,0,0.348525,0,0,0};

	for(int nCh = 0; nCh < nChannel; nCh++){
		for(int NVar = 0; NVar < nVariable; NVar++){
			for(int NStep = 0; NStep < StepNum; NStep++){
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

				//plotpad_[NVar][NStep][nCh]->SetLogy();
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

					if(nMC==1){//tt-others
						histo_MonteCal[NVar][NStep][nCh][nMC] = new TH1F(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
						if(nCh==0)tree[nMC]->Project(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Variable[NVar],Step_Cut[NStep]+Channel_Cut[nCh]+TCut_base+weight_cut+tt_others[nCh]+Advanced_cut[NStep]);
						if(nCh!=0)tree[nMC]->Project(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Variable[NVar],Step_Cut[NStep]+Channel_Cut[nCh]+TCut_base+weight_cut+tt_others[nCh]+Advanced_cut[NStep]+Form("&& gen_partonMode==%d",nCh));
					}

					if(nMC>1){//etc (except tt-powheg)
						histo_MonteCal[NVar][NStep][nCh][nMC] = new TH1F(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
						tree[nMC]->Project(Form("histo_MonteCal_%d_%d_%d_%d",NVar,NStep,nCh,nMC),Variable[NVar],Step_Cut[NStep]+Channel_Cut[nCh]+TCut_base+weight_cut+Advanced_cut[NStep]);
					}

					//histo_MonteCal[NVar][NStep][nCh][nMC]->SetLineWidth(2);
					/*if(nMC == 0){//tt-signal(visible)
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
					}*/
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

				/////////////////////////////////////////////// MonteCals ///////////////////////////////////////////////////

				for(int nMC = 0; nMC < nMonteCal; nMC++){
					histo_MonteCal[NVar][NStep][nCh][nMC]->Scale(MonteCal_xsec[nMC]*lumi/totevents[nMC]);
				}

				double MonteCal_ev = 0;
				double SingleTop_ev = 0;
				double Diboson_ev = 0;
				double Zgamma_ev = 0;

				for(int nMC = 0; nMC < nMonteCal; nMC++){//singleTop, Diboson, Z-gamma
					if(nMC >= 3 && nMC <= 4){
						histo_SingleTop[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);
					}
					if(nMC >= 5 && nMC <= 7){
						histo_Diboson[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);
					}
					if( nMC >= 8 && nMC <= 9){
						histo_Zr[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);
					}
				}

				/*histo_SingleTop[NVar][NStep][nCh]->SetLineColor(STop_c);
				  histo_SingleTop[NVar][NStep][nCh]->SetFillColor(STop_c);
				  histo_SingleTop[NVar][NStep][nCh]->SetMarkerColor(STop_c);
				  histo_Diboson[NVar][NStep][nCh]->SetLineColor(Diboson_c);
				  histo_Diboson[NVar][NStep][nCh]->SetFillColor(Diboson_c);
				  histo_Diboson[NVar][NStep][nCh]->SetMarkerColor(Diboson_c);
				  histo_Zr[NVar][NStep][nCh]->SetLineColor(Z_pshy_c);
				  histo_Zr[NVar][NStep][nCh]->SetFillColor(Z_pshy_c);
				  histo_Zr[NVar][NStep][nCh]->SetMarkerColor(Z_pshy_c);*/

				/*for(int nMC = 0; nMC < nMonteCal; nMC++){
				  if(nMC==8)l_[NVar][NStep][nCh]->AddEntry(histo_Zr[NVar][NStep][nCh],Legend_Name[nMC], "lp");
				  if(nMC==6)l_[NVar][NStep][nCh]->AddEntry(histo_Diboson[NVar][NStep][nCh],Legend_Name[nMC], "lp");
				  if(nMC==4)l_[NVar][NStep][nCh]->AddEntry(histo_SingleTop[NVar][NStep][nCh],Legend_Name[nMC], "lp");
				  if(nMC==2)l_[NVar][NStep][nCh]->AddEntry(histo_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "lp");
				  if(nMC==1)l_[NVar][NStep][nCh]->AddEntry(histo_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "lp");
				  if(nMC==0)l_[NVar][NStep][nCh]->AddEntry(histo_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "lp");
				  }*/

				//--------------------------------------------Print-----------------------------------------------

				int Int_MonteCal[nMonteCal] = {0,};
				int Int_SingleTop = 0;
				int Int_Diboson = 0;
				int Int_Zgamma = 0;
				int total_1 = 0;
				int total_2 = 0;
				int total;
				int bkg;

				for(int nMC = 0; nMC < nMonteCal; nMC++){
					if(nMC<=2){
						MonteCal_ev = histo_MonteCal[NVar][NStep][nCh][nMC]->GetBinContent(nbin[NVar]+1);
						Int_MonteCal[nMC] = histo_MonteCal[NVar][NStep][nCh][nMC]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" yield : "<<Int_MonteCal[nMC]<<", err : "<<sqrt(MonteCal_ev)<<endl;
						total_1 += Int_MonteCal[nMC];
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

				total_2 = Int_SingleTop+Int_Diboson+Int_Zgamma;
				total = total_1+total_2;
				bkg = total-Int_MonteCal[0];

				cout<<""<<endl;
				cout<<"bkg : "<<bkg<<endl;
				cout<<"total : "<<total<<endl;
				cout<<""<<endl;

				////////////////////////////////////////////////// RealData ///////////////////////////////////////////////////

				histo_RealData[NVar][NStep][nCh] = new TH1F(Form("histo_RealData_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				histo_RealData[NVar][NStep][nCh]->SetLineColor(data_c);

				//histo_RealData[NVar][NStep][nCh]->SetLineStyle(2);

				for(int nReal = 0; nReal < nRealData; nReal++){
					histo_RealData[NVar][NStep][nCh]->Add(histo_nRealData[NVar][NStep][nCh][nReal]);
				}

				histo_MC[NVar][NStep][nCh] = new TH1F(Form("histo_MC_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				for(int nMC = 0; nMC < nMonteCal; nMC++){
					histo_MC[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);
				}

				histo_MC[NVar][NStep][nCh]->SetLineWidth(2);

				double nev_Data = 1;
				nev_Data = histo_RealData[NVar][NStep][nCh]->Integral(1,nbin[NVar]+1);
				histo_RealData[NVar][NStep][nCh]->Scale(1/nev_Data);

				double nev_MC = 1;
				nev_MC = histo_MC[NVar][NStep][nCh]->Integral(1,nbin[NVar]+1);
				histo_MC[NVar][NStep][nCh]->Scale(1/nev_MC);

				/*hs[NVar][NStep][nCh] = new THStack(Form("hs_%d_%d_%d",NVar,NStep,nCh),Form(""));
				  for(int nMC = 0; nMC < nMonteCal; nMC++){
				  if(nMC >= 0 && nMC <= 2)hs[NVar][NStep][nCh]->Add(histo_MonteCal[NVar][NStep][nCh][nMC]);//MC
				  }

				  hs[NVar][NStep][nCh]->Add(histo_SingleTop[NVar][NStep][nCh]);
				  hs[NVar][NStep][nCh]->Add(histo_Diboson[NVar][NStep][nCh]);
				  hs[NVar][NStep][nCh]->Add(histo_Zr[NVar][NStep][nCh]);*/

				double revents = 0;
				revents += histo_RealData[NVar][NStep][nCh]->GetEntries();

				double ev = 0;
				ev = histo_RealData[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
				histo_RealData[NVar][NStep][nCh]->SetBinError(nbin[NVar]+1,sqrt(ev));
				cout<<"data : "<<revents<<endl;
				cout<<""<<endl;
				cout<<""<<endl;

				histo_RealData[NVar][NStep][nCh]->SetLineColor(1);
				histo_RealData[NVar][NStep][nCh]->SetLineWidth(1);
				histo_RealData[NVar][NStep][nCh]->SetFillStyle(3001);
				histo_RealData[NVar][NStep][nCh]->SetFillColor(14);
				histo_RealData[NVar][NStep][nCh]->SetMarkerStyle(20);
				histo_RealData[NVar][NStep][nCh]->SetMarkerSize(1.2);

				l_[NVar][NStep][nCh]->AddEntry(histo_RealData[NVar][NStep][nCh],"Data ", "lp");

				/*double ymax = 0;
				  ymax = hs[NVar][NStep][nCh]->GetMaximum();
				  hs[NVar][NStep][nCh]->SetMaximum(ymax*100);
				  hs[NVar][NStep][nCh]->SetMinimum(ymin[NVar]);
				  hs[NVar][NStep][nCh]->Draw();

				  histo_RealData[NVar][NStep][nCh]->Draw("same");
				  canv_[NVar][NStep][nCh]->Modified();

				  lt1.DrawLatex(xx_1,yy_1,Channel_txt[nCh]+"_"+Step_txt[NStep]);
				  lt2.DrawLatex(x_1,y_1,"CMS");
				  lt3.DrawLatex(x_2,y_2,"Preliminary");
				  lt4.DrawLatex(tx,ty,"35.9 fb^{-1}, #sqrt{s} = 13 TeV");
				  l_[NVar][NStep][nCh]->Draw();*/

				const int n = histo_RealData[NVar][NStep][nCh]->GetNbinsX();//============================================>Variable
				double DataBin_ev[n];
				for(int i = 0; i < n; i++){
					DataBin_ev[i] = histo_RealData[NVar][NStep][nCh]->GetBinContent(i);
				}

				double MCBin_ev[n];
				for(int i = 0; i < n; i++){
					MCBin_ev[i] = histo_MC[NVar][NStep][nCh]->GetBinContent(i);
				}

				double RatioBin_ev[n];
				for(int i = 0; i < n; i++){
					RatioBin_ev[i] = DataBin_ev[i]/MCBin_ev[i];
					//cout<<"FixedRatio: "<<RatioBin_ev[i]<<",    Data: "<<DataBin_ev[i]<<",    MC: "<<MCBin_ev[i]<<endl;
					//cout<<RatioBin_ev[i]<<endl;
				}
				//cout<<""<<endl;
				//cout<<""<<endl;

				////////////////////////////////////////////PUeventReweighting////////////////////////////////////////

				histo_nReweight_SingleTop[NVar][NStep][nCh] = new TH1F(Form("histo_nReweight_SingleTop_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				histo_nReweight_Diboson[NVar][NStep][nCh] = new TH1F(Form("histo_nReweight_Diboson_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				histo_nReweight_Zr[NVar][NStep][nCh] = new TH1F(Form("histo_nReweight_Zr_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);

				for(int tr = 0; tr < nMonteCal; tr++){
					histo_nReweight_MonteCal[NVar][NStep][nCh][tr] = new TH1F(Form("histo_nReweight_%d_%d_%d_%d",NVar,NStep,nCh,tr),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
					cout<<"Reweighted tree event, "<<Sample_name[tr]<<": "<<tree[tr]->GetEntries()<<endl;
					for(int nev = 0; nev < tree[tr]->GetEntries(); nev++){

						if(tr==0 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==1 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==2 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==3 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==4 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==5 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==6 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==7 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==8 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}
						if(tr==9 && NVar==0) {tree[tr]->GetEntry(nev); single_cut_var[NVar][tr] = nvertex;}

						if(dilep == NULL ) continue;
						if(pseudottbar == NULL) continue;
						if(pseudojet1 == NULL) continue;
						if(pseudojet2 == NULL) continue;
						if(lep1 == NULL) continue;
						if(lep2 == NULL) continue;

						if(nCh==0) if(!(channel==1)) continue;
						if(nCh==1) if(!(channel==2)) continue;
						if(nCh==2) if(!(channel==3)) continue;
						if(nCh==3) if(!(channel==1 || channel==2 || channel==3)) continue;


						if(!(tri!=0&&filtered==1&&is3lep==2)) continue;
						if(tr==0&&nCh==0) if(!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==channel && gen_partonMode!=0)) continue;//tt-signal
						if(tr==0&&nCh!=0) if(!(gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel)) continue;
						if(tr==1&&nCh==0) if((gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==channel && gen_partonMode!=0)) continue;//tt-others
						if(tr==1&&nCh!=0) if((gen_partonChannel==2 && gen_partonMode==gen_pseudoChannel && gen_partonMode==channel)) continue;
						double PUeventReweight = 1;
						PUeventReweight = puweight*tri;

						/*if(nCh==0 && NStep==0) if(!(PUeventReweight = Reweight_nvertex_ch0_s1[nvertex-1])) continue;
						  if(nCh==0 && NStep==1) if(!(PUeventReweight = Reweight_nvertex_ch0_s2[nvertex-1])) continue;
						  if(nCh==0 && NStep==2) if(!(PUeventReweight = Reweight_nvertex_ch0_s3[nvertex-1])) continue;
						  if(nCh==0 && NStep==3) if(!(PUeventReweight = Reweight_nvertex_ch0_s4[nvertex-1])) continue;
						  if(nCh==0 && NStep==4) if(!(PUeventReweight = Reweight_nvertex_ch0_s5[nvertex-1])) continue;*/

						if(NStep==0 && step>=1)histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->Fill(single_cut_var[NVar][tr],PUeventReweight);
						if(NStep==1 && step>=2)histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->Fill(single_cut_var[NVar][tr],PUeventReweight);
						if(NStep==2 && step>=3)histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->Fill(single_cut_var[NVar][tr],PUeventReweight);
						if(NStep==3 && step>=4)histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->Fill(single_cut_var[NVar][tr],PUeventReweight);
						if(NStep==4 && step>=5)histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->Fill(single_cut_var[NVar][tr],PUeventReweight);
						//cout<<"single_cut_var: "<<single_cut_var[0]<<endl;
					}

					if(tr == 0){//tt-signal(visible)
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetLineColor(ttsignal_c);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetFillColor(ttsignal_c);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetFillStyle(1001);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetLineWidth(2);
					}
					if(tr == 1){//tt-others
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetLineColor(ttothers_c);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetFillColor(ttothers_c);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetFillStyle(1001);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetLineWidth(2);
					}
					if(tr == 2){//w+jets
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetLineColor(wjets_c);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetFillColor(wjets_c);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetFillStyle(1001);
						histo_nReweight_MonteCal[NVar][NStep][nCh][tr]->SetLineWidth(2);
					}
				}

				double Reweight_MonteCal_ev = 0;
				double Reweight_SingleTop_ev = 0;
				double Reweight_Diboson_ev = 0;
				double Reweight_Zgamma_ev = 0;
				int Reweight_Int_MonteCal[nMonteCal] = {0,};
				int Reweight_Int_SingleTop = 0;
				int Reweight_Int_Diboson = 0;
				int Reweight_Int_Zgamma = 0;

				for(int nMC = 0; nMC < nMonteCal; nMC++){
					histo_nReweight_MonteCal[NVar][NStep][nCh][nMC]->Scale(MonteCal_xsec[nMC]*lumi/totevents[nMC]);
				}

				for(int nMC = 3; nMC < nMonteCal; nMC++){//singleTop, Diboson, Z-gamma
					if(nMC >= 3 && nMC <= 4){
						histo_nReweight_SingleTop[NVar][NStep][nCh]->Add(histo_nReweight_MonteCal[NVar][NStep][nCh][nMC]);
					}

					if(nMC >= 5 && nMC <= 7){
						histo_nReweight_Diboson[NVar][NStep][nCh]->Add(histo_nReweight_MonteCal[NVar][NStep][nCh][nMC]);
					}

					if( nMC >= 8 && nMC <= 9){
						histo_nReweight_Zr[NVar][NStep][nCh]->Add(histo_nReweight_MonteCal[NVar][NStep][nCh][nMC]);
					}
				}

				cout<<""<<endl;
				for(int nMC = 0; nMC < nMonteCal; nMC++){
					if(nMC<=2){
						Reweight_MonteCal_ev = histo_nReweight_MonteCal[NVar][NStep][nCh][nMC]->GetBinContent(nbin[NVar]+1);
						Reweight_Int_MonteCal[nMC] = histo_nReweight_MonteCal[NVar][NStep][nCh][nMC]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" Reweighted yield: "<<Reweight_Int_MonteCal[nMC]<<", err : "<<sqrt(Reweight_MonteCal_ev)<<endl;
					}
					if(nMC==4){
						Reweight_SingleTop_ev = histo_nReweight_SingleTop[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
						Reweight_Int_SingleTop = histo_nReweight_SingleTop[NVar][NStep][nCh]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" Reweighted yield : "<<Reweight_Int_SingleTop<<", err : "<<sqrt(Reweight_SingleTop_ev)<<endl;
					}
					if(nMC==6){
						Reweight_Diboson_ev = histo_nReweight_Diboson[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
						Reweight_Int_Diboson = histo_nReweight_Diboson[NVar][NStep][nCh]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" Reweighted yield : "<<Reweight_Int_Diboson<<", err : "<<sqrt(Reweight_Diboson_ev)<<endl;
					}
					if(nMC==8){
						Reweight_Zgamma_ev = histo_nReweight_Zr[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
						Reweight_Int_Zgamma = histo_nReweight_Zr[NVar][NStep][nCh]->Integral(1,nbin[NVar]+1);
						cout<<Legend_Name[nMC]<<" Reweighted yield : "<<Reweight_Int_Zgamma<<", err : "<<sqrt(Reweight_Zgamma_ev)<<endl;
					}
				}

				histo_nReweight_SingleTop[NVar][NStep][nCh]->SetLineColor(STop_c);
				histo_nReweight_SingleTop[NVar][NStep][nCh]->SetFillColor(STop_c);
				histo_nReweight_SingleTop[NVar][NStep][nCh]->SetFillStyle(1001);
				histo_nReweight_SingleTop[NVar][NStep][nCh]->SetLineWidth(2);
				histo_nReweight_Diboson[NVar][NStep][nCh]->SetLineColor(Diboson_c);
				histo_nReweight_Diboson[NVar][NStep][nCh]->SetFillColor(Diboson_c);
				histo_nReweight_Diboson[NVar][NStep][nCh]->SetLineWidth(2);
				histo_nReweight_Diboson[NVar][NStep][nCh]->SetFillStyle(1001);
				histo_nReweight_Zr[NVar][NStep][nCh]->SetLineColor(Z_pshy_c);
				histo_nReweight_Zr[NVar][NStep][nCh]->SetFillColor(Z_pshy_c);
				histo_nReweight_Zr[NVar][NStep][nCh]->SetLineWidth(2);
				histo_nReweight_Zr[NVar][NStep][nCh]->SetFillStyle(1001);

				for(int nMC = 0; nMC < nMonteCal; nMC++){
					if(nMC==8)l_[NVar][NStep][nCh]->AddEntry(histo_nReweight_Zr[NVar][NStep][nCh],Legend_Name[nMC], "f");
					if(nMC==6)l_[NVar][NStep][nCh]->AddEntry(histo_nReweight_Diboson[NVar][NStep][nCh],Legend_Name[nMC], "f");
					if(nMC==4)l_[NVar][NStep][nCh]->AddEntry(histo_nReweight_SingleTop[NVar][NStep][nCh],Legend_Name[nMC], "f");
					if(nMC==2)l_[NVar][NStep][nCh]->AddEntry(histo_nReweight_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "f");
					if(nMC==1)l_[NVar][NStep][nCh]->AddEntry(histo_nReweight_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "f");
					if(nMC==0)l_[NVar][NStep][nCh]->AddEntry(histo_nReweight_MonteCal[NVar][NStep][nCh][nMC],Legend_Name[nMC], "f");
				}

				//---------------------------------------------------------

				histo_nReweight_Data[NVar][NStep][nCh] = new TH1F(Form("histo_nReweight_Data_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				histo_nReweight_Data[NVar][NStep][nCh]->SetLineColor(data_c);

				for(int nReal = 0; nReal < nRealData; nReal++){
					histo_nReweight_Data[NVar][NStep][nCh]->Add(histo_nRealData[NVar][NStep][nCh][nReal]);
				}

				/*double Reweight_SingleTop[n];
				  double Reweight_Diboson[n];
				  double Reweight_Zr[n];
				  double Reweight_ttsignal[n];
				  double Reweight_ttothers[n];
				  double Reweight_wjets[n];
				  double Reweight_totMC[n];
				  double Reweight_DataBin[n];

				  double Reweight_RatioBin[n];*/

				/*for(int i = 0; i < n; i++){
				  if(i==61){
				  Reweight_SingleTop[i] = histo_nReweight_SingleTop[NVar][NStep][nCh]->GetBinContent(i);
				  Reweight_Diboson[i] = histo_nReweight_Diboson[NVar][NStep][nCh]->GetBinContent(i);
				  Reweight_Zr[i] = histo_nReweight_Zr[NVar][NStep][nCh]->GetBinContent(i);
				  Reweight_ttsignal[i] = histo_nReweight_MonteCal[NVar][NStep][nCh][0]->GetBinContent(i);
				  Reweight_ttothers[i] = histo_nReweight_MonteCal[NVar][NStep][nCh][1]->GetBinContent(i);
				  Reweight_wjets[i] = histo_nReweight_MonteCal[NVar][NStep][nCh][2]->GetBinContent(i);
				  Reweight_DataBin[i] = histo_nReweight_Data[NVar][NStep][nCh]->GetBinContent(i);

				  Reweight_tot_MC[i] = Reweight_SingleTop[i]+Reweight_Diboson[i]+Reweight_Zr[i]+Reweight_ttsignal[i]+Reweight_ttothers[i]+Reweight_wjets[i];

				  Reweight_RatioBin[i] = Reweight_DataBin[i]/Reweight_totMC[i];

				  cout<<"ReweightBin: "<<Reweight_RatioBin_ev[i]<<",    ReweightData: "<<Reweight_DataBin_ev[i]<<",    MC: "<<Reweight_MCBin_ev[i]<<endl;
				  }
				  }*/

				double ev_Re = 0;
				ev_Re = histo_nReweight_Data[NVar][NStep][nCh]->GetBinContent(nbin[NVar]+1);
				histo_nReweight_Data[NVar][NStep][nCh]->SetBinError(nbin[NVar]+1,sqrt(ev_Re));
				cout<<"data : "<<revents<<endl;
				cout<<""<<endl;

				histo_nReweight_Data[NVar][NStep][nCh]->SetLineColor(1);
				histo_nReweight_Data[NVar][NStep][nCh]->SetLineWidth(1);
				histo_nReweight_Data[NVar][NStep][nCh]->SetMarkerStyle(20);
				histo_nReweight_Data[NVar][NStep][nCh]->SetMarkerSize(1.2);

				//--------------------------------------------------------

				plotpad_[NVar][NStep][nCh]->cd();
				double ymax = 0;
				ymax = histo_nReweight_Data[NVar][NStep][nCh]->GetMaximum();
				histo_nReweight_Data[NVar][NStep][nCh]->SetMaximum(ymax*1.3);
				histo_nReweight_Data[NVar][NStep][nCh]->GetYaxis()->SetTitle(Ytitle[NVar]);
				histo_nReweight_Data[NVar][NStep][nCh]->SetMinimum(ymin[NVar]);
				histo_nReweight_Data[NVar][NStep][nCh]->Draw();

				hs[NVar][NStep][nCh] = new THStack(Form("hs_%d_%d_%d",NVar,NStep,nCh),Form(""));

				for(int nMC = 0; nMC < nMonteCal-7; nMC++){
					hs[NVar][NStep][nCh]->Add(histo_nReweight_MonteCal[NVar][NStep][nCh][nMC]);//MC
				}

				hs[NVar][NStep][nCh]->Add(histo_nReweight_SingleTop[NVar][NStep][nCh]);
				hs[NVar][NStep][nCh]->Add(histo_nReweight_Diboson[NVar][NStep][nCh]);
				hs[NVar][NStep][nCh]->Add(histo_nReweight_Zr[NVar][NStep][nCh]);

				hs[NVar][NStep][nCh]->Draw("same");

				canv_[NVar][NStep][nCh]->Modified();

				lt1.DrawLatex(xx_1,yy_1,Channel_txt[nCh]+"_"+Step_txt[NStep]);
				lt2.DrawLatex(x_1,y_1,"CMS");
				lt3.DrawLatex(x_2,y_2,"Preliminary");
				lt4.DrawLatex(tx,ty,"35.9 fb^{-1}, #sqrt{s} = 13 TeV");
				l_[NVar][NStep][nCh]->Draw();

				/////////////////////////////////////////////// Ratio plot //////////////////////////////////////////

				histo_nReweight_MC[NVar][NStep][nCh] = new TH1F(Form("histo_nReweight_MC_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				for(int nMC = 0; nMC < nMonteCal; nMC++){
					histo_nReweight_MC[NVar][NStep][nCh]->Add(histo_nReweight_MonteCal[NVar][NStep][nCh][nMC]);
				}

				histo_nReweight_MC[NVar][NStep][nCh]->SetLineWidth(2);

				histo_Ratio[NVar][NStep][nCh] = new TH1F(Form("histo_Ratio_%d_%d_%d",NVar,NStep,nCh),Form(""),nbin[NVar],xmin[NVar],xmax[NVar]);
				histo_Ratio[NVar][NStep][nCh]->Divide(histo_nReweight_Data[NVar][NStep][nCh],histo_nReweight_MC[NVar][NStep][nCh],1,1,"b");

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
				histo_Ratio[NVar][NStep][nCh]->SetAxisRange(0.5,2,"y");
				histo_Ratio[NVar][NStep][nCh]->Draw("e");
				/*auto rp = new TRatioPlot(histo_MC[NVar][NStep][nCh],histo_RealData[NVar][NStep][nCh]);
				  rp->Draw();
				  canv_[NVar][NStep][nCh]->Update();*/
				canv_[NVar][NStep][nCh]->cd();
				canv_[NVar][NStep][nCh]->SaveAs(Save_dir+Variable[NVar]+"_"+Channel_txt[nCh]+"_"+Step_txt[NStep]+".png");
			}
		}
	}
	cout<<"13 TeV"<<endl;
}
