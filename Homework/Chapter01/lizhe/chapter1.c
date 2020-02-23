      1 #include <math.h>
      2 void tree(){
      3   const Double_t D=500.;//cm, distance between target and the scin.(Center)
      4   const Double_t L=100.;//cm, half length of the scin.
      5   const Double_t dD=5.;//cm, thickness of the scin.
      6   const Double_t TRes=1.;//ns, time resolution(FWHM) of the scintillator.
      7   const Double_t Lambda=380.;//cm, attenuation lenght of the scin.
      8   const Double_t QRes=0.1;//relative energy resolution(FWHM) of the scin.
      9   const Double_t Vsc=7.5;//ns/cm, speed of light in the scin.
     10   const Double_t En0=100;//MeV, average neutron energy
     11   const Double_t EnRes=50.;//MeV, energy spread of neutron(FWHM)
     12   const Double_t Eg0=1;//MeV, gamma energy
     13   const Double_t Rg=0.3;//ratio of gamma,ratio of neutron 1-Rg
     14   
     15 
     16   Double_t tofc;
     17   Double_t qx;//能量值算出的位置
     18   Double_t ce;//利用绝对标定的方法求出的能量  
     19  //1. 定义新ROOT文件，声明新的Tree
     20   TFile *opf=new TFile("tree.root","recreate");//新文件tree.root，指针 *opf
     21   TTree *opt=new TTree("tree","tree structure");//新tree，指针 *opt
     22 
     23  //2. 声明在tree结构中定义需要的变量分支
     24   Double_t x;//入射位置
     25   Double_t e;//能量
     26   int pid;    //粒子种类，n:pid=1,g:pid=0
     27   Double_t tof, ctof;//TOF:粒子实际飞行时间，cTOF：计算得到的TOF
     28   Double_t tu, td;
     29   Double_t qu, qd;
     30 
     31   Double_t tu_off=5.5;//time offset，//PMT的度越时间+电缆上的传输时间
     32   Double_t td_off=20.4;//time offset
     33 
     34   //3. 将变量地址添加到tree结构中
     35     //第一个参数为变量名称，第二个为上面定义的变量地址，第三个为变量的类型说明，D>        表示Double_t。
     36   opt->Branch("tofc", &tofc, "tofc/D");
     37   opt->Branch("ce", &ce, "ce/D");
     38   opt->Branch("x", &x, "x/D");
    39   opt->Branch("e", &e, "e/D");
     40   opt->Branch("tof", &tof, "tof/D");
     41   opt->Branch("ctof",&ctof,"ctof/D");
     42   opt->Branch("pid", &pid, "pid/I");
     43   opt->Branch("tu", &tu, "tu/D");
     44   opt->Branch("td", &td, "td/D");
     45   opt->Branch("qu", &qu, "qu/D");
     46   opt->Branch("qd", &qd, "qd/D");
     47   opt->Branch("qx", &qx, "qx/D");
     48 
     49 
     50 // histogram，ROOT文件中除了TTree结构外，还可存储histogram，graph等
     51   TH1D *hctof=new TH1D("hctof","neutron time of flight",1000,0,200);
     52   TRandom3 *gr=new TRandom3(0);//声明随机数
     53 
     54   //4. 循环，计算变量的值，逐事件往tree结构添加变量值。
     55   for(int i=0;i<100000;i++){
     56     x=gr->Uniform(-L, L);//均匀入射在中子探测器表面.
     57    // cout<<x<<endl;
     58     Double_t Dr=D+gr->Uniform(-0.5,0.5)*dD;//粒子在探测器厚度范围内均匀产生光信号
     59     Double_t d=TMath::Sqrt(Dr*Dr+x*x);//m, flight path
     60     if(gr->Uniform() < Rg) { //判断为gamma入射
     61        pid=0;
     62        e=Eg0;
     63        tof=3*(d*0.01);
     64     }
     65     else {  //neutron
     66         pid=1;
     67         e=gr->Gaus(En0, EnRes/2.35); // neutron
     68         tof=72./TMath::Sqrt(e)*(d*0.01);//ns
     69     }
     70     tu=tof+(L-x)/Vsc+gr->Gaus(0,TRes/2.35)+tu_off;
     71     td=tof+(L+x)/Vsc+gr->Gaus(0,TRes/2.35)+td_off;
     72     ctof=(tu+td)/2.;//simplified calculation.
     73    // cout<<ctof<<endl;
     74     hctof->Fill(ctof);
     75     Double_t q0=e*gr->Uniform();//energy of recoil proton in plas. 0-En
     76     qu=q0*TMath::Exp(-(L-x)/Lambda);
     77     qu=gr->Gaus(qu,qu*QRes/2.35);
     78     qd=q0*TMath::Exp(-(L+x)/Lambda);
     79     qd=gr->Gaus(qd,qd*QRes/2.35);
     80     qx=(Lambda/2)*log(qu/qd);
     81     opt->Fill();//5.将计算好的变量值填到Tree中
     82   }
     83   // 6.将数据写入root文件中
     84   hctof->Write();
     85   TCanvas *c1 = new TCanvas("c1","c1");
     86   hctof->Draw();
     87 
     88   ////tu-td
     89 
     90   TH1D *tdiff=new TH1D("tdiff","td-tu",140,-20,50);
     91   TCanvas *c2=new TCanvas("c2","c2");
     92 
     93   Long64_t nentries=opt->GetEntries();//得到事件总数
     94   for(Long64_t jentry=0; jentry<nentries; jentry++) {//对每个事件进行遍历
     95     opt->GetEntry(jentry);
     96     tdiff->Fill(td-tu);  // if(ng==1) tx->Fill(:tu-td), 只写入满足给定条件的事件          
     97   }
     98   tdiff->Draw();
     99   c2->Draw();
    100 
    101   ////寻找边界
    102 
    103   TCanvas *c3 = new TCanvas("c3","c3");
    104 
    105   TH1D *dtd=new TH1D("dtd","dt/dx",141,-20.25,50.25);
    106   Int_t i1=0;
    107 for(int i1=1;i1<tdiff->GetNbinsX();i1++) {
    108     Double_t df=tdiff->GetBinContent(i1+1)-tdiff->GetBinContent(i1);
    109     dtd->Fill(tdiff->GetBinLowEdge(i1+1),df);
    110 }
    111 dtd->Sumw2(0);
    112 dtd->Draw();
    113 dtd->Fit("gaus","","",-14,-9);//txl
    114 
    115   TF1 *f1 = new TF1("f1","[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)",39.5,43);
    116 f1->SetParameter(0,-350);
    117 f1->SetParameter(1,41.5);
    118 f1->SetParameter(2,0.5);
    119 dtd->Fit("f1","R+");
    120 dtd->Draw();
    121 c3->Draw();
    122 
    123 //计算中子能力
    124 TH2D *hgtofx=new TH2D("hgtofx","hgtofx",100,-120,120,100,39,48);
    125 TH1D *hgctof=new TH1D("hgctof","hgctof",100,39,48);
    126 
    127 Int_t jentry=0;
    128 
    129 for(Long64_t jentry=0; jentry<nentries; jentry++) {//对每个事件进行遍历
    130     opt->GetEntry(jentry);
    131     Double_t tx=3.749*(td-tu-14.932);
    132     if(ctof>40&& ctof<42.5) {
    133         hgtofx->Fill(tx,ctof);
    134         if(abs(tx)<5) hgctof->Fill(ctof);//gamma hits the center of the det.
    135     }
    136   }
    137 TCanvas *c4 = new TCanvas("c4","c4");
    138 hgtofx->Draw("cloz");
    139 c4->Draw();
    140 
    141 TCanvas *c5 = new TCanvas("c5","c5");
    142 hgctof->Draw();
    143 hgctof->Fit("gaus");
    144 c5->Draw();
    145 
    146 TH2D *hgtofcx=new TH2D("hgtofcx","corrected TOF",100,-120,120,100,0,50);
    147 TH1D *htofc=new TH1D("htofc","htof",200,0,100);
    148 
    149 Int_t jentry1=0;
    150 
    151 for(Long64_t jentry1=0; jentry1<nentries; jentry1++) {//对每个事件进行遍历
    152     opt->GetEntry(jentry1);
    153     Double_t tx=3.749*(td-tu-14.932);
    154     Double_t d=TMath::Sqrt(502.5*502.5+tx*tx);
    155     tofc=((ctof-24.5454)/d)*100.;//normalized to 100cm
    156     ce = 72*72/(tofc*tofc);//单位ns
    157     hgtofcx->Fill(tx,tofc);//gamma hits the center of the det.
    158     htofc->Fill(tofc);
    159   }
    160 hgtofcx->Draw("colz");//tofc与x之间无关联
    161 TCanvas *c6 = new TCanvas("c6","c6");
    162 c6->Draw();
    163 
    164 TCanvas *c7 = new TCanvas("c7","c7");
    165 c7->SetLogy();
    166 htofc->Draw();//修正后的飞行时间谱。
    167 c7->Draw();
    168   opt->Write();
    169   opf->Close();
    170 }
