#include <math.h>
void tree(){  
  const Double_t D=500.;//cm, distance between target and the scin.(Center)
  const Double_t L=100.;//cm, half length of the scin.
  const Double_t dD=5.;//cm, thickness of the scin.
  const Double_t TRes=1.;//ns, time resolution(FWHM) of the scintillator.
  const Double_t Lambda=380.;//cm, attenuation lenght of the scin.
  const Double_t QRes=0.1;//relative energy resolution(FWHM) of the scin. 
  const Double_t Vsc=7.5;//ns/cm, speed of light in the scin.
  const Double_t En0=100;//MeV, average neutron energy
  const Double_t EnRes=50.;//MeV, energy spread of neutron(FWHM)
  const Double_t Eg0=1;//MeV, gamma energy  
  const Double_t Rg=0.3;//ratio of gamma,ratio of neutron 1-Rg 
  

  Double_t tofc;
  Double_t qx;//能量值算出的位置
  Double_t ce;//利用绝对标定的方法求出的能量
 //1. 定义新ROOT文件，声明新的Tree 
  TFile *opf=new TFile("tree.root","recreate");//新文件tree.root，指针 *opf
  TTree *opt=new TTree("tree","tree structure");//新tree，指针 *opt

 //2. 声明在tree结构中定义需要的变量分支
  Double_t x;//入射位置
  Double_t e;//能量
  int pid;    //粒子种类，n:pid=1,g:pid=0
  Double_t tof, ctof;//TOF:粒子实际飞行时间，cTOF：计算得到的TOF
  Double_t tu, td;
  Double_t qu, qd;

  Double_t tu_off=5.5;//time offset，//PMT的度越时间+电缆上的传输时间
  Double_t td_off=20.4;//time offset

  //3. 将变量地址添加到tree结构中
    //第一个参数为变量名称，第二个为上面定义的变量地址，第三个为变量的类型说明，D表示Double_t。
  opt->Branch("tofc", &tofc, "tofc/D");
  opt->Branch("ce", &ce, "ce/D");
  opt->Branch("x", &x, "x/D");
  opt->Branch("e", &e, "e/D");
  opt->Branch("tof", &tof, "tof/D");
  opt->Branch("ctof",&ctof,"ctof/D");
  opt->Branch("pid", &pid, "pid/I");
  opt->Branch("tu", &tu, "tu/D");
  opt->Branch("td", &td, "td/D");
  opt->Branch("qu", &qu, "qu/D"); 
  opt->Branch("qd", &qd, "qd/D");  
  opt->Branch("qx", &qx, "qx/D");


// histogram，ROOT文件中除了TTree结构外，还可存储histogram，graph等
  TH1D *hctof=new TH1D("hctof","neutron time of flight",1000,0,200);
  TRandom3 *gr=new TRandom3(0);//声明随机数

  //4. 循环，计算变量的值，逐事件往tree结构添加变量值。
  for(int i=0;i<100000;i++){
    x=gr->Uniform(-L, L);//均匀入射在中子探测器表面.
   // cout<<x<<endl;
    Double_t Dr=D+gr->Uniform(-0.5,0.5)*dD;//粒子在探测器厚度范围内均匀产生光信号
    Double_t d=TMath::Sqrt(Dr*Dr+x*x);//m, flight path
    if(gr->Uniform() < Rg) { //判断为gamma入射
       pid=0;
       e=Eg0;
       tof=3*(d*0.01);
    }
    else {  //neutron
        pid=1;
        e=gr->Gaus(En0, EnRes/2.35); // neutron
        tof=72./TMath::Sqrt(e)*(d*0.01);//ns
    }
    tu=tof+(L-x)/Vsc+gr->Gaus(0,TRes/2.35)+tu_off;
    td=tof+(L+x)/Vsc+gr->Gaus(0,TRes/2.35)+td_off;
    ctof=(tu+td)/2.;//simplified calculation.
   // cout<<ctof<<endl;
    hctof->Fill(ctof);
    Double_t q0=e*gr->Uniform();//energy of recoil proton in plas. 0-En
    qu=q0*TMath::Exp(-(L-x)/Lambda);
    qu=gr->Gaus(qu,qu*QRes/2.35);
    qd=q0*TMath::Exp(-(L+x)/Lambda);
    qd=gr->Gaus(qd,qd*QRes/2.35);
    qx=(Lambda/2)*log(qu/qd);
    opt->Fill();//5.将计算好的变量值填到Tree中
  }
  // 6.将数据写入root文件中
  hctof->Write();
  TCanvas *c1 = new TCanvas("c1","c1");
  hctof->Draw();
  
  ////tu-td
  
  TH1D *tdiff=new TH1D("tdiff","td-tu",140,-20,50);  
  TCanvas *c2=new TCanvas("c2","c2");

  Long64_t nentries=opt->GetEntries();//得到事件总数
  for(Long64_t jentry=0; jentry<nentries; jentry++) {//对每个事件进行遍历
    opt->GetEntry(jentry);
    tdiff->Fill(td-tu);  // if(ng==1) tx->Fill(:tu-td), 只写入满足给定条件的事件      
  }
  tdiff->Draw();
  c2->Draw();

  ////寻找边界

  TCanvas *c3 = new TCanvas("c3","c3");

  TH1D *dtd=new TH1D("dtd","dt/dx",141,-20.25,50.25);
  Int_t i1=0;
for(int i1=1;i1<tdiff->GetNbinsX();i1++) {
    Double_t df=tdiff->GetBinContent(i1+1)-tdiff->GetBinContent(i1);
    dtd->Fill(tdiff->GetBinLowEdge(i1+1),df);
}
dtd->Sumw2(0);
dtd->Draw();
dtd->Fit("gaus","","",-14,-9);//txl

  TF1 *f1 = new TF1("f1","[0]*TMath::Exp(-0.5*((x-[1])/[2])^2)",39.5,43);
f1->SetParameter(0,-350);
f1->SetParameter(1,41.5);
f1->SetParameter(2,0.5);
dtd->Fit("f1","R+");
dtd->Draw();
c3->Draw();

//计算中子能力
TH2D *hgtofx=new TH2D("hgtofx","hgtofx",100,-120,120,100,39,48);
TH1D *hgctof=new TH1D("hgctof","hgctof",100,39,48);

Int_t jentry=0;

for(Long64_t jentry=0; jentry<nentries; jentry++) {//对每个事件进行遍历
    opt->GetEntry(jentry);
    Double_t tx=3.749*(td-tu-14.932);
    if(ctof>40&& ctof<42.5) { 
        hgtofx->Fill(tx,ctof);
        if(abs(tx)<5) hgctof->Fill(ctof);//gamma hits the center of the det.
    }
  }
TCanvas *c4 = new TCanvas("c4","c4");
hgtofx->Draw("cloz");
c4->Draw();

TCanvas *c5 = new TCanvas("c5","c5");
hgctof->Draw();
hgctof->Fit("gaus");
c5->Draw();

TH2D *hgtofcx=new TH2D("hgtofcx","corrected TOF",100,-120,120,100,0,50);
TH1D *htofc=new TH1D("htofc","htof",200,0,100);

Int_t jentry1=0;

for(Long64_t jentry1=0; jentry1<nentries; jentry1++) {//对每个事件进行遍历
    opt->GetEntry(jentry1);
    Double_t tx=3.749*(td-tu-14.932);
    Double_t d=TMath::Sqrt(502.5*502.5+tx*tx);
    tofc=((ctof-24.5454)/d)*100.;//normalized to 100cm
    ce = 72*72/(tofc*tofc);//单位ns
    hgtofcx->Fill(tx,tofc);//gamma hits the center of the det.
    htofc->Fill(tofc);
  }
hgtofcx->Draw("colz");//tofc与x之间无关联
TCanvas *c6 = new TCanvas("c6","c6");
c6->Draw();

TCanvas *c7 = new TCanvas("c7","c7");
c7->SetLogy();
htofc->Draw();//修正后的飞行时间谱。
c7->Draw();
  opt->Write();
  opf->Close();
}
