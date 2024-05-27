close all
clear

%load in the WT parameters that were used as baseline
load('Fit_HUES8WT_CpGsOnly_Chr1_3517.mat','Parameters','DistLength');
p=Parameters;
k1 = p(1);
k2=p(2);
k3=p(3);
k4=p(4);
k6=p(6);
k8=p(8);
k9=p(9);
k11=p(11);
k13=p(13);
k14=p(14);
DLm=DistLength(6);
DLu=DistLength(9);

FS=24;

Nt=24; %cell cycle time
HyperVal=0.8; %cutoff for hypermeth
HypoVal=0.2; %cutoff for hypometh
NCpG=27;

k2t=k2+k14;
k4t=k4+k13;
WTW=1;%k1*k3;
WTE=(k2t*k4t)/(k1*k3);
WTCoE=k9*k11/(k2t*k4t);
WTCoW=k6*k8/(k1*k3);


load Keepn
load KeepED50
load KeepCoW
load KeepCoE
yax=KeepCoW(:,1);
xax=KeepCoE(1,:);


nmin=3; %min Hill coeff for colorscale
nmax=5.5;
emin=0.01; %min ED50 for colorsacle
emax=0.5; 
KeepFlag=zeros(size(KeepED50));
KeepFlag(KeepED50==0)=1;
nmin=min(Keepn(KeepFlag==0));
nmax=max(Keepn(KeepFlag==0));
Keepn(KeepFlag>0)=nan;
KeepED50(KeepFlag>0)=nan;

figure(4)
%subplot(2,3,1)
subplot(1,3,1)
pcolor(xax,yax,Keepn)
colorbar
shading flat
colorbar
title('Hill Coefficient')
xlabel('Net Co. Demethylation')
ylabel('Net Co. Methylation')
colormap(parula)
caxis([nmin nmax])
axis square
hold on
a3=plot(WTCoE,WTCoW,'p','MarkerSize',30,'MarkerFaceColor','m','MarkerEdgeColor','k')
set(gca,'FontSize',FS,'FontName','Times')


%subplot(2,3,2)
figure(5)
subplot(1,3,1)
pcolor(xax,yax,KeepED50)
colorbar
shading flat
colorbar
title('ED50')
xlabel('Net Co. Demethylation')
ylabel('Net Co. Methylation')
hold on
a3=plot(WTCoE,WTCoW,'p','MarkerSize',30,'MarkerFaceColor','m','MarkerEdgeColor','k')
colormap(parula)
caxis([emin emax])
axis square
set(gca,'FontSize',FS,'FontName','Times')
%set(gca,'Position',getpos)

%print -dpng ParamSweepFig_CoMeth

%MaKe FIg 2- the net methylation
load KeepnN
load KeepED50N
load KeepKN 
load KeepNormW 
load KeepRelE 

yax=KeepNormW(:,1);
xax=KeepRelE(9,:);

KeepFlag=zeros(size(KeepED50N));
KeepFlag(KeepED50N==0)=1;
nmin=min(KeepnN(KeepFlag==0));
nmax=max(KeepnN(KeepFlag==0));
KeepED50N(KeepFlag>0)=nan;
KeepnN(KeepFlag>0)=nan;

figure(4)
%subplot(2,3,3)
subplot(1,3,2)
pcolor(xax,yax,KeepnN)
colorbar
shading flat
colorbar
title('Hill Coefficient')
xlabel('Relative Net Demethylation')
ylabel('Normalized Net Methylation')
colormap(parula)
caxis([nmin nmax])
axis square
hold on

WTW=1;%yax(9);
a3=plot(WTE,WTW,'p','MarkerSize',30,'MarkerFaceColor','m','MarkerEdgeColor','k')
set(gca,'FontSize',FS,'FontName','Times')
%getpos=get(gca,'Position');

%subplot(2,3,4)
figure(5)
subplot(1,3,2)
pcolor(xax,yax,KeepED50N)
shading flat
colorbar
title('ED50')
xlabel('Relative Net Demethylation')
ylabel('Normalized Net Methylation')
hold on
%a1=plot(CoEax,KeepEqW,'-y','LineWidth',2)
%a2=plot(CoEax,KeepEqW2,'-m','LineWidth',2)
a3=plot(WTE,WTW,'p','MarkerSize',30,'MarkerFaceColor','m','MarkerEdgeColor','k')
colormap(parula)
caxis([emin emax])
axis square
set(gca,'FontSize',FS,'FontName','Times')
%set(gca,'Position',getpos)
%print -dpng ParamSweepFig_NetMeth

%make FIg 3 - distance lengthscale

load KeepnD 
load KeepED50D 
load KeepKD 
load WriteAxD 
load EraseAxD 
nmax=10;

KeepFlag=zeros(size(KeepED50D));
KeepFlag(KeepED50D==0)=1;
nmin=min(KeepnD(KeepFlag==0));
nmax=max(KeepnD(KeepFlag==0));
KeepED50D(KeepFlag>0)=nan;
KeepnD(KeepFlag>0)=nan;

figure(4)
%subplot(2,3,5)
subplot(1,3,3)
pcolor(EraseAxD,WriteAxD,KeepnD)
colorbar
shading flat
colorbar
title('Hill Coefficient')
xlabel('Length Co. Demethylation (bp)')
ylabel('Length Co. Methylation (bp)')
colormap(parula)
caxis([nmin nmax])
axis square
hold on
WTWD=DLm;%WriteAxD(getwta);
WTED=DLu;%EraseAxD(getwtb);
a3=plot(WTED,WTWD,'p','MarkerSize',30,'MarkerFaceColor','m','MarkerEdgeColor','k')
set(gca,'FontSize',FS,'FontName','Times')
%getpos=get(gca,'Position');


%subplot(2,3,6)
figure(5)
subplot(1,3,3)
pcolor(EraseAxD,WriteAxD,KeepED50D)
colorbar
shading flat
colorbar
title('ED50')
xlabel('Length Co. Demethylation (bp)')
ylabel('Length Co. Methylation (bp)')
hold on
a3=plot(WTED,WTWD,'p','MarkerSize',30,'MarkerFaceColor','m','MarkerEdgeColor','k')
colormap(parula)
caxis([emin emax])
axis square
set(gca,'FontSize',FS,'FontName','Times')
%set(gca,'Position',getpos)
print -dpng ParamSweepFig_Dist


%uncomment below to plot some representative parameter sets

% Parameters=[k1,k2,k3,k4,0,k6,0,k8,k9,0,k11,0,k13,k14];
% DLu=DL;
% DLm=DL;
% DistLength=[0,0,0,0,DLm,DLm,DLm,DLm,DLu,DLu,DLu,DLu];
% JustPlotCurves_NoData(Parameters,DistLength,4)
% alphau=1.55;
% alpham=1;
% k9=k9*alphau;
% k11=k11*alphau;
% k6=k6*alpham;
% k8=k8*alpham;
% Parameters=[k1,k2,k3,k4,0,k6,0,k8,k9,0,k11,0,k13,k14];
% DLm=100;
% DLu=DL;
% DistLength=[0,0,0,0,DLm,DLm,DLm,DLm,DLu,DLu,DLu,DLu];
% 
% JustPlotCurves_NoData(Parameters,DistLength,4)
% % NCoE=k9*k11/(k2t*k4t);
% % NCoW=k6*k8/(k1*k3);
% % figure(1)
% % subplot(2,1,1)
% % a4=plot(NCoE,NCoW,'p','MarkerSize',30,'MarkerFaceColor','c','MarkerEdgeColor','k')
% 
% % Parameters=[k1,k2,k3,k4,0,k6,0,k8,k9,0,k11,0,k13,k14];
% % k9=k9*3;
% % Parameters=[k1,k2,k3,k4,0,k6,0,k8,k9,0,k11,0,k13,k14];
% % DLm=170;
% % %DLu=50;
% % DistLength=[0,0,0,0,DLm,DLm,DLm,DLm,DLu,DLu,DLu,DLu];
% % 
% % JustPlotCurves_NoData(Parameters,DistLength,4)