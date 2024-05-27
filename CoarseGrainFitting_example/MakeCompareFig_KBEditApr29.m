clear
close all
clc
FS=28;
MS=800;

load Fit_HUES8WT_CpGsOnly_Chr1_3517.mat
Parameters
DistLength
[NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);
%DataName='HUES8_TKO_CpGsOnly_Chr1';
DataName='HUES8WT_CpGsOnly_Chr1';
%DataName='HUES8_DKO_P6_CpGsOnly_Chr1';
figtitle='HUES8WT';

%%
fn=['Compare_' figtitle];
figind=2;
[SSD] = JustPlotCompare(DataName,Parameters,DistLength,fn,figtitle,figind)
NetMeth=1;
DLm=DistLength(6);
DLu=DistLength(9);

figure(1)
subplot(1,3,1)
scatter(RelNetDeMeth,NetMeth,MS,'p','MarkerFaceColor','m','MarkerEdgeColor','k')
hold on
axis square
xlabel('Relative Net Demethylation')
ylabel('Normalized Net Methylation')
set(gca,'FontSize',FS)


subplot(1,3,2)
scatter(CoE,CoW,MS,'p','MarkerFaceColor','m','MarkerEdgeColor','k')
hold on
axis square
xlabel('Net Co. Demethylation')
ylabel('Net Co. Methylation')
set(gca,'FontSize',FS)

subplot(1,3,3)
scatter(DLu,DLm,MS,'p','MarkerFaceColor','m','MarkerEdgeColor','k')
hold on
axis square
set(gca,'FontSize',FS)
xlabel('Length. Co. Demethylation (bp)')
ylabel('Length. Co. Methylation (bp)')

%now load all the DKO fits, find the good ones for comparison
listing=dir('Fit_HUES8_DKO*.mat');
sz=size(listing);
NF=sz(1); %number of files to load
for ii=1:NF
    afn=listing(ii).name;
    load(afn,'Parameters','DistLength','SSD')
    KeepSSD(ii)=SSD;
    KeepPar(ii,:)=Parameters;
    KeepDL(ii,:)=DistLength;
end
[minSSD,keepi]=min(KeepSSD)
for ii=1:numel(keepi)
    afn=listing(ii).name;
    ind=keepi(ii);
    afn=listing(ind).name;
    load(afn,'Parameters','DistLength','SSD')

    %DataName='HUES8_TKO_CpGsOnly_Chr1';
    %DataName='HUES8WT_CpGsOnly_Chr1';
    DataName='HUES8_DKO_P6_CpGsOnly_Chr1';
    figtitle='HUES8 DNMT3A/B';
    fn=['Compare_DKO'];
    figind=3;
    [SSD] = JustPlotCompare(DataName,Parameters,DistLength,fn,figtitle,figind)

    [NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);
    DLm=DistLength(6);
    DLu=DistLength(9);
    NetMeth=1;

    figure(1)
    subplot(1,3,1)
    scatter(RelNetDeMeth,NetMeth,MS,'o','MarkerFaceColor','c','MarkerEdgeColor','k')
    hold on
    axis square

    subplot(1,3,2)
    scatter(CoE,CoW,MS,'o','MarkerFaceColor','c','MarkerEdgeColor','k')
    hold on

    subplot(1,3,3)
    scatter(DLu,DLm,MS,'o','MarkerFaceColor','c','MarkerEdgeColor','k')
    hold on
end

%now load all the TKO fits, find the good ones for comparison
KeepSSD=[];
KeepPar=[];
KeepDL=[];
listing=dir('Fit_HUES8_TKO*.mat');
sz=size(listing);
NF=sz(1); %number of files to load
for ii=1:NF
    afn=listing(ii).name;
    load(afn,'Parameters','DistLength','SSD')
    KeepSSD(ii)=SSD;
    KeepPar(ii,:)=Parameters;
    KeepDL(ii,:)=DistLength;
end

[minSSD,keepi]=min(KeepSSD)
for ii=1:numel(keepi)
    afn=listing(ii).name;
    ind=keepi(ii);
    afn=listing(ind).name;
    load(afn,'Parameters','DistLength','SSD')
    [NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);

    DataName='HUES8_TKO_CpGsOnly_Chr1';
    %DataName='HUES8WT_CpGsOnly_Chr1';
    %DataName='HUES8_DKO_P6_CpGsOnly_Chr1';
    figtitle='HUES8 TETx';
    fn=['Compare_TKO'];
    figind=4;
    [SSD] = JustPlotCompare(DataName,Parameters,DistLength,fn,figtitle,figind)

    DLm=DistLength(6);
    DLu=DistLength(9);
    NetMeth=1;

    figure(1)
    subplot(1,3,1)
    scatter(RelNetDeMeth,NetMeth,MS,'s','MarkerFaceColor','g','MarkerEdgeColor','k')
    hold on

    subplot(1,3,2)
    scatter(CoE,CoW,MS,'s','MarkerFaceColor','g','MarkerEdgeColor','k')
    hold on

    subplot(1,3,3)
    scatter(DLu,DLm,MS,'s','MarkerFaceColor','g','MarkerEdgeColor','k')
    hold on

end
%print -dpng CompareFitsFig

%%
% figure(2)
% set(gca,'FontName','Times')
% grid on
% yticks(0:0.25:1)
% figure(3)
% set(gca,'FontName','Times')
% grid on
% yticks(0:0.25:1)
% figure(4)
% set(gca,'FontName','Times')
% grid on
% yticks(0:0.25:1)
figure(1000)
for i =1:3
	subplot(1,3,i)
	set(gca,'FontName','Times')
	grid on
	ylim([0 1.05])
	yticks(0:0.25:1)
end

figure(1)
subplot(1,3,1)
set(gca,'FontName','Times')
grid on
subplot(1,3,2)
ylim([1.4 2.7])
yticks(1.4:0.2:2.6)
set(gca,'FontName','Times')
grid on
subplot(1,3,3)
ylim([35 52.5])
yticks(35:5:50)
l=legend('HUES8','DNMT3A/B','TETx');
l.FontSize=FS;
set(gca,'FontName','Times','FontSize',FS)
grid on
