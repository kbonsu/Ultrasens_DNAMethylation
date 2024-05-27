clear
close all
clc
FS=28;
MS=800;

%% Modified to Only output Collab. model
load Fit_HUES8WT_CpGsOnly_Chr1_3517.mat
Parameters;
DistLength;
[NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);
%DataName='HUES8_TKO_CpGsOnly_Chr1';
DataName='HUES8WT_CpGsOnly_Chr1';
%DataName='HUES8_DKO_P6_CpGsOnly_Chr1';
figtitle='HUES8WT Collaborative';

%%
fn=['Compare_' figtitle];
figind=3;
[SSD] = JustPlotCompare(DataName,Parameters,DistLength,fn,figtitle,figind)


%% Modified to Only output standard model
load Fit_HUES8WT_CpGsOnly_Chr1_3517.mat
Parameters;
Parameters(5:12)=0;
DistLength;
[NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);
%DataName='HUES8_TKO_CpGsOnly_Chr1';
DataName='HUES8WT_CpGsOnly_Chr1';
%DataName='HUES8_DKO_P6_CpGsOnly_Chr1';
figtitle='HUES8WT Standard';

%%
fn=['Compare_' figtitle];
figind=2;
[SSD] = JustPlotCompare(DataName,Parameters,DistLength,fn,figtitle,figind)
%%
figure(1000)
for i=1:2
	subplot(1,2,i)
	set(gca,'FontName','Times','FontSize',FS+4)
	grid on
	ylim([0 1.05])
	yticks(0:0.25:1)
	xticks(0:0.1:0.5)
end