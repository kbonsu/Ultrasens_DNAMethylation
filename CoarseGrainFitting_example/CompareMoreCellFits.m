%load the file that has the initial parameter set
close all
clear

Sc=1.02; %how much bigger than min SSD to allow
load Fit_HUES8WT_CpGsOnly_Chr1_7656.mat
Parameters
DistLength
[NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);

DLm=DistLength(6);
DLu=DistLength(9);

figure(1)
subplot(1,4,1)
scatter(RelNetDeMeth,NetMeth,'pk')
hold on

subplot(1,4,2)
scatter(CoE,CoW,'pk')
hold on

subplot(1,4,3)
scatter(DLu,DLm,'pk')
hold on

subplot(1,4,4)
scatter(NormCoE,NormCoW,'pk')
hold on

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
minSSD=min(KeepSSD)
keepi=find(KeepSSD<=(Sc*minSSD));
for ii=1:numel(keepi)
    afn=listing(ii).name;
    ind=keepi(ii);
    afn=listing(ind).name;
    load(afn,'Parameters','DistLength','SSD')
    [NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);
    DLm=DistLength(6);
    DLu=DistLength(9);

    figure(1)
    subplot(1,4,1)
    scatter(RelNetDeMeth,NetMeth,'m')
    hold on

    subplot(1,4,2)
    scatter(CoE,CoW,'m')
    hold on

    subplot(1,4,3)
    scatter(DLu,DLm,'m')
    hold on

    subplot(1,4,4)
    scatter(NormCoE,NormCoW,'m')
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
minSSD=min(KeepSSD)
keepi=find(KeepSSD<=(Sc*minSSD));
for ii=1:numel(keepi)
    afn=listing(ii).name;
    ind=keepi(ii);
    afn=listing(ind).name;
    load(afn,'Parameters','DistLength','SSD')
    [NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);
    DLm=DistLength(6);
    DLu=DistLength(9);

    figure(1)
    subplot(1,4,1)
    scatter(RelNetDeMeth,NetMeth,'g')
    hold on

    subplot(1,4,2)
    scatter(CoE,CoW,'g')
    hold on

    subplot(1,4,3)
    scatter(DLu,DLm,'g')
    hold on

    subplot(1,4,4)
    scatter(NormCoE,NormCoW,'g')
    hold on
end


%now load all the WT fits, find the good ones for comparison
KeepSSD=[];
KeepPar=[];
KeepDL=[];
listing=dir('Fit_HUES8WT*.mat');
sz=size(listing);
NF=sz(1); %number of files to load
for ii=1:NF
    afn=listing(ii).name;
    load(afn,'Parameters','DistLength','SSD')
    KeepSSD(ii)=SSD;
    KeepPar(ii,:)=Parameters;
    KeepDL(ii,:)=DistLength;
end
minSSD=min(KeepSSD)
keepi=find(KeepSSD<=(Sc*minSSD));
for ii=1:numel(keepi)
    afn=listing(ii).name;
    ind=keepi(ii);
    afn=listing(ind).name;
    load(afn,'Parameters','DistLength','SSD')
    [NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(Parameters,DistLength);
    DLm=DistLength(6);
    DLu=DistLength(9);

    figure(1)
    subplot(1,4,1)
    scatter(RelNetDeMeth,NetMeth,'k')
    hold on

    subplot(1,4,2)
    scatter(CoE,CoW,'k')
    hold on

    subplot(1,4,3)
    scatter(DLu,DLm,'k')
    hold on

    subplot(1,4,4)
    scatter(NormCoE,NormCoW,'k')
    hold on
end



