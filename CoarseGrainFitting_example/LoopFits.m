clear
close all
NLoops=10;
for ii=1:NLoops
    ii

%load the file that has the initial parameter set
load Fit_HUES8WT_CpGsOnly_Chr1_7656.mat
%add some randomness
Parameters
DistLength

DataName='HUES8_TKO_CpGsOnly_Chr1';
%DataName='HUES8WT_CpGsOnly_Chr1';
%DataName='HUES8_DKO_P6_CpGsOnly_Chr1';

    %[p,SSD] = Fitting_Tues(Parameters,DistLength,DataName)
     [p,SSD] = Fitting_Thurs(Parameters,DistLength,DataName)

end

%save KeepSSD KeepSSD
%save KeepPar KeepPar
