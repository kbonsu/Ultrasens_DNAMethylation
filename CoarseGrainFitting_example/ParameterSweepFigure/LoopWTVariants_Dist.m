clear
close all
load('Fit_HUES8WT_CpGsOnly_Chr1_3517.mat','Parameters','DistLength');
p=Parameters;

Beta=1;%sqrt(WTCoErase/WTCoWrite);

Inner=[0.3:0.05:1.7];
%Inner=[0.5:0.1:1.4];
vals=Inner;
%vals=[0,Inner];
%vals=[0,0.5,0.9,0.95,1,1.05,1.1,1.5,2,4];

alphas=vals;%sqrt(vals); %this is how much to scale CoErase
betas=alphas*Beta; %how much to scale CoWrite

for jj=1:numel(betas)
    distm=betas(jj);
    for kk=1:numel(alphas)
        distu=alphas(kk)
        [n,K,ED50,flag,DLE,DLW] = WTVariants_Dist(distm,distu,p,DistLength);

        KeepnD(jj,kk)=n;
        KeepED50D(jj,kk)=ED50;
        KeepKD(jj,kk)=K;
        KeepCoE(jj,kk)=DLE;
        KeepCoW(jj,kk)=DLW;
        KeepFlag(jj,kk)=flag;
    end

end

WriteAxD=KeepCoW(:,1);
EraseAxD=KeepCoE(1,:);

nmin=min(KeepnD(KeepFlag==0));
nmax=max(KeepnD(KeepFlag==0));
KeepED50D(KeepFlag>0)=nan;
KeepnD(KeepFlag>0)=nan;

save KeepnD KeepnD
save KeepED50D KeepED50D
save KeepKD KeepKD
save WriteAxD WriteAxD
save EraseAxD EraseAxD

