clear
close all

load('Fit_HUES8WT_CpGsOnly_Chr1_3517.mat','Parameters','DistLength');

Beta=1;%sqrt(WTCoErase/WTCoWrite);

Inner=[0.3:0.05:2];
vals=Inner;
alphas=sqrt(vals); %this is how much to scale CoErase
betas=alphas*Beta; %how much to scale CoWrite


for jj=1:numel(alphas)
    alpham=betas(jj)
    for kk=1:numel(alphas)
    alphau=alphas(kk)

[n,K,ED50,flag,CoE,CoW,CoE2,CoW2] = WTVariants3(alpham,alphau,Parameters,DistLength);
Keepn(jj,kk)=n;
KeepED50(jj,kk)=ED50;
KeepK(jj,kk)=K;
KeepCoE(jj,kk)=CoE;
KeepCoW(jj,kk)=CoW;
KeepCoE2(jj,kk)=CoE2;
KeepCoW2(jj,kk)=CoW2;
KeepFlag(jj,kk)=flag;
    end

end
yax=KeepCoW(:,1);
xax=KeepCoE(1,:);
yax2=KeepCoW2(:,1);
xax2=KeepCoE2(1,:);

nmin=min(Keepn(KeepFlag==0));
nmax=max(Keepn(KeepFlag==0));
Keepn(KeepFlag>0)=nan;
KeepED50(KeepFlag>0)=nan;


save Keepn Keepn
save KeepED50 KeepED50
save KeepK KeepK
save KeepCoW KeepCoW
save KeepCoE KeepCoE






