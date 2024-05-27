clear
close all


load('Fit_HUES8WT_CpGsOnly_Chr1_3517.mat','Parameters','DistLength');
Beta=1;

%Inner=[0.3:0.05:2];
Inner=[0.6:.05:1.4];
%vals=[0,0.5,0.9,0.95,1,1.05,1.1,1.5,2,4];
%vals=[0,Inner];
vals=Inner;
alphas=sqrt(vals); %this is how much to scale Net Meth
betas=[alphas*Beta]; %how much to scale NetDeMeth

for jj=1:numel(alphas)
    alpha=alphas(jj)
    for kk=1:numel(betas)
        beta=betas(kk)
        [n,K,ED50,flag,NormW,RelE] = WTVariants_NetMeth(alpha,beta,Parameters,DistLength);

        KeepnN(jj,kk)=n;
        KeepED50N(jj,kk)=ED50;
        KeepKN(jj,kk)=K;
        KeepFlag(jj,kk)=flag;
        KeepNormW(jj,kk)=NormW;
        KeepRelE(jj,kk)=RelE;
    end
end
yax=KeepNormW(:,1);
xax=KeepRelE(1,:);

nmin=min(KeepnN(KeepFlag==0));
nmax=max(KeepnN(KeepFlag==0));
KeepED50N(KeepFlag>0)=nan;
KeepnN(KeepFlag>0)=nan;

save KeepnN KeepnN
save KeepED50N KeepED50N
save KeepKN KeepKN
save KeepNormW KeepNormW
save KeepRelE KeepRelE
