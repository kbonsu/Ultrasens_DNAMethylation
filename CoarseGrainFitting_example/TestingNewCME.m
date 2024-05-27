function [PVecMSM,MBins,PVec,Prob_ind,IndCpGp] = TestingNewCME(NCpG,Parameters,CpGPositions,DistLength)
% Nt=24; %cell cycle time
% HyperVal=0.8; %cutoff for hypermeth
% HypoVal=0.2; %cutoff for hypometh
% NCpG=27;
% 
% %Calculate the standard model parameters
% %assume we know k2, k3
% k3 = 2;
% k2=0.1;
% k13=1/Nt;
% k14=1/(2*Nt);
% %make sure that k3*ph>k13*pm! Or things break--> negative value for k4
% 
% k1=0.151;
% k4=0.2083;
% DL=32;
% %distance parameters
% DistLength=[0,0,0,0,DL,DL,DL,DL,DL,DL,DL,DL];
% 
% CoWrite=3.7;%4.8;%6.3921;%5.3667;%1;%5.4;%0.06;%0.0381;
% CoErase=4.5;%.4;%5.8;%8.2785;%6.5432;%1.65;%6.5;%0.35;%0.2207;
% Parameters=[k1,k2,k3,k4,0,CoWrite,0,CoWrite,CoErase,0,CoErase,0,k13,k14];

k=Parameters;
Dists=pdist(CpGPositions');
DistMat=squareform(Dists);
DistMat=DistMat+eye(NCpG)*1E6; %add a big number to the diagonal to silence self interaction


%enumerate the types--this determines the size of the coarse-grained master
%eqn
ct=0;
for nm=0:NCpG
    for nh = 0:NCpG-nm
        nu=NCpG-nm-nh;
        ct=ct+1;
        Types(ct,:)=[nu,nh,nm];
    end
end
NMacro=ct; %number of macrostates (= # unique types)

%this loop pre-calculates the field exerted by the various sites. Currently
%it assumes that the distance function decays exponentially with decay
%constant IL. IL is obtained from the inverse of the input parameters DistLength
% However, IL can vary for the different collaborative reactions.
%In the limit that DistLength->\infty, this is the same as assuming that all sites
%contribute equally (i.e., no distance dependence).

%find out how many unique distance values there are, for which the fields
%need to be computed
collabinds=[5:12];
[Vals,ia,ic]=unique(DistLength(collabinds));
DistFuncAssign=zeros(1,12);
DistFuncAssign(collabinds)=ic;
%DistFuncAssign contains the distance function assignments for each
%reaction. e.g., the distance for the 8th reaction is
%Vals(DistFuncAssign(8))

IL=1./Vals; %Contains the inverse lengths. decay constant = inverse length
NumFuncs=numel(IL); %the number of different distance functions specified
%MCpG=ceil(NCpG/2);
fu_ar=zeros(NMacro,NumFuncs); %field is calculated depending on macrostate,
%and which reaction number
fh_ar=fu_ar; fm_ar = fu_ar;
for M=1:NMacro
    nu=Types(M,1);
    nh=Types(M,2);
    nm=Types(M,3);
    % Multi(M)=factorial(NCpG)/(factorial(nu)*factorial(nh)*factorial(nm)); %the multiplicity of the type
    Probs(M,:)=Types(M,:)/NCpG; %the probability of pu,ph,pm within this type
    Meth(M)=nh*0.5+nm; %the methylation level of this type
    pu=Probs(M,1);
    ph=Probs(M,2);
    pm=Probs(M,3);
    for DL=1:NumFuncs %this loops over the collaborative reactions distance functions
        ExpMat=exp(-IL(DL)*DistMat);
        fu_ar(M,DL) = pu*sum(ExpMat(:))/NCpG;
        fh_ar(M,DL) = ph*sum(ExpMat(:))/NCpG;
        fm_ar(M,DL) = pm*sum(ExpMat(:))/NCpG;
    end
end

%Standard Model for 1 CpG (independent of the rest of the island
StdMod=[-k(1), k(2) + k(14), 0;
    k(1), -(k(2)+k(14)+k(3)), k(4) + k(13);
    0, k(3), -(k(4) + k(13))];

%solve the steady state probability for the single INDEPENDENT site
% (without the collaborative reactions. This would be the expected [pu,ph,pm] for
% a site that is very far away from other CpGs (not interacting)
[VV,DD]=eigs(StdMod,3,1E-8);
Prob_ind=VV(:,1)/sum(VV(:,1)); %the probabilities pu,ph,pm for the independent site

MSM=zeros(NMacro);
RxnStoich=[-1,1,0; %u ->h  %K1=RM1(2,1); %the u -> h reaction
    1,-1,0; %h->u   %K2=RM1(1,2); %h ->u
    0,-1,1; %h-> m    %K3=RM1(3,2); %h->m
    0,1,-1]; %m->h  %K4=RM1(2,3); %m->h

for M=1:NMacro
    nu=Types(M,1);
    nh=Types(M,2);
    nm=Types(M,3);

    fu9=fu_ar(M,DistFuncAssign(9)); %set the field --contingent on M and the collab rxn
    fu11=fu_ar(M,DistFuncAssign(11));
    fh5=fh_ar(M,DistFuncAssign(5));
    fh7=fh_ar(M,DistFuncAssign(7));
    fh10=fh_ar(M,DistFuncAssign(10));
    fh12=fh_ar(M,DistFuncAssign(12));
    fm6=fm_ar(M,DistFuncAssign(6));
    fm8=fm_ar(M,DistFuncAssign(8));

    Collab=[-k(7)*fh7-k(8)*fm8, k(9)*fu9+k(10)*fh10, 0;
        k(7)*fh7 + k(8)*fm8, -k(9)*fu9 - k(10)*fh10 - k(5)*fh5 - k(6)*fm6, k(12)*fh12 + k(11)*fu11;
        0, k(5)*fh5 + k(6)*fm6, -k(12)*fh12 - k(11)*fu11];

    RM1=StdMod+Collab; %small rate matrix for single site

    %run through the 4 types of reactions
    %%Rxn 1
    rx = 1;
    NewType=Types(M,:)+RxnStoich(rx,:);
    if prod(NewType>=0) %makes sure it's an allowed reaction
        % rate1=0; %initialize the rate for this type of reaction
        [val,ind]=ismember(NewType,Types,'rows'); %figure out which type it lands on

        rate1=RM1(2,1); %the rate of this reaction if it occurs at this site
        MSM(ind,M)=MSM(ind,M)+nu*rate1; %insert the computed rate in the rate matrix
    end

    %%Rxn 2
    rx = 2;
    NewType=Types(M,:)+RxnStoich(rx,:);
    if prod(NewType>=0) %makes sure it's an allowed reaction
        [val,ind]=ismember(NewType,Types,'rows'); %figure out which type it lands on

        rate1=RM1(1,2); %the rate of this reaction if it occurs at this site
        MSM(ind,M)=MSM(ind,M)+nh*rate1; %insert the computed rate in the rate matrix
    end

    %%Rxn 3
    rx = 3;
    NewType=Types(M,:)+RxnStoich(rx,:);
    if prod(NewType>=0) %makes sure it's an allowed reaction
        [val,ind]=ismember(NewType,Types,'rows'); %figure out which type it lands on

        rate1=RM1(3,2); %the rate of this reaction if it occurs at this site
        MSM(ind,M)=MSM(ind,M)+nh*rate1; %insert the computed rate in the rate matrix
    end

    %%Rxn 4
    rx = 4;
    NewType=Types(M,:)+RxnStoich(rx,:);
    if prod(NewType>=0) %makes sure it's an allowed reaction
        [val,ind]=ismember(NewType,Types,'rows'); %figure out which type it lands on

        rate1=RM1(2,3); %the rate of this reaction if it occurs at this site
        MSM(ind,M)=MSM(ind,M)+nm*rate1; %insert the computed rate in the rate matrix
    end
end

% %%Making columns in transition matrix sum to 0
MSM=MSM-diag(sum(MSM));
%get the steady state probability
[VV,DD]=eigs(MSM,3,1E-8);
PVec=VV(:,1)/sum(VV(:,1));

% %Now we make an even more coarse-grained pr27bability vector
%set up more coarse grained bins
MBins=[0:0.5:NCpG]; %the macrostate bins

NMicro=numel(PVec); %the old "macrostates" became microstates
NMacro=numel(MBins);
Chi=zeros(NMicro,NMacro);

for bin = 1:numel(MBins)
    getinds=find(Meth==MBins(bin));
    Chi(getinds,bin)=1; %membership matrix
end

%Now calculate the MSM (the coarse grained rate matrix)--this is what we
%are trying to approximate without enumerating all the states!
Chi_Wtd=Chi.*repmat(PVec,1,NMacro);
PVecMSM=sum(Chi_Wtd); %coarse grained prob sums the contributions from each microstate

%plot(MBins/NCpG,PVecMSM)
%hold on
%
% SumChiWtd=sum(Chi_Wtd,1);
% NormWtChi=Chi_Wtd./repmat(SumChiWtd,NMicro,1);
% %Calculate the MSM directly from the rate matrix
% MSMBin=Chi'*MSM*NormWtChi; %this coarser version of MSM in terms of Meth bins --for validation
%We want to recover the pu,ph,pm for individual CpGs from PVecMSM (which is

%the probability of total methylation states
putot=sum(Types(:,1).*PVec)/NCpG;
phtot=sum(Types(:,2).*PVec)/NCpG;
pmtot=sum(Types(:,3).*PVec)/NCpG;
IndCpGp=[putot,phtot,pmtot];



