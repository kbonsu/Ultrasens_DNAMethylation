function [SSD] = JustPlotCompare(DataName,Parameters,DistLength,fn,FTitle,fi)
LW=2;
MS=10;
FS=18;
HyperVal=0.8; %cutoff for hypermeth
HypoVal=0.2; %cutoff for hypometh
NCpG=27; %number of CpGs in the model

DataFile=['Save_' DataName '.mat'];
load(DataFile,'DataStruct')

[ModelStruct] = Run_Model_MorePoints(Parameters,DistLength);
SSD = Get_Error(ModelStruct,DataStruct);

dens=ModelStruct.densityvals;
Hyper=ModelStruct.Hyper;
Hypo=ModelStruct.Hypo;
Inter=ones(size(Hyper))-Hyper-Hypo;

lsMark={':ro',':bo',':ko';'-r','-b','-k'};


%figure(fi)

% %Edit to Create subplot figures
figure(1000)
subplot(1,3,(fi-1))

% plot(dens,Hyper,':o','MarkerSize',MS,'MarkerFaceColor','r')
% hold on
% plot(dens,Hypo,':o','MarkerSize',MS,'MarkerFaceColor','b')
% a1=plot(dens,Inter,':o','MarkerSize',MS,'MarkerFaceColor','k')
% axis square

plot(dens,Hyper,':ro','MarkerSize',MS,'MarkerFaceColor','r')
hold on
plot(dens,Hypo,':bo','MarkerSize',MS,'MarkerFaceColor','b')
a1=plot(dens,Inter,':ko','MarkerSize',MS,'MarkerFaceColor','k')
axis square

% kin=1;
% plot(dens,Hyper,lsMark{kin,1},'LineWidth',MS,'MarkerFaceColor','r')
% hold on
% plot(dens,Hypo,lsMark{kin,2},'LineWidth',MS,'MarkerFaceColor','b')
% a1=plot(dens,Inter,lsMark{kin,3},'LineWidth',MS,'MarkerFaceColor','k')
% axis square

dens=DataStruct.densityvals;
Hyper=DataStruct.Hyper;
Hypo=DataStruct.Hypo;
Inter=ones(numel(Hyper),1)-Hyper-Hypo;

plot(dens,Hyper,'-r','LineWidth',LW)
hold on
plot(dens,Hypo,'-b','LineWidth',LW)
a2=plot(dens,Inter,'-k','LineWidth',LW)

% kinn=2;
% plot(dens,Hyper,lsMark{kinn,1},'LineWidth',LW)
% hold on
% plot(dens,Hypo,lsMark{kinn,2},'LineWidth',LW)
% a2=plot(dens,Inter,lsMark{kinn,3},'LineWidth',LW)

title(FTitle)
txt='Hyper (CpGMe>0.8)';
posx=0.225;
posy=0.60;
text(posx,posy,txt,'FontSize',FS,'Color','r')
txt='Hypo (CpGMe<0.2)';
posy=0.70;
text(posx,posy,txt,'FontSize',FS,'Color','b')

text(posx,0.5,['SSD =' num2str(SSD)],'FontSize',FS,'FontName','Times')
xlim([0 0.525])
ylim([0 1])
ylabel('Fraction of CpGs')
xlabel('Normalized CpG Density')
lpos=get(legend,'Position');

%if fi==4
    legend([a1 a2], 'Model','Data','Position',[lpos(1) 0.3 lpos(3) lpos(4)])
%else b=gca;
%    legend(b,'off')
%end
set(gca,'FontSize',FS)

print(fn,'-dpng')


    function SSD = Fit_Model(p,DataStruct,Parameters)
        Parameters(4)=p(1);
        Parameters(6)=p(2);
        Parameters(8)=p(3);
        Parameters(9)=p(4);
        Parameters(11)=p(5);

        [ModelStruct] = Run_Model(Parameters,DistLength);
        SSD = Get_Error(ModelStruct,DataStruct);
    end


    function SSD = Get_Error(ModelStruct,DataStruct)
        ModelMat=[ModelStruct.Hyper(:),ModelStruct.Hypo(:)];
        %we need to interpolate the experimental data to make it match the
        %model output
        xq=ModelStruct.densityvals;
        Hypervq=interp1(DataStruct.densityvals,DataStruct.Hyper,xq);
        Hypovq=interp1(DataStruct.densityvals,DataStruct.Hypo,xq);
        DataMat=[Hypervq(:),Hypovq(:)];

        Diffs=ModelMat-DataMat;
        SSD=sum(Diffs(:).^2);
    end

    function [ModelStruct] = Run_Model(Parameters,DistLength)
        ds=[2,4,5,6,7,8,9,11,13,17,26,110];
        %ds=[2,4,6,8,11,17,110];
        ds=fliplr(ds);

        for loopd=1:numel(ds)
            d=ds(loopd);
            CpGPositions=[1:d:NCpG*d];
            Densities=CpGDensities_Function(CpGPositions,50);
            MeanDens(loopd)=mean(Densities);

            %call the function that computes the P(NetMeth) for the coarse-grained
            %approximate model
            [PVecMSM,MBins,PVec,Prob_ind,IndCpGp] = TestingNewCME(NCpG,Parameters,CpGPositions,DistLength);

            %make sure PVec is properly positive and normalized
            PVecMSM(PVecMSM<0)=0;
            PVecMSM=PVecMSM/sum(PVecMSM);
            MethRatio=MBins/NCpG;
            HyperInds=find(MethRatio>HyperVal);
            HypoInds=find(MethRatio<HypoVal);
            Hyper(loopd)=sum(PVecMSM(HyperInds));
            Hypo(loopd)=sum(PVecMSM(HypoInds));
        end

        ModelStruct.densityvals=MeanDens;
        ModelStruct.Hyper=Hyper;
        ModelStruct.Hypo=Hypo;
    end

    function [ModelStruct] = Run_Model_MorePoints(Parameters,DistLength)
        ds=[2,3,4,5,6,7,8,9,11,13,17,26,110];
        ds=fliplr(ds);

        for loopd=1:numel(ds)
            d=ds(loopd);
            CpGPositions=[1:d:NCpG*d];
            Densities=CpGDensities_Function(CpGPositions,50);
            MeanDens(loopd)=mean(Densities);

            %call the function that computes the P(NetMeth) for the coarse-grained
            %approximate model
            [PVecMSM,MBins,PVec,Prob_ind,IndCpGp] = TestingNewCME(NCpG,Parameters,CpGPositions,DistLength);

            %make sure PVec is properly positive and normalized
            PVecMSM(PVecMSM<0)=0;
            PVecMSM=PVecMSM/sum(PVecMSM);
            MethRatio=MBins/NCpG;
            HyperInds=find(MethRatio>HyperVal);
            HypoInds=find(MethRatio<HypoVal);
            Hyper(loopd)=sum(PVecMSM(HyperInds));
            Hypo(loopd)=sum(PVecMSM(HypoInds));
        end

        ModelStruct.densityvals=MeanDens;
        ModelStruct.Hyper=Hyper;
        ModelStruct.Hypo=Hypo;
    end

    function [NetMeth,RelNetDeMeth,CoW,CoE,NormCoW,NormCoE]=GetParMetrics(p,DL)
        k1=p(1);
        k2=p(2);
        k3=p(3);
        k4=p(4);
        k5=p(5);
        k6=p(6);
        k7=p(7);
        k8=p(8);
        k9=p(9);
        k10=p(10);
        k11=p(11);
        k12=p(12);
        k13=p(13);
        k14=p(14);
        DLm=DL(6); %this assumes all m reactions (5-8) have same lengthscale
        DLu=DL(9); %assumes all u reactions (9-12) have same lengthscale
        k2t=k2+k14;
        k4t=k4+k13;
        NetMeth=k1*k3;
        RelNetDeMeth=k2t*k4t/(k1*k3);
        CoW=k6*k8/(k1*k3);
        CoE=k9*k11/(k2t*k4t);
        x=10; %inter-CpG distance
        fm=exp(-x/DLm);
        fu=exp(-x/DLu);
        CoWPart=fm*(k8*k3+k6*k1+fm*k6*k8);
        CoEPart=fu*(k9*k4t+k11*k2t+fu*k9*k11);
        NormCoW=CoWPart/(k1*k3+CoWPart);
        NormCoE=CoEPart/(k2*k4+CoEPart);

    end
end
