function [p,SSD] = Fitting_Thurs(Par0,DL0,DataName)

%Par0 - initial parameters for fit
%DL0 - initial distance lengthscale parameters for fit
%DataName - points to right dataset to be fit
close all
HyperVal=0.8; %cutoff for hypermeth
HypoVal=0.2; %cutoff for hypometh
NCpG=27; %number of CpGs in the model


DataFile=['Save_' DataName '.mat'];
load(DataFile,'DataStruct')
k1=Par0(1);
k2=Par0(2);
k3=Par0(3);
k13=Par0(13);
k14=Par0(14);

p0(1)=Par0(4);
p0(2)=Par0(6);
p0(3)=Par0(8);
p0(4)=Par0(9);
p0(5)=Par0(11);

%add in a little randomness
p0=p0+1E-3*randn(size(p0));

%initialize the parameters
Parameters=[k1,k2,k3,p0(1),0,p0(2),0,p0(3),p0(4),0,p0(5),0,k13,k14];
DistLength=DL0;
fun=@(p)Fit_Model(p,DataStruct,Parameters)

mns=ones(size(p0))*1E-3;
mxs=[3,25,25,25,25,75,75];
options=optimoptions('fmincon','Display','iter');
nonlcon=[];

%run the fit
p = fmincon(fun,p0,[],[],[],[],mns,mxs,nonlcon,options);
Parameters=[k1,k2,k3,p(1),0,p(2),0,p(3),p(4),0,p(5),0,k13,k14]

tag=randi(100000); %generate a random integer tag. This will be a unique identifier for the output
%from the fitting
fn=['FitTh_' DataName '_' num2str(tag)];
ffn=[fn '.png']; %figure filename
afn=[fn '.mat']; %array filename

[SSD] = PlotCompare(DataStruct,Parameters,DistLength,ffn);
save(afn,'Parameters','DistLength','SSD')



    function [SSD] = PlotCompare(DataStruct,Parameters,DistLength,fn)
        [ModelStruct] = Run_Model_MorePoints(Parameters,DistLength);
        SSD = Get_Error(ModelStruct,DataStruct);

        dens=ModelStruct.densityvals;
        Hyper=ModelStruct.Hyper;
        Hypo=ModelStruct.Hypo;
        Inter=ones(size(Hyper))-Hyper-Hypo;

        figure(2)
        plot(dens,Hyper,'-or')
        hold on
        plot(dens,Hypo,'-ob')
        plot(dens,Inter,'-ok')

        dens=DataStruct.densityvals;
        Hyper=DataStruct.Hyper;
        Hypo=DataStruct.Hypo;
        Inter=ones(numel(Hyper),1)-Hyper-Hypo;

        plot(dens,Hyper,'--r')
        hold on
        plot(dens,Hypo,'--b')
        plot(dens,Inter,'--k')

        text(0.4,0.5,['SSD =' num2str(SSD)])

        print(fn,'-dpng')
    end


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
