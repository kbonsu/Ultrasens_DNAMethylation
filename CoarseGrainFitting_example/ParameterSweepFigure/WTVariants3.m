function [n,K,ED50,flag,CoE,CoW,CoE2,CoW2] = WTVariants3(alpham,alphau,p,DistLength)

Nt=24; %cell cycle time
HyperVal=0.8; %cutoff for hypermeth
HypoVal=0.2; %cutoff for hypometh
NCpG=27;

k1 = p(1);
k2=p(2);
k3=p(3);
k4=p(4);
k6=p(6);
k8=p(8);
k9=p(9);
k11=p(11);
k13=p(13);
k14=p(14);

k2t=k2+k14;
k4t=k4+k13;

k9=k9*alphau;
k11=k11*alphau;
k6=k6*alpham;
k8=k8*alpham;
%2 different ways of quantifying relative CoE, CoW
CoE=k9*k11/(k2t*k4t);
CoW=k6*k8/(k1*k3);
CoW2=(k8*k3+k6*k1+k6*k8)/(k1*k3);
CoE2=(k9*k4t+k11*k2t+k9*k11)/(k2t*k4t);

Parameters=[k1,k2,k3,k4,0,k6,0,k8,k9,0,k11,0,k13,k14];

[n,K,ED50,flag]=JustParams(Parameters,DistLength);

    function [n,K,ED50] = JustPlot(Parameters,DistLength)
        [ModelStruct] = Run_Model(Parameters,DistLength);

        dens=ModelStruct.densityvals;
        Hyper=ModelStruct.Hyper;
        Hypo=ModelStruct.Hypo;
        Inter=ones(size(Hyper))-Hyper-Hypo;

        figure(1)
        plot(dens,Hyper,'-or')
        hold on
        plot(dens,Hypo,'-ob')
        plot(dens,Inter,'-ok')

        [n,K] = GetHillParameters(dens,Hypo);
        ED50 = GetED50(Hypo,Hyper,dens);
    end

    function [n,K,ED50,flag] = JustParams(Parameters,DistLength)
        [ModelStruct] = Run_Model(Parameters,DistLength);

        dens=ModelStruct.densityvals;
        Hyper=ModelStruct.Hyper;
        Hypo=ModelStruct.Hypo;
        Inter=ones(size(Hyper))-Hyper-Hypo;

        [n,K] = GetHillParameters(dens,Hypo);
        [ED50,flag] = GetED50(Hypo,Hyper,dens);
    end



    function [SSD] = PlotCompare(DataStruct,Parameters,DistLength)
        [ModelStruct] = Run_Model(Parameters,DistLength);
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

        print -dpng CompareFig
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
        ds=[2,3,4,5,6,7,8,9,11,13,17,26,200];
        %ds=[2,4,6,8,11,13,17,26,110];
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

    function [ED50,flag] = GetED50(Hypo,Hyper,xax)

        finex=linspace(min(xax),max(xax),300);
        Hypofine=interp1(xax,Hypo,finex);
        Hyperfine=interp1(xax,Hyper,finex);
        AAA = (Hypofine-Hyperfine);
        [v,inn]=min(abs(AAA));
        ED50=finex(inn);
        if v>0.05
            flag=1; %flag if it doesn't cross
        else
            flag=0;
        end
    end

    function [n,K] = GetHillParameters(x,y)
        %the log transformation method
        newy=log10(y./(1-y));
        newx=log10(x);

        %find the halfway point (similar to ED50)
        [v,i]=min(abs(y-0.5));

        %take just the neighboring points in vicinity
        fx=newx(2:end);%newx(i-3:i+3);
        fy=newy(2:end);%newy(i-3:i+3);

        %fit a line to the reduced, log-transformed data
        P = polyfit(fx,fy,1);

        %recover the Hill parameters (see equations at top for explanation)
        n=P(1);
        K=10^-(P(2)/P(1));
        H=(K^n)./(x.^n+K^n);
    end

    function EqualCoWrite = GetEqualCoWrite(d)
        CpGPositions=[1:d:NCpG*d];
        Densities=CpGDensities_Function(CpGPositions,50);
        MeanDens=mean(Densities);
        Dists=pdist(CpGPositions');
        DistMat=squareform(Dists);
        DistMat=DistMat+eye(NCpG)*1E6; %add a big number to the diagonal to silence self interaction
        ExpMat=exp(-IL*DistMat);
        %this calculates the average field
        f = 0.5*sum(ExpMat(:))/NCpG;
        %f=exp(-IL*d);
        NetDeMeth=(k2t+k9*f)*(k4t+k11*f);
        r1=lam*f^2;
        r2=f*k3+lam*f*k1;
        r3=-NetDeMeth+k1*k3;
        rr=[r1 r2 r3];
        rs=roots(rr);
        newk8=rs(rs>0);
        newk6=newk8*lam;
        EqualCoWrite=newk6(1)*newk8(1);
    end

end
