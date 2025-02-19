clear all; close all,clc;

L = 10000;
Le=L/100;
p = 5;
type=1;
MEE_Length=10;
lamda_RLS = 0.99;
lamda_RLS_cta=0.99;
lamda_RMCC = 0.98;
lamda_RMCC_cta = 0.98;
lamda_RMEE=0.995;
lamda_New_ATC=0.993;
lamda_New_CTA=0.993;
lamda_RBF_ATC=0.97;
lamda_RBF_CTA=0.97;

orders=2;
orders1=4;
num=16;   %number of node
%% Initial
load('adjacencyMatrix.mat');
load('a.mat');
load('r.mat');
load('z.mat');
MIU_EM=zeros(num,orders);
SIGMA_EM=zeros(num,orders);
WEIGHT_EM=zeros(num,orders);

UU=[];
for i=1:num
    UU=[UU randn(p,L)];

end
wo1 = randn(p,1);
wo = [ kron(wo1, ones(1,L)) ];
tic
P=50; %% Independent Simulation
m=100;
%% Mean by EM
MIU_EM(1,:)=[0,0];
MIU_EM(2,:)=[0,0];
MIU_EM(3,:)=[0,0];
MIU_EM(4,:)=[0,0];

MIU_EM(5,:)=[-0.27,0.31];
MIU_EM(6,:)=[-0.18,0.31];
MIU_EM(7,:)=[-0.12,0.23];
MIU_EM(8,:)=[-0.16,0.27];

MIU_EM(9,:)=[-0.5,0.5];
MIU_EM(10,:)=[-1,1];
MIU_EM(11,:)=[-1.5,1.5];
MIU_EM(12,:)=[-2,2];

MIU_EM(13,:)=[-2.8,3.3];
MIU_EM(14,:)=[-4.5,5.3];
MIU_EM(15,:)=[-3.3,3.1];
MIU_EM(16,:)=[-4.6,5.1];
%% Sigma by EM
SIGMA_EM(1,:)=[0.1,5];
SIGMA_EM(2,:)=[0.5,9.2];
SIGMA_EM(3,:)=[0.1,11.3];
SIGMA_EM(4,:)=[0.5,10.4];

SIGMA_EM(5,:)=[1.1,2.1];
SIGMA_EM(6,:)=[1.2,2.3];
SIGMA_EM(7,:)=[0.6,1];
SIGMA_EM(8,:)=[0.9,1.8];

SIGMA_EM(9,:)=[0.05,0.05];
SIGMA_EM(10,:)=[0.25,0.25];
SIGMA_EM(11,:)=[0.5,0.5];
SIGMA_EM(12,:)=[0.75,0.75];

SIGMA_EM(13,:)=[2.2,0.1];
SIGMA_EM(14,:)=[2.1,0.1];
SIGMA_EM(15,:)=[2.9,0.1];
SIGMA_EM(16,:)=[3.2,0.1];
%% Weight by EM
WEIGHT_EM(1,:)=[0.9,0.1];
WEIGHT_EM(2,:)=[0.89,0.11];
WEIGHT_EM(3,:)=[0.84,0.16];
WEIGHT_EM(4,:)=[0.9,0.1];

WEIGHT_EM(5,:)=[0.5,0.5];
WEIGHT_EM(6,:)=[0.5,0.5];
WEIGHT_EM(7,:)=[0.53,0.47];
WEIGHT_EM(8,:)=[0.6,0.4];

WEIGHT_EM(9,:)=[0.52,0.48];
WEIGHT_EM(10,:)=[0.5,0.5];
WEIGHT_EM(11,:)=[0.53,0.47];
WEIGHT_EM(12,:)=[0.46,0.54];

WEIGHT_EM(13,:)=[0.51,0.49];
WEIGHT_EM(14,:)=[0.53,0.47];
WEIGHT_EM(15,:)=[0.5,0.5];
WEIGHT_EM(16,:)=[0.45,0.55];

MIU_Rand=MIU_EM;
SIGMA_Rand=SIGMA_EM;
WEIGHT_Rand=WEIGHT_EM;

for mm = 1 : P

    un1=zeros(p,num);
    dd=zeros(num,L);
    Fai=zeros(p,num);
    Fai_QPD=zeros(p,num);
    Fai_cta_QPD=Fai;
    Fai_cta=Fai;
    FaiR=zeros(p,num);
    FaiR_cta=Fai;
    Fai_RMC=Fai;
    Fai_RLS=Fai;
    Fai_MCC=Fai;
    Fai_RBF=Fai;
    Fai_RBF_CTA=Fai;
    %% Noise
    %% Impulse Noise
    for n=1:L
        Epsilon=rand;
        if Epsilon<0.9
            VV1(n)=normrnd(0,0.1);
        else
            VV1(n)=normrnd(0,5);
        end
    end

    for n=1:L
        Epsilon=rand;
        if Epsilon<0.9
            VV2(n)=normrnd(0,0.5);
        else
            VV2(n)=normrnd(0,5);
        end
    end

    for n=1:L
        Epsilon=rand;
        if Epsilon<0.9
            VV3(n)=normrnd(0,0.1);
        else
            VV3(n)=normrnd(0,10);
        end
    end

    for n=1:L
        Epsilon=rand;
        if Epsilon<0.9
            VV4(n)=normrnd(0,0.5);
        else
            VV4(n)=normrnd(0,10);
        end
    end

    v(1,:)=VV1;
    v(2,:)=VV2;
    v(3,:)=VV3;
    v(4,:)=VV4;

    %% Skewed Noise
    v(5,:)=raylrnd(8,1,L)*0.1;
    v(6,:)=gamrnd(3,2,1,L)*0.1;
    v(7,:)=poissrnd(6,1,L)*0.1;
    v(8,:)=chi2rnd(6,1,L)*0.1;
    %% Uniform Noise
    v(9,:)=2*rand(1,L)-1;
    v(10,:)=4*rand(1,L)-2;
    v(11,:)=6*rand(1,L)-3;
    v(12,:)=8*rand(1,L)-4;
    %% Multi-peak Noise
    for n=1:L
        Epsilon=rand;
        if Epsilon<0.5
            VV13(n)=normrnd(-3,2);
        else
            VV13(n)=normrnd(3,0.1);
        end
    end

    mu=-3;
    sigma=4;
    b=sigma/sqrt(2);
    mu1=3;
    sigma1=0.01;
    b1=sigma1/sqrt(2);
    weight_laplace=[0.5 0.5];

    for n=1:L
        Epsilon=rand;
        if Epsilon<weight_laplace(1)
            aa=rand(1,1)-0.5;
            VV14(n)=mu-b*sign(aa).*log(1-2*abs(aa));
        else
            aa=rand(1,1)-0.5;
            VV14(n)=mu1-b1*sign(aa).*log(1-2*abs(aa));
        end
    end

    for n=1:L
        Epsilon=rand;
        if Epsilon<0.5
            VV15(n)=normrnd(-5,3);
        else
            VV15(n)=normrnd(5,0.1);
        end
    end

    mu=-5;
    sigma=9;
    b=sigma/sqrt(2);
    mu1=5;
    sigma1=0.01;
    b1=sigma1/sqrt(2);
    weight_laplace=[0.5 0.5];

    for n=1:L
        Epsilon=rand;
        if Epsilon<weight_laplace(1)
            aa=rand(1,1)-0.5;
            VV16(n)=mu-b*sign(aa).*log(1-2*abs(aa));
        else
            aa=rand(1,1)-0.5;
            VV16(n)=mu1-b1*sign(aa).*log(1-2*abs(aa));
        end
    end
    v(13,:)=VV13;
    v(14,:)=VV14;
    v(15,:)=VV15;
    v(16,:)=VV16;

    for i=1:num
        vv(i,:)=v(i,:)-mean(v(i,:));
        Sig(i)=var (vv(i,1:L));
    end
    for i=1:num
        for ii = 1 : L
            dd(i,ii) = wo(:,ii)' * UU(:,(i-1)*L+ii) + vv(i,ii);
        end
    end
    C=normalizedAdjacencyMatrix;
    sigma1=2; %DRMCC
    sigma2=4; %DRMEE
    %% Iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ATC
    w_ATC=randn(p,num);
    w_ATC_RLS=w_ATC;
    w_ATC_RLS_BPD=w_ATC;
    w_ATC_RLS_QPD=w_ATC;
    w_ATC_RLS_RBF=w_ATC;
    w_MEE=w_ATC;
    w_CTA_MEE=w_ATC;
    w_ATC_RMC=w_ATC;
    w_ATC_MCC=w_ATC;
    for cnt=1:num
        Pn_RLS_BPD(:,:,cnt) = eye(p);
    end

    Pn_RLS=Pn_RLS_BPD;
    Pn_RMC=Pn_RLS_BPD;
    Pn_RBF=Pn_RLS_BPD;
    %% CTA
    w_CTA=w_ATC;
    w_CTA_RLS=w_ATC;
    w_CTA_RLS_BPD=w_ATC;
    w_CTA_RLS_QPD=w_ATC;
    w_CTA_RLS_RBF=w_ATC;
    w_CTA_RMC=w_ATC;
    w_CTA_MEE=w_ATC;
    w_M_RLS=w_ATC;
    w_M_CTA_RLS=w_ATC;
    Pn_RLS_BPD_cta=Pn_RLS_BPD;
    Pn_RLS_cta=Pn_RLS_BPD;
    Pn_RMC_cta=Pn_RLS_BPD;
    Pn_RLS_QPD_cta=Pn_RLS_BPD;
    Pn_RBF_CTA=Pn_RLS_BPD;
    for ii = 1 : L
        AddR=zeros(p,num);
        AddR_QPD=zeros(p,num);
        Add_ATC_RLS=zeros(p,num);
        Add_ATC_RMC=zeros(p,num);
        Add_ATC_MCC=zeros(p,num);
        Add_ATC_RBF=zeros(p,num);
        dn = dd(:,ii);
        for cnt=1:num
            Err_ATC_RLS(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS(:,cnt));
            Err_ATC_RMC(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RMC(:,cnt))' *  (wo(:,ii)  - w_ATC_RMC(:,cnt));
            Err_ATC_RLS_BPD(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_BPD(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_BPD(:,cnt));
            Err_ATC_RLS_RBF(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_RBF(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_RBF(:,cnt));
            un1(:,cnt)=UU(:,(cnt-1)*L+ii);
            %% ATC-RLS
            ek_ATC_RLS(cnt)=dn(cnt)-w_ATC_RLS(:,cnt)'*un1(:,cnt);
            kn_RLS(:,:,cnt) = Pn_RLS(:,:,cnt) * un1(:,cnt) / (lamda_RLS + un1(:,cnt)' * Pn_RLS(:,:,cnt) * un1(:,cnt)) ;
            Pn_RLS(:,:,cnt) = 1/lamda_RLS * ( Pn_RLS(:,:,cnt) - kn_RLS(:,:,cnt) *un1(:,cnt)' * Pn_RLS(:,:,cnt));
            F_ATC_RLS(:,cnt) = w_ATC_RLS(:,cnt) +kn_RLS(:,:,cnt) * ek_ATC_RLS(cnt);
            %% ATC-RMC
            ek_ATC_RMC(cnt)=dn(cnt)-w_ATC_RMC(:,cnt)'*un1(:,cnt);
            kn_RMC(:,:,cnt) = Pn_RMC(:,:,cnt) * un1(:,cnt) / ( exp(ek_ATC_RMC(cnt)^2/2/sigma1^2)*lamda_RMCC + un1(:,cnt)' * Pn_RMC(:,:,cnt) * un1(:,cnt) );
            Pn_RMC(:,:,cnt) = 1/lamda_RMCC * ( Pn_RMC(:,:,cnt) - kn_RMC(:,:,cnt) *un1(:,cnt)' * Pn_RMC(:,:,cnt));
            F_ATC_RMC(:,cnt) = w_ATC_RMC(:,cnt) +kn_RMC(:,:,cnt) * ek_ATC_RMC(cnt);
            %% ATC-RLS-BPD
            ek_Rand1(cnt)=dn(cnt)-w_ATC_RLS_BPD(:,cnt)'*un1(:,cnt);
            for k=1:orders
                P_Rand(ii,k+orders*(cnt-1) )=exp(-1*(ek_Rand1(cnt)-MIU_Rand(cnt,k))^2/(2*SIGMA_Rand(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand(cnt,k));
            end
            for k=1:orders
                V_Rand(ii,k+orders*(cnt-1) )=WEIGHT_Rand(cnt,k)*P_Rand(ii,k+orders*(cnt-1))/(WEIGHT_Rand(cnt,:)*P_Rand(ii,orders*(cnt-1)+1:orders*(cnt-1)+orders)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders
                R1=R1+V_Rand(ii,k+orders*(cnt-1))/SIGMA_Rand(cnt,k)^2;
                R2=R2+V_Rand(ii,k+orders*(cnt-1))*MIU_Rand(cnt,k)/SIGMA_Rand(cnt,k)^2;
            end
            kn(:,:,cnt) = Pn_RLS_BPD(:,:,cnt) * un1(:,cnt) / ( lamda_New_ATC+ R1*un1(:,cnt)' * Pn_RLS_BPD(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_BPD(:,:,cnt) = 1/lamda_New_ATC * ( Pn_RLS_BPD(:,:,cnt) - R1*kn(:,:,cnt) * un1(:,cnt)' * Pn_RLS_BPD(:,:,cnt));
            FaiR(:,cnt) = w_ATC_RLS_BPD(:,cnt) +kn(:,:,cnt) * (R1*ek_Rand1(cnt)-R2);
            %% ATC-RLS-RBF
            ek_Rand(cnt) = dn(cnt) - w_ATC_RLS_RBF(:,cnt)' * un1(:,cnt);
            if ii<Le
                Fai_RBF(:,cnt)=F_ATC_RLS(:,cnt);
            else
                R1=0;
                R2=0;
                for jj = 1 : m
                    R1=R1+r(cnt,jj) * a(cnt,jj)^2* exp( -( ek_Rand(cnt)-c(cnt,jj) )^2 * a(cnt,jj)^2) ;
                    R2=R2+r(cnt,jj) * a(cnt,jj)^2* exp( -( ek_Rand(cnt)-c(cnt,jj) )^2 * a(cnt,jj)^2) *c(cnt,jj);
                end
                kn_RBF(:,:,cnt)  = Pn_RBF(:,:,cnt)  * un1(:,cnt) / ( lamda_RBF_ATC+ R1*un1(:,cnt)' * Pn_RBF(:,:,cnt) * un1(:,cnt) );
                Pn_RBF(:,:,cnt) = 1/lamda_RBF_ATC * ( Pn_RBF(:,:,cnt)  - R1*kn_RBF(:,:,cnt) *un1(:,cnt)' * Pn_RBF(:,:,cnt) );
                Fai_RBF(:,cnt) = w_ATC_RLS_RBF(:,cnt) + kn_RBF(:,:,cnt) * (R1*ek_Rand(cnt)-R2);
            end
        end
        for cnt=1:num
            for jj=1:num
                AddR(:,cnt)=C(cnt,jj)*FaiR(:,jj)+AddR(:,cnt);
                AddR_QPD(:,cnt)=C(cnt,jj)* Fai_QPD(:,jj)+AddR_QPD(:,cnt);
                Add_ATC_RMC(:,cnt)=C(cnt,jj)*F_ATC_RMC(:,jj)+ Add_ATC_RMC(:,cnt);
                Add_ATC_RLS(:,cnt)=C(cnt,jj)*F_ATC_RLS(:,jj)+ Add_ATC_RLS(:,cnt);
                Add_ATC_RBF(:,cnt)=C(cnt,jj)*Fai_RBF(:,jj)+ Add_ATC_RBF(:,cnt);
            end
            w_ATC_RLS_BPD(:,cnt) = AddR(:,cnt);
            w_ATC_RMC(:,cnt) =Add_ATC_RMC(:,cnt);
            w_ATC_RLS(:,cnt) =Add_ATC_RLS(:,cnt);
            w_ATC_RLS_RBF(:,cnt) =Add_ATC_RBF(:,cnt);
        end

        %% CTA
        AddR_cta=zeros(p,num);
        AddR_QPD_cta=zeros(p,num);
        Add_RLS_cta=zeros(p,num);
        Add_RMC_cta=zeros(p,num);
        Add_RLS_RBF_CTA=zeros(p,num);
        dn = dd(:,ii);
        for cnt=1:num
            for jj=1:num
                AddR_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS_BPD(:,jj)+AddR_cta(:,cnt);
                Add_RMC_cta(:,cnt)=C(cnt,jj)*w_CTA_RMC(:,jj)+ Add_RMC_cta(:,cnt);
                Add_RLS_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS(:,jj)+ Add_RLS_cta(:,cnt);
                Add_RLS_RBF_CTA(:,cnt)=C(cnt,jj)*w_CTA_RLS_RBF(:,jj)+ Add_RLS_RBF_CTA(:,cnt);
            end
            FaiR_cta(:,cnt) = AddR_cta(:,cnt);
            Fai_cta_QPD(:,cnt) = AddR_QPD_cta(:,cnt);
            Fai_RMC(:,cnt) =Add_RMC_cta(:,cnt);
            Fai_RLS(:,cnt) =Add_RLS_cta(:,cnt);
            if ii<Le
                Fai_RBF_CTA(:,cnt)=Fai_RLS(:,cnt);
            else
                Fai_RBF_CTA(:,cnt) =Add_RLS_RBF_CTA(:,cnt);
            end
        end
        for cnt=1:num
            Err_CTA_RLS(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS(:,cnt));
            Err_CTA_RMC(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RMC(:,cnt))' *  (wo(:,ii)  - w_CTA_RMC(:,cnt));
            Err_CTA_RLS_BPD(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_BPD(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_BPD(:,cnt));
            Err_CTA_RLS_RBF(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_RBF(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_RBF(:,cnt));
            un1(:,cnt)=UU(:,(cnt-1)*L+ii);
            %% CTA-RLS
            ek_CTA_RLS(cnt)=dn(cnt)-Fai_RLS(:,cnt)'*un1(:,cnt);
            kn_RLS_cta(:,:,cnt)= Pn_RLS_cta(:,:,cnt) * un1(:,cnt) / (lamda_RLS_cta + un1(:,cnt)' * Pn_RLS_cta(:,:,cnt) * un1(:,cnt)) ;
            Pn_RLS_cta(:,:,cnt) = 1/lamda_RLS_cta * ( Pn_RLS_cta(:,:,cnt) - kn_RLS_cta(:,:,cnt) *un1(:,cnt)' * Pn_RLS_cta(:,:,cnt));
            w_CTA_RLS(:,cnt) = Fai_RLS(:,cnt) +kn_RLS_cta(:,:,cnt) * ek_CTA_RLS(cnt);
            %% CTA-RMC
            ek_CTA_RMC(:,:,cnt)=dn(cnt)-Fai_RMC(:,cnt)'*un1(:,cnt);
            kn_RMC_cta(:,:,cnt) = Pn_RMC_cta(:,:,cnt) * un1(:,cnt) / ( exp(ek_CTA_RMC(cnt)^2/2/sigma1^2)*lamda_RMCC_cta + un1(:,cnt)' * Pn_RMC_cta(:,:,cnt) * un1(:,cnt) );
            Pn_RMC_cta(:,:,cnt) = 1/lamda_RMCC_cta * ( Pn_RMC_cta(:,:,cnt) - kn_RMC_cta(:,:,cnt)*un1(:,cnt)' * Pn_RMC_cta(:,:,cnt));
            w_CTA_RMC(:,cnt) = Fai_RMC(:,cnt) +kn_RMC_cta(:,:,cnt) * ek_CTA_RMC(cnt);
            %% CTA-RLS-BPD
            ek_Rand1_cta(cnt)=dn(cnt)-FaiR_cta(:,cnt)'*un1(:,cnt);
            for k=1:orders
                P_Rand_cta(ii,k+orders*(cnt-1) )=exp(-1*(ek_Rand1_cta(cnt)-MIU_Rand(cnt,k))^2/(2*SIGMA_Rand(cnt,k)^2))/(sqrt(2*pi)*SIGMA_Rand(cnt,k));
            end
            for k=1:orders
                V_Rand_cta(ii,k+orders*(cnt-1) )=WEIGHT_Rand(cnt,k)*P_Rand_cta(ii,k+orders*(cnt-1))/(WEIGHT_Rand(cnt,:)*P_Rand_cta(ii,orders*(cnt-1)+1:orders*(cnt-1)+orders)'+eps);
            end
            R1=0;
            R2=0;
            for k=1:orders
                R1=R1+V_Rand_cta(ii,k+orders*(cnt-1))/SIGMA_Rand(cnt,k)^2;
                R2=R2+V_Rand_cta(ii,k+orders*(cnt-1))*MIU_Rand(cnt,k)/SIGMA_Rand(cnt,k)^2;
            end
            kn_cta(:,:,cnt) = Pn_RLS_BPD_cta(:,:,cnt) * un1(:,cnt) / ( lamda_New_CTA+ R1*un1(:,cnt)' * Pn_RLS_BPD_cta(:,:,cnt) * un1(:,cnt) );
            Pn_RLS_BPD_cta(:,:,cnt)= 1/lamda_New_CTA * ( Pn_RLS_BPD_cta(:,:,cnt) - R1*kn_cta(:,:,cnt)* un1(:,cnt)' * Pn_RLS_BPD_cta(:,:,cnt));
            w_CTA_RLS_BPD(:,cnt) = FaiR_cta(:,cnt) +kn_cta(:,:,cnt) * (R1*ek_Rand1_cta(cnt)-R2);
            %% CTA-RLS-RBF
            ek_Rand_cta(cnt) = dn(cnt) - Fai_RBF_CTA(:,cnt)' * un1(:,cnt);

            R1=0;
            R2=0;
            for jj = 1 : m
                R1=R1+r(cnt,jj) * a(cnt,jj)^2* exp( -( ek_Rand_cta(cnt)-c(cnt,jj) )^2 * a(cnt,jj)^2) ;
                R2=R2+r(cnt,jj) * a(cnt,jj)^2* exp( -( ek_Rand_cta(cnt)-c(cnt,jj) )^2 * a(cnt,jj)^2) *c(cnt,jj);
            end
            kn_RBF_cta(:,:,cnt)  = Pn_RBF_CTA(:,:,cnt)  * un1(:,cnt) / ( lamda_RBF_CTA+ R1*un1(:,cnt)' * Pn_RBF_CTA(:,:,cnt) * un1(:,cnt) );
            Pn_RBF_CTA(:,:,cnt) = 1/lamda_RBF_CTA * ( Pn_RBF_CTA(:,:,cnt)  - R1*kn_RBF_cta(:,:,cnt) *un1(:,cnt)' * Pn_RBF_CTA(:,:,cnt) );
            w_CTA_RLS_RBF(:,cnt) = Fai_RBF_CTA(:,cnt) + kn_RBF_cta(:,:,cnt) * (R1*ek_Rand_cta(cnt)-R2);

        end
    end
    %% ATC RMEE, algorithm
    for cnt=1:num
        PL(:,:,cnt) = eye(p)*1;
    end
    eL = zeros( MEE_Length, 1);
    PHI = zeros( MEE_Length,  MEE_Length);
    for ii =  MEE_Length : L
        Add_RMEE=zeros(p,num);
        F_RMEE= w_M_RLS;
        for cnt=1:num
            Err_RLS_MEE(mm,ii,cnt) = (wo(:,ii)  - w_M_RLS(:,cnt))' * (wo(:,ii)  - w_M_RLS(:,cnt));
            for jj = 1 :  MEE_Length
                eL(jj) = dd(cnt,ii - MEE_Length + jj) - w_M_RLS(:,cnt)' * UU(:,(cnt-1)*L+ii - MEE_Length + jj);
            end
            u0 = UU(:,(cnt-1)*L+ii);
            e0 = eL(MEE_Length);
            phi0 = 0;
            for kk = 2: MEE_Length
                ek = eL( MEE_Length-kk+1);
                phi0 = phi0 + lamda_RMEE^(kk-1) * exp(- (e0-ek).^2/sigma2^2/2);
            end
            KL(:,:,cnt)= PL(:,:,cnt) * u0 * inv ( lamda_RMEE^2 / phi0 + u0' * PL(:,:,cnt) * u0 );
            PL(:,:,cnt) = 1/lamda_RMEE^2 * (PL(:,:,cnt) - KL(:,:,cnt) * u0' * PL(:,:,cnt));
            F_RMEE(:,cnt) =F_RMEE(:,cnt) +KL(:,:,cnt) * e0;
        end
        for cnt=1:num
            for jj=1:num
                Add_RMEE(:,cnt)=C(cnt,jj)*F_RMEE(:,jj)+Add_RMEE(:,cnt);
            end
            w_M_RLS(:,cnt) = Add_RMEE(:,cnt);
        end
    end
    %% CTA RMEE, algorithm
    for cnt=1:num
        PL1(:,:,cnt) = eye(p)*1;
    end
    eL1 = zeros( MEE_Length, 1);
    PHI = zeros( MEE_Length,  MEE_Length);
    for ii =  MEE_Length : L
        Add_CTA_RMEE=zeros(p,num);

        for cnt=1:num
            for jj=1:num
                Add_CTA_RMEE(:,cnt)=C(cnt,jj)*w_M_CTA_RLS(:,jj)+Add_CTA_RMEE(:,cnt);
            end
            w_M_CTA_RLS(:,cnt) = Add_CTA_RMEE(:,cnt);
        end
        for cnt=1:num
            Err_RLS_CTA_MEE(mm,ii,cnt) = (wo(:,ii)  - w_M_CTA_RLS(:,cnt))' * (wo(:,ii)  - w_M_CTA_RLS(:,cnt));
            for jj = 1 :  MEE_Length
                eL1(jj) = dd(cnt,ii - MEE_Length + jj) - w_M_CTA_RLS(:,cnt)' * UU(:,(cnt-1)*L+ii - MEE_Length + jj);
            end
            u0 = UU(:,(cnt-1)*L+ii);
            e0 = eL1(MEE_Length);
            phi0 = 0;
            for kk = 2: MEE_Length
                ek = eL1( MEE_Length-kk+1);
                phi0 = phi0 + lamda_RMEE^(kk-1) * exp(- (e0-ek).^2/sigma2^2/2);
            end
            KL1(:,:,cnt) = PL1(:,:,cnt) * u0 * inv ( lamda_RMEE^2 / phi0 + u0' * PL1(:,:,cnt) * u0 );
            PL1(:,:,cnt) = 1/lamda_RMEE^2 * (PL1(:,:,cnt) - KL1(:,:,cnt)* u0' * PL1(:,:,cnt));
            w_M_CTA_RLS(:,cnt) =w_M_CTA_RLS(:,cnt) +KL1(:,:,cnt) * e0;
        end
    end
    disp(mm);
end


Err_ATC_RLS1=mean(Err_ATC_RLS,3);
Err_ATC_RMC1=mean(Err_ATC_RMC,3);
Err_ATC_RLS_BPD1=mean(Err_ATC_RLS_BPD,3);
Err_ATC_RMEE1=mean(Err_RLS_MEE,3);
Err_ATC_RLS_RBF1=mean(Err_ATC_RLS_RBF,3);

Err_CTA_RLS1=mean(Err_CTA_RLS,3);
Err_CTA_RMC1=mean(Err_CTA_RMC,3);
Err_CTA_RLS_BPD1=mean(Err_CTA_RLS_BPD,3);
Err_CTA_RMEE1=mean( Err_RLS_CTA_MEE,3);
Err_CTA_RLS_RBF1=mean(Err_CTA_RLS_RBF,3);
for cnt=1:num
    Err_TH_MCC(cnt) =mean(Err_ATC_RMC(:,L,cnt));
    Err_TH_CTA_RMC(cnt) =mean(Err_CTA_RMC(:,L,cnt));
    Err_TH_MSE(cnt) = mean(Err_ATC_RLS(:,L,cnt));
    Err_TH_CTA_RLS(cnt) = mean(Err_CTA_RLS(:,L,cnt));
    Err_TH_ATC_RMEE(cnt) = mean(Err_RLS_MEE(:,L,cnt));
    Err_TH_RMEE(cnt) = mean(Err_RLS_CTA_MEE(:,L,cnt));
    Err_TH_RLS_New(cnt) = mean(Err_ATC_RLS_BPD(:,L,cnt));
    Err_TH_CTA_RLS_BPD(cnt) = mean(Err_CTA_RLS_BPD(:,L,cnt));
    Err_TH_ATC_RBF(cnt) = mean(Err_ATC_RLS_RBF(:,L,cnt));
    Err_TH_CTA_RBF(cnt) = mean(Err_CTA_RLS_RBF(:,L,cnt));
end

toc

Err1=mean(Err_CTA_RLS1);
Err2=mean(Err_ATC_RLS1);
Err3=mean(Err_CTA_RMC1);
Err4=mean(Err_ATC_RMC1);
Err5=mean(Err_CTA_RMEE1);
Err6=mean(Err_ATC_RMEE1);
Err7=mean(Err_CTA_RLS_BPD1);
Err8=mean(Err_ATC_RLS_BPD1);
Err9=mean(Err_CTA_RLS_RBF1);
Err10=mean(Err_ATC_RLS_RBF1);

figure(1),hold on;
plot(10* log10(Err1),'LineWidth',0.1,'color','m');
plot(10* log10(Err2),'LineWidth',0.1,'color','b');
plot(10* log10(Err3),'LineWidth',0.1,'color',[0.6 0.1 0.5]);
plot(10* log10(Err4),'LineWidth',0.1,'color','r');
plot(10* log10(Err5),'LineWidth',0.1,'color','y');
plot(10* log10(Err6),'LineWidth',0.1,'color','k');
plot(10*log10(Err7),'LineWidth',0.1,'color','g');
plot(10*log10(Err8),'LineWidth',0.1,'color','c');
plot(10*log10(Err9),'LineWidth',0.1,'color',[0.3 0.6 0.9]);
plot(10*log10(Err10),'LineWidth',0.1,'color',[0.9 0.6 0.3]);
legend('CTA-RLS','ATC-RLS','CTA-RMCC(\sigma=2)','ATC-RMCC(\sigma=2)','CTA-RMEE(\sigma=4,L=10)','ATC-RMEE(\sigma=4,L=10)','CTA-RMLF(K=2)','ATC-RMLF(K=2)','Proposed CTA-RBF(M=50)','Proposed ATC-RBF(M=50)');xlabel('Iterations');ylabel('MSD (dB)');
ylim([-31,35]);
xlim([1,L]);
box on;grid on;
hold off

figure(2),hold on
plot(10*log10((Err_TH_CTA_RLS)),'-m*','LineWidth',1.5);
plot(10*log10((Err_TH_MSE)),'-bp','LineWidth',1.5);
plot(10*log10((Err_TH_CTA_RMC)),'-<','LineWidth',1.5,'color',[0.6 0.1 0.5]);
plot(10*log10((Err_TH_MCC)),'-rs','LineWidth',1.5);
plot(10*log10((Err_TH_RMEE)),'-yv','LineWidth',1.5);
plot(10*log10(( Err_TH_ATC_RMEE)),'-ko','LineWidth',1.5);
plot(10*log10((Err_TH_CTA_RLS_BPD)),'-g^','LineWidth',1.5);
plot(10*log10((Err_TH_RLS_New)),'-cx','LineWidth',1.5);
plot(10*log10(( Err_TH_CTA_RBF)),'-h','LineWidth',1.5,'color',[0.3 0.6 0.9]);
plot(10*log10(( Err_TH_ATC_RBF)),'-d','LineWidth',1.5,'color',[0.9 0.6 0.3]);
legend('CTA-RLS','ATC-RLS','CTA-RMCC(\sigma=2)','ATC-RMCC(\sigma=2)','CTA-RMEE(\sigma=4,L=10)','ATC-RMEE(\sigma=4,L=10)','CTA-RMLF(K=2)','ATC-RMLF(K=2)','Proposed CTA-RBF(M=50)','Proposed ATC-RBF(M=50)');
xlabel('Node number,n');ylabel('Steady State MSD (dB)');
ylim([-33,30]);
xlim([1,16]);
legend('NumColumns', 2)
box on;grid on;
hold off

