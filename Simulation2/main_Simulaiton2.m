clear all; close all,clc;

L = 10000;
Le=L/100;
p = 5;
type=1;

lamda_RLS = 0.99;
lamda_RLS_cta=0.99;
lamda_RBF_ATC=0.97;
lamda_RBF_CTA=0.97;
lamda_RBF_ATC0=0.977;
lamda_RBF_CTA0=0.977;
lamda_RBF_ATC2=0.963;
lamda_RBF_CTA2=0.963;

num=16;   %number of node
%% Initial
load('adjacencyMatrix.mat');

load('a.mat');
load('r.mat');
load('z.mat');
load('a0.mat');
load('r0.mat');
load('z0.mat');
load('a2.mat');
load('r2.mat');
load('z2.mat');

UU=[];
for i=1:num
    UU=[UU randn(p,L)];
end
wo1 = randn(p,1);
wo = [ kron(wo1, ones(1,L)) ];
tic
P=500; %% Independent Simulation
m0=40;
m2=100;
m=50;

for mm = 1 : P

    un1=zeros(p,num);
    dd=zeros(num,L);
    Fai=zeros(p,num);
    Fai_cta=Fai;
    FaiR=zeros(p,num);
    FaiR_cta=Fai;
    Fai_RMC=Fai;
    Fai_RLS=Fai;
    Fai_MCC=Fai;
    Fai_RBF0=Fai;
    Fai_RBF=Fai;
    Fai_RBF_CTA0=Fai;
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

    %% Iteration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ATC
    w_ATC=randn(p,num);
    w_ATC_RLS=w_ATC;
    w_ATC_RLS_RBF0=w_ATC;
    w_ATC_RLS_RBF2=w_ATC;
    w_ATC_RLS_RBF=w_ATC;

    for cnt=1:num
        Pn_RLS_BPD(:,:,cnt) = eye(p);
    end

    Pn_RLS=Pn_RLS_BPD;
    Pn_RMC=Pn_RLS_BPD;
    Pn_RBF0=Pn_RLS_BPD;
    Pn_RBF2=Pn_RLS_BPD;
    Pn_RBF=Pn_RLS_BPD;
    %% CTA
    w_CTA=w_ATC;
    w_CTA_RLS=w_ATC;
    w_CTA_RLS_RBF0=w_ATC;
    w_CTA_RLS_RBF2=w_ATC;
    w_CTA_RLS_RBF=w_ATC;
    w_CTA_RMC=w_ATC;
    w_CTA_MEE=w_ATC;
    w_M_RLS=w_ATC;
    w_M_CTA_RLS=w_ATC;

    Pn_RLS_cta=Pn_RLS_BPD;
    Pn_RMC_cta=Pn_RLS_BPD;
    Pn_RBF_CTA0=Pn_RLS_BPD;
    Pn_RBF_CTA2=Pn_RLS_BPD;
    Pn_RBF_CTA=Pn_RLS_BPD;
    for ii = 1 : L
        AddR=zeros(p,num);
        AddR_QPD=zeros(p,num);
        Add_ATC_RLS=zeros(p,num);
        Add_ATC_RMC=zeros(p,num);
        Add_ATC_MCC=zeros(p,num);
        Add_ATC_RBF0=zeros(p,num);
        Add_ATC_RBF2=zeros(p,num);
        Add_ATC_RBF=zeros(p,num);
        dn = dd(:,ii);
        for cnt=1:num
            Err_ATC_RLS(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS(:,cnt));
            Err_ATC_RLS_RBF0(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_RBF0(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_RBF0(:,cnt));
            Err_ATC_RLS_RBF2(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_RBF2(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_RBF2(:,cnt));
            Err_ATC_RLS_RBF(mm,ii,cnt) = (wo(:,ii)  -  w_ATC_RLS_RBF(:,cnt))' *  (wo(:,ii)  - w_ATC_RLS_RBF(:,cnt));
            un1(:,cnt)=UU(:,(cnt-1)*L+ii);
            %% ATC-RLS
            ek_ATC_RLS(cnt)=dn(cnt)-w_ATC_RLS(:,cnt)'*un1(:,cnt);
            kn_RLS(:,:,cnt) = Pn_RLS(:,:,cnt) * un1(:,cnt) / (lamda_RLS + un1(:,cnt)' * Pn_RLS(:,:,cnt) * un1(:,cnt)) ;
            Pn_RLS(:,:,cnt) = 1/lamda_RLS * ( Pn_RLS(:,:,cnt) - kn_RLS(:,:,cnt) *un1(:,cnt)' * Pn_RLS(:,:,cnt));
            F_ATC_RLS(:,cnt) = w_ATC_RLS(:,cnt) +kn_RLS(:,:,cnt) * ek_ATC_RLS(cnt);
            %% ATC-RLS-RBF0
            ek_Rand0(cnt) = dn(cnt) - w_ATC_RLS_RBF0(:,cnt)' * un1(:,cnt);
            if ii<Le
                Fai_RBF0(:,cnt)=F_ATC_RLS(:,cnt);
            else
                R11=0;
                R22=0;
                for jj = 1 : m0
                    R11=R11+r0(cnt,jj) * a0(cnt,jj)^2* exp( -( ek_Rand0(cnt)-c0(cnt,jj) )^2 * a0(cnt,jj)^2) ;
                    R22=R22+r0(cnt,jj) * a0(cnt,jj)^2* exp( -( ek_Rand0(cnt)-c0(cnt,jj) )^2 * a0(cnt,jj)^2) *c0(cnt,jj);
                end
                kn_RBF0(:,:,cnt)  = Pn_RBF0(:,:,cnt)  * un1(:,cnt) / ( lamda_RBF_ATC0+ R11*un1(:,cnt)' * Pn_RBF0(:,:,cnt) * un1(:,cnt) );
                Pn_RBF0(:,:,cnt) = 1/lamda_RBF_ATC0 * ( Pn_RBF0(:,:,cnt)  - R11*kn_RBF0(:,:,cnt) *un1(:,cnt)' * Pn_RBF0(:,:,cnt) );
                Fai_RBF0(:,cnt) = w_ATC_RLS_RBF0(:,cnt) + kn_RBF0(:,:,cnt) * (R11*ek_Rand0(cnt)-R22);
            end
            %% ATC-RLS-RBF2
            ek_Rand2(cnt) = dn(cnt) - w_ATC_RLS_RBF2(:,cnt)' * un1(:,cnt);
            if ii<Le
                Fai_RBF2(:,cnt)=F_ATC_RLS(:,cnt);
            else
                R1=0;
                R2=0;
                for jj = 1 : m2
                    R1=R1+r2(cnt,jj) * a2(cnt,jj)^2* exp( -( ek_Rand2(cnt)-c2(cnt,jj) )^2 * a2(cnt,jj)^2) ;
                    R2=R2+r2(cnt,jj) * a2(cnt,jj)^2* exp( -( ek_Rand2(cnt)-c2(cnt,jj) )^2 * a2(cnt,jj)^2) *c2(cnt,jj);
                end
                kn_RBF2(:,:,cnt)  = Pn_RBF2(:,:,cnt)  * un1(:,cnt) / ( lamda_RBF_ATC2+ R1*un1(:,cnt)' * Pn_RBF2(:,:,cnt) * un1(:,cnt) );
                Pn_RBF2(:,:,cnt) = 1/lamda_RBF_ATC2 * ( Pn_RBF2(:,:,cnt)  - R1*kn_RBF2(:,:,cnt) *un1(:,cnt)' * Pn_RBF2(:,:,cnt) );
                Fai_RBF2(:,cnt) = w_ATC_RLS_RBF2(:,cnt) + kn_RBF2(:,:,cnt) * (R1*ek_Rand2(cnt)-R2);
            end
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

                Add_ATC_RLS(:,cnt)=C(cnt,jj)*F_ATC_RLS(:,jj)+ Add_ATC_RLS(:,cnt);
                Add_ATC_RBF0(:,cnt)=C(cnt,jj)*Fai_RBF0(:,jj)+ Add_ATC_RBF0(:,cnt);
                Add_ATC_RBF2(:,cnt)=C(cnt,jj)*Fai_RBF2(:,jj)+ Add_ATC_RBF2(:,cnt);
                Add_ATC_RBF(:,cnt)=C(cnt,jj)*Fai_RBF(:,jj)+ Add_ATC_RBF(:,cnt);
            end
            w_ATC_RLS(:,cnt) =Add_ATC_RLS(:,cnt);
            w_ATC_RLS_RBF(:,cnt) =Add_ATC_RBF(:,cnt);
            w_ATC_RLS_RBF0(:,cnt) =Add_ATC_RBF0(:,cnt);
            w_ATC_RLS_RBF2(:,cnt) =Add_ATC_RBF2(:,cnt);
        end
        %% CTA
        AddR_cta=zeros(p,num);
        Add_RLS_cta=zeros(p,num);
        Add_RMC_cta=zeros(p,num);
        Add_RLS_RBF_CTA0=zeros(p,num);
        Add_RLS_RBF_CTA2=zeros(p,num);
        Add_RLS_RBF_CTA=zeros(p,num);
        dn = dd(:,ii);
        for cnt=1:num
            for jj=1:num
                Add_RLS_cta(:,cnt)=C(cnt,jj)*w_CTA_RLS(:,jj)+ Add_RLS_cta(:,cnt);
                Add_RLS_RBF_CTA0(:,cnt)=C(cnt,jj)*w_CTA_RLS_RBF0(:,jj)+ Add_RLS_RBF_CTA0(:,cnt);
                Add_RLS_RBF_CTA2(:,cnt)=C(cnt,jj)*w_CTA_RLS_RBF2(:,jj)+ Add_RLS_RBF_CTA2(:,cnt);
                Add_RLS_RBF_CTA(:,cnt)=C(cnt,jj)*w_CTA_RLS_RBF(:,jj)+ Add_RLS_RBF_CTA(:,cnt);
            end
            Fai_RLS(:,cnt) =Add_RLS_cta(:,cnt);
            if ii<Le
                Fai_RBF_CTA0(:,cnt)=Fai_RLS(:,cnt);
                Fai_RBF_CTA2(:,cnt)=Fai_RLS(:,cnt);
                Fai_RBF_CTA(:,cnt)=Fai_RLS(:,cnt);
            else
                Fai_RBF_CTA0(:,cnt) =Add_RLS_RBF_CTA0(:,cnt);
                Fai_RBF_CTA2(:,cnt) =Add_RLS_RBF_CTA2(:,cnt);
                Fai_RBF_CTA(:,cnt) =Add_RLS_RBF_CTA(:,cnt);
            end
        end
        for cnt=1:num
            Err_CTA_RLS(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS(:,cnt));
            Err_CTA_RLS_RBF0(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_RBF0(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_RBF0(:,cnt));
            Err_CTA_RLS_RBF2(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_RBF2(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_RBF2(:,cnt));
            Err_CTA_RLS_RBF(mm,ii,cnt) = (wo(:,ii)  -  w_CTA_RLS_RBF(:,cnt))' *  (wo(:,ii)  - w_CTA_RLS_RBF(:,cnt));
            un1(:,cnt)=UU(:,(cnt-1)*L+ii);
            %% CTA-RLS
            ek_CTA_RLS(cnt)=dn(cnt)-Fai_RLS(:,cnt)'*un1(:,cnt);
            kn_RLS_cta(:,:,cnt)= Pn_RLS_cta(:,:,cnt) * un1(:,cnt) / (lamda_RLS_cta + un1(:,cnt)' * Pn_RLS_cta(:,:,cnt) * un1(:,cnt)) ;
            Pn_RLS_cta(:,:,cnt) = 1/lamda_RLS_cta * ( Pn_RLS_cta(:,:,cnt) - kn_RLS_cta(:,:,cnt) *un1(:,cnt)' * Pn_RLS_cta(:,:,cnt));
            w_CTA_RLS(:,cnt) = Fai_RLS(:,cnt) +kn_RLS_cta(:,:,cnt) * ek_CTA_RLS(cnt);
            %% CTA-RLS-RBF0
            ek_Rand_cta0(cnt) = dn(cnt) - Fai_RBF_CTA0(:,cnt)' * un1(:,cnt);
            R1=0;
            R2=0;
            for jj = 1 : m0
                R1=R1+r0(cnt,jj) * a0(cnt,jj)^2* exp( -( ek_Rand_cta0(cnt)-c0(cnt,jj) )^2 * a0(cnt,jj)^2) ;
                R2=R2+r0(cnt,jj) * a0(cnt,jj)^2* exp( -( ek_Rand_cta0(cnt)-c0(cnt,jj) )^2 * a0(cnt,jj)^2) *c0(cnt,jj);
            end
            kn_RBF_cta0(:,:,cnt)  = Pn_RBF_CTA0(:,:,cnt)  * un1(:,cnt) / ( lamda_RBF_CTA0+ R1*un1(:,cnt)' * Pn_RBF_CTA0(:,:,cnt) * un1(:,cnt) );
            Pn_RBF_CTA0(:,:,cnt) = 1/lamda_RBF_CTA0 * ( Pn_RBF_CTA0(:,:,cnt)  - R1*kn_RBF_cta0(:,:,cnt) *un1(:,cnt)' * Pn_RBF_CTA0(:,:,cnt) );
            w_CTA_RLS_RBF0(:,cnt) = Fai_RBF_CTA0(:,cnt) + kn_RBF_cta0(:,:,cnt) * (R1*ek_Rand_cta0(cnt)-R2);
            %% CTA-RLS-RBF2
            ek_Rand_cta2(cnt) = dn(cnt) - Fai_RBF_CTA2(:,cnt)' * un1(:,cnt);
            R1=0;
            R2=0;
            for jj = 1 : m2
                R1=R1+r2(cnt,jj) * a2(cnt,jj)^2* exp( -( ek_Rand_cta2(cnt)-c2(cnt,jj) )^2 * a2(cnt,jj)^2) ;
                R2=R2+r2(cnt,jj) * a2(cnt,jj)^2* exp( -( ek_Rand_cta2(cnt)-c2(cnt,jj) )^2 * a2(cnt,jj)^2) *c2(cnt,jj);
            end
            kn_RBF_cta2(:,:,cnt)  = Pn_RBF_CTA2(:,:,cnt)  * un1(:,cnt) / ( lamda_RBF_CTA2+ R1*un1(:,cnt)' * Pn_RBF_CTA2(:,:,cnt) * un1(:,cnt) );
            Pn_RBF_CTA2(:,:,cnt) = 1/lamda_RBF_CTA2 * ( Pn_RBF_CTA2(:,:,cnt)  - R1*kn_RBF_cta2(:,:,cnt) *un1(:,cnt)' * Pn_RBF_CTA2(:,:,cnt) );
            w_CTA_RLS_RBF2(:,cnt) = Fai_RBF_CTA2(:,cnt) + kn_RBF_cta2(:,:,cnt) * (R1*ek_Rand_cta2(cnt)-R2);
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
    disp(mm);
end

Err_ATC_RLS1=mean(Err_ATC_RLS,3);
Err_ATC_RLS_RBF00=mean(Err_ATC_RLS_RBF0,3);
Err_ATC_RLS_RBF22=mean(Err_ATC_RLS_RBF2,3);
Err_ATC_RLS_RBF1=mean(Err_ATC_RLS_RBF,3);

Err_CTA_RLS1=mean(Err_CTA_RLS,3);
Err_CTA_RLS_RBF00=mean(Err_CTA_RLS_RBF0,3);
Err_CTA_RLS_RBF22=mean(Err_CTA_RLS_RBF2,3);
Err_CTA_RLS_RBF1=mean(Err_CTA_RLS_RBF,3);
for cnt=1:num
    Err_TH_MSE(cnt) = mean(Err_ATC_RLS(:,L,cnt));
    Err_TH_CTA_RLS(cnt) = mean(Err_CTA_RLS(:,L,cnt));
    Err_TH_ATC_RBF0(cnt) = mean(Err_ATC_RLS_RBF0(:,L,cnt));
    Err_TH_CTA_RBF0(cnt) = mean(Err_CTA_RLS_RBF0(:,L,cnt));
    Err_TH_ATC_RBF2(cnt) = mean(Err_ATC_RLS_RBF2(:,L,cnt));
    Err_TH_CTA_RBF2(cnt) = mean(Err_CTA_RLS_RBF2(:,L,cnt));
    Err_TH_ATC_RBF(cnt) = mean(Err_ATC_RLS_RBF(:,L,cnt));
    Err_TH_CTA_RBF(cnt) = mean(Err_CTA_RLS_RBF(:,L,cnt));
end

toc
gap=400;
x= 1:gap:L;

Err1=mean(Err_CTA_RLS1);
Err2=mean(Err_ATC_RLS1);

Err9=mean(Err_CTA_RLS_RBF1);
Err10=mean(Err_ATC_RLS_RBF1);
Err11=mean(Err_CTA_RLS_RBF00);
Err12=mean(Err_ATC_RLS_RBF00);
Err13=mean(Err_CTA_RLS_RBF22);
Err14=mean(Err_ATC_RLS_RBF22);

figure(1),hold on;

plot(10* log10(Err1),'LineWidth',2,'color','m');
plot(10* log10(Err2),'LineWidth',2,'color','b');

plot(10*log10(Err11),'LineWidth',2,'color','c');
plot(10*log10(Err12),'LineWidth',2,'color','g');

plot(10*log10(Err9),'LineWidth',2,'color',[0.3 0.6 0.9]);
plot(10*log10(Err10),'LineWidth',2,'color',[0.9 0.6 0.3]);

plot(10*log10(Err13),'LineWidth',2,'color','k');
plot(10*log10(Err14),'LineWidth',2,'color','y');

legend('CTA-RLS','ATC-RLS','Proposed CTA-RBF(M=30)','Proposed ATC-RBF(M=30)','Proposed CTA-RBF(M=50)','Proposed ATC-RBF(M=50)','Proposed CTA-RBF(M=100)','Proposed ATC-RBF(M=100)');xlabel('Iterations');ylabel('MSD (dB)');
ylim([-31,20]);
xlim([1,L]);
box on;grid on;
hold off

figure(2),hold on
plot(10*log10((Err_TH_CTA_RLS)),'-m*','LineWidth',1.5);
plot(10*log10((Err_TH_MSE)),'-bp','LineWidth',1.5);

plot(10*log10(( Err_TH_CTA_RBF0)),'-s','LineWidth',1.5,'color','c');
plot(10*log10(( Err_TH_ATC_RBF0)),'-v','LineWidth',1.5,'color','g');

plot(10*log10(( Err_TH_CTA_RBF)),'-h','LineWidth',1.5,'color',[0.3 0.6 0.9]);
plot(10*log10(( Err_TH_ATC_RBF)),'-d','LineWidth',1.5,'color',[0.9 0.6 0.3]);

plot(10*log10(( Err_TH_CTA_RBF2)),'-^','LineWidth',1.5,'color','k');
plot(10*log10(( Err_TH_ATC_RBF2)),'-x','LineWidth',1.5,'color','y');

legend('CTA-RLS','ATC-RLS','Proposed CTA-RBF(M=30)','Proposed ATC-RBF(M=30)','Proposed CTA-RBF(M=50)','Proposed ATC-RBF(M=50)','Proposed CTA-RBF(M=100)','Proposed ATC-RBF(M=100)');
xlabel('Node number,n');ylabel('Steady State MSD (dB)');
ylim([-31,-12]);
xlim([1,16]);
legend('NumColumns', 2)
box on;grid on;
hold off

