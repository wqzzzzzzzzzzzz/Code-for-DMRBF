clear all, close all ,clc
%%
%% Noise
L=10000;
num=16;
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
end
for k=1:num
    [fx(k,:),xi(k,:)] = ksdensity(v(k,:),'NumPoints',L);
end
% fx=fx/sum(fx).*100;

%% RBF NN fitting
m=100;
Len_xi = L;
Num = L;
for k=1:num

    a(k,:) =rand(1,m)*2-1;
    c(k,:)=randn(1,m);
    r(k,:) = rand(1,m);
    z = zeros(m,1);
    y(k,:) = zeros(1,Len_xi);
    mu = 0.005;
    for kk =1: Num
        ii = randi(Len_xi);
        xk = xi(k,ii);
        dk = fx(k,ii);
        for jj=1:m
            z(jj) = exp(-(xk-c(k,jj))^2 * a(k,jj)^2);
        end
        y(k,ii) = r(k,:)* z ;
        e = dk - y(k,ii);
        r(k,:) = r(k,:)  + mu * e * z';
        for jj=1:m
            c(k,jj) = c(k,jj) + mu * e *r(k,jj) * z(jj) * (xk-c(k,jj));
            a(k,jj) =a(k,jj) - mu * e* r(k,jj) * z(jj) * (xk-c(k,jj))^2 *a(k,jj);
        end
        Err(k,kk) = norm(fx(k,:)-y(k,:));
    end

    for ii=1 : Len_xi
        xk = xi(k,ii);
        for jj=1:m
            z(jj) = exp(-(xk-c(k,jj))^2 * a(k,jj)^2);
        end
        y_rbf(k,ii) = r(k,:)* z ;
    end
    figure(k)
  
    hold on;
    plot(xi(k,:),fx(k,:),'r');
    plot(xi(k,:),y_rbf(k,:),'k');
    legend('拟合','RBF');

    disp(k)
end

save('a.mat', 'a');
save('r.mat', 'r');
save('z.mat', 'c');