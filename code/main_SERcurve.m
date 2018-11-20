%%-------------------------------------------------
% Simulation for SER performance of DAMP (vs SNR)
%%-------------------------------------------------

clear;
addpath subfunctions;

N=1000; % number of unknown variables
delta=0.8; % measurement ratio
M=round(N*delta); % number of observations
nIteration=200; % number of iterations
nSample=10; % number of samples
arrSNR=0:2.5:30; % array for signal-to-noise ratio

%% measurement matrix (i.i.d. or correlated or DCT)
strMat='i.i.d.';
% strMat='correlated';
% strMat='DCT';

%% distribution of the unknown variable (possibly sparse vector with (-1, 0, 1))
p0=0.2; % probability of 0
arrP=[(1-p0)/2 p0 (1-p0)/2];
arrR=[-1 0 1];
if p0<=1/3
  Q3_opt=0;
  Q2_opt=0;
else
  Q=-20:0.001:20;
  arrG3=1/2*(1-3*p0)*(normpdf(Q)+Q.*normcdf(Q))+p0*Q;
  [~,index_Q3opt]=min(abs(arrG3));
  Q3_opt=Q(index_Q3opt)
  arrG2=-1/2*(1-3*p0)*(normpdf(Q)+Q.*normcdf(Q))+1/2*(1-p0)*Q;
  [~,index_Q2opt]=min(abs(arrG2));
  Q2_opt=Q(index_Q2opt)
end
arrQ_opt=[-1e10 Q2_opt Q3_opt 1e10];
arrQ_orig=[-1 -p0 p0 1];

%% cumulative distribution
L=length(arrR);
matOne=ones(L,L);
arrCDF=arrP*triu(matOne);

%% parameter for correlated channel
if strcmp(strMat,'correlated')
  lambda=1;
  d_t=0.5*lambda;
  d_r=0.5*lambda;
  Delta_t=zeros(N,N);
  for i=1:N
      Delta_t(i,:)=abs((i-1):(-1):(i-N));
  end
  Delta_r=zeros(M,M);
  for i=1:M
      Delta_r(i,:)=abs((i-1):(-1):(i-M));
  end
  Phi_t=besselj(0,Delta_t*2*pi*d_t/lambda);
  Phi_r=besselj(0,Delta_r*2*pi*d_r/lambda);
end

%% DCT matrix
if strcmp(strMat,'DCT')
  D=myDCTMatrix(N);
end

rng('shuffle');

arrSumSERcurve_STDAMP=zeros(1,length(arrSNR));
arrSumSERcurve_BODAMP=zeros(1,length(arrSNR));
arrSumSERcurve_SOAV_orig=zeros(1,length(arrSNR));
arrSumSERcurve_SOAV_opt=zeros(1,length(arrSNR));
for sampleIndex=1:nSample
  disp(['  sampleIndex=' num2str(sampleIndex)]);
  %% problem settings
  % unknown discrete-valued vector
  b_rand=rand(N,1);
  b=ones(N,1)*arrR(1);
  for valueIndex=2:L
    b(b_rand>=arrCDF(valueIndex-1))=arrR(valueIndex);
  end
  % measurement matrix
  if strcmp(strMat,'i.i.d.')
    A=randn(M,N)/sqrt(M);
  elseif strcmp(strMat,'correlated')
    A=randn(M,N)/sqrt(M);
    A=Phi_r^(1/2)*A*Phi_t^(1/2);
  elseif strcmp(strMat,'DCT')
    I=eye(N);
    S=I(randperm(N,M),:); % selection matrix
    A=S*D;
    A=A/sqrt(delta); %scaling for DMAP
  end

  for SNRIndex=1:length(arrSNR)
    % additive noise vector
    SNR=arrSNR(SNRIndex);
    sigma2_w=N*(1-p0)/(M*10^(SNR/10));
    if strcmp(strMat,'DCT')
      % scaling for DAMP
      sigma2_w=sigma2_w/delta;
    end
    v=randn(M,1)*sqrt(sigma2_w);
    % linear measurements
    y=A*b+v;
    
    %% soft thresholding DAMP
    [x_hat_STDAMP,~]=STDAMP(y,A,delta,arrQ_opt,arrP,arrR,nIteration,b);
    SER_STDAMP=nnz(quantize(x_hat_STDAMP,arrR)-b)/N;
    arrSumSERcurve_STDAMP(SNRIndex)=arrSumSERcurve_STDAMP(SNRIndex)+SER_STDAMP;

    %% Bayes optimal DAMP
    [x_hat_BODAMP,~]=BODAMP(y,A,delta,arrP,arrR,nIteration,b);
    SER_BODAMP=nnz(quantize(x_hat_BODAMP,arrR)-b)/N;
    arrSumSERcurve_BODAMP(SNRIndex)=arrSumSERcurve_BODAMP(SNRIndex)+SER_BODAMP;

    %% SOAV optimization (original)
    alpha=10; % parameter
    gamma=alpha*norm(A)^(2);
    x_hat_SOAV_orig=SOAV_BT(y,A,arrQ_orig,arrR,alpha,gamma,nIteration);
    SER_SOAV_orig=nnz(quantize(x_hat_SOAV_orig,arrR)-b)/N;
    arrSumSERcurve_SOAV_orig(SNRIndex)=arrSumSERcurve_SOAV_orig(SNRIndex)+SER_SOAV_orig;

    %% SOAV optimization (proposed)
    x_hat_SOAV_opt=SOAV_BT(y,A,arrQ_opt,arrR,alpha,gamma,nIteration);
    SER_SOAV_opt=nnz(quantize(x_hat_SOAV_opt,arrR)-b)/N;
    arrSumSERcurve_SOAV_opt(SNRIndex)=arrSumSERcurve_SOAV_opt(SNRIndex)+SER_SOAV_opt;

  end
end
arrSERcurve_STDAMP=arrSumSERcurve_STDAMP/nSample;
arrSERcurve_BODAMP=arrSumSERcurve_BODAMP/nSample;
arrSERcurve_SOAV_orig=arrSumSERcurve_SOAV_orig/nSample;
arrSERcurve_SOAV_opt=arrSumSERcurve_SOAV_opt/nSample;

%% Display results
close all;
figure;
h=semilogy(arrSNR,arrSERcurve_STDAMP,'--^','LineWidth',1,'MarkerSize',8);
hold on;
h=semilogy(arrSNR,arrSERcurve_BODAMP,'--^','LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
h=semilogy(arrSNR,arrSERcurve_SOAV_orig,'-s','LineWidth',1,'MarkerSize',8,'Color',[0.4660 0.6740 0.1880]);
h=semilogy(arrSNR,arrSERcurve_SOAV_opt,'-s','LineWidth',1,'MarkerSize',8,'Color',[0.3010 0.7450 0.9330],'MarkerFaceColor',[0.3010 0.7450 0.9330]);
grid on;
xlabel('SNR (dB)','Interpreter','latex');
ylabel('SER','Interpreter','latex');
objLegend=legend('STDAMP','BODAMP','SOAV (original)','SOAV (proposed)');
objLegend.Interpreter='latex';
objLegend.Location='southwest';
objLegend.FontSize=16;
axis([0 arrSNR(length(arrSNR)) 1e-4 1]);
fig=gca;
fig.FontSize=18;
fig.TickLabelInterpreter='latex';

saveas(h, 'SERcurve.eps', 'epsc');
