%%-------------------------------------------
% Simulation for MSE of DAMP (vs iteration)
%%-------------------------------------------

clear;
addpath subfunctions;

N=1000; % number of unknown variables
delta=0.8; % measurement ratio
M=round(N*delta); % number of observations
nIteration=50; % number of iterations
nSample=10; % number of samples
SNR=30; % signal-to-noise ratio

%% distribution of the unknown variable (possibly sparse vector with (-1, 0, 1))
p0=0.2
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
  arrG2=-1/2*(1-3*p)*(normpdf(Q)+Q.*normcdf(Q))+1/2*(1-p0)*Q;
  [~,index_Q2opt]=min(abs(arrG2));
  Q2_opt=Q(index_Q2opt)
end
arrQ_opt=[-1e10 Q2_opt Q3_opt 1e10];

%% cumulative distribution
L=length(arrR);
matOne=ones(L,L);
arrCDF=arrP*triu(matOne);

%% noise variance
sigma2_w=N*(1-p0)/(M*10^(SNR/10));

%% State evolution
disp('State evolution');
arrSE_STDAMP=zeros(1,nIteration);
arrSE_STDAMP(1)=arrP*(arrR.^(2))';
arrSE_BODAMP=zeros(1,nIteration);
arrSE_BODAMP(1)=arrP*(arrR.^(2))';
for iterationIndex=2:nIteration
  disp(['iterationIndex=' num2str(iterationIndex)]);
  arrSE_STDAMP(iterationIndex)=Psi(arrSE_STDAMP(iterationIndex-1)+delta*sigma2_w,delta,arrP,arrR,arrQ_opt);
  arrSE_BODAMP(iterationIndex)=Psi_opt(arrSE_BODAMP(iterationIndex-1)+delta*sigma2_w,delta,arrP,arrR);
end

rng('shuffle');

arrSumMSE_STDAMP=zeros(1,nIteration);
arrSumMSE_BODAMP=zeros(1,nIteration);
for sampleIndex=1:nSample
  disp(['  sampleIndex=' num2str(sampleIndex)]);
  % unknown discrete-valued vector
  b_rand=rand(N,1);
  b=ones(N,1)*arrR(1);
  for valueIndex=2:L
    b(b_rand>=arrCDF(valueIndex-1))=arrR(valueIndex);
  end
  % measurement matrix
  A=randn(M,N)/sqrt(M);
  % additive noise vector
  v=randn(M,1)*sqrt(sigma2_w);
  % linear measurements
  y=A*b+v;
  
  %% soft thresholding DAMP
  [~,arrMSE_STDAMP]=STDAMP(y,A,delta,arrQ_opt,arrP,arrR,nIteration,b);
  arrSumMSE_STDAMP=arrSumMSE_STDAMP+arrMSE_STDAMP;
  
  %% Bayes optimal DAMP
  [~,arrMSE_BODAMP]=BODAMP(y,A,delta,arrP,arrR,nIteration,b);
  arrSumMSE_BODAMP=arrSumMSE_BODAMP+arrMSE_BODAMP;
  
end
arrMSE_STDAMP=arrSumMSE_STDAMP/nSample;
arrMSE_BODAMP=arrSumMSE_BODAMP/nSample;

%% plot results
close all;
figure;
h=semilogy(1:nIteration,arrMSE_STDAMP,'^','LineWidth',1,'MarkerSize',8);
hold on;
h=semilogy(1:nIteration,arrMSE_BODAMP,'o','LineWidth',1,'MarkerSize',8);
grid on;
h=semilogy(1:nIteration,arrSE_STDAMP,'--k','LineWidth',1,'MarkerSize',8);
h=semilogy(1:nIteration,arrSE_BODAMP,'-k','LineWidth',1,'MarkerSize',8);

xlabel('number of iterations');
ylabel('MSE');
objLegend=legend('STDAMP (empirical)','BODAMP (empirical)','STDAMP (theoretical)','BODAMP (theoretical)');
objLegend.Interpreter='latex';
objLegend.Location='northeast';
objLegend.FontSize=16;
fig=gca;
fig.FontSize=18;
fig.TickLabelInterpreter='latex';
fig.XLabel.Interpreter='latex';
fig.YLabel.Interpreter='latex';
axis([0 nIteration 1e-4 1]);

saveas(h, ['MSE.eps'], 'epsc');
