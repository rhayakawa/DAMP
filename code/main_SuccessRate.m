%%------------------------------------------------------------
% Simulation for success rate of DAMP (vs measurement ratio)
%%------------------------------------------------------------

clear;
addpath subfunctions;

N=100; % number of unknown variables
nIteration=300; % number of iterations
nSample=30; % number of samples

%% distribution of the unknown variable (binary with (0, 1))
p1=0.1;
p2=1-p1;
arrP=[p1 p2];
arrR=[0 1];

%% array for observation ratio
% threshold: delta=0.256 for p1=0.1
arrDelta=[0.06:0.02:0.7];

%% optimal Q
Q=-20:0.001:20;
arrG2=p1*(-normpdf(Q)+Q.*(1-normcdf(Q)))+p2*(normpdf(Q)+Q.*normcdf(Q));
[G2_min,index_opt]=min(abs(arrG2));
Q2_opt=Q(index_opt);
arrQ_opt=[-1e10 Q2_opt 1e10];
arrQ_orig=[-1 p1-p2 1];

%% cumulative distribution
L=length(arrR);
matOne=ones(L,L);
arrCDF=arrP*triu(matOne);

rng('shuffle');

arrNumSuccess_STDAMP=zeros(1,length(arrDelta));
arrNumSuccess_BODAMP=zeros(1,length(arrDelta));
arrNumSuccess_SOAV_orig=zeros(1,length(arrDelta));
arrNumSuccess_SOAV_opt=zeros(1,length(arrDelta));
for deltaIndex=1:length(arrDelta)
  delta=arrDelta(deltaIndex);
  disp(['delta=' num2str(delta)]);
  M=round(N*delta); % measurement ratio
  
  for sampleIndex=1:nSample
    %% problem settings
    % unknown discrete-valued vector
    b_rand=rand(N,1);
    b=ones(N,1)*arrR(1);
    for valueIndex=2:L
      b(b_rand>=arrCDF(valueIndex-1))=arrR(valueIndex);
    end
    % measurement matrix
    A=randn(M,N)/sqrt(M);
    % linear measurements
    y=A*b;
    
    %% soft thresholding DAMP
    [x_hat_STDAMP,~]=STDAMP(y,A,delta,arrQ_opt,arrP,arrR,nIteration,b);
    numError_STDAMP=nnz(quantize(x_hat_STDAMP,arrR)-b);
    if numError_STDAMP==0
      arrNumSuccess_STDAMP(deltaIndex)=arrNumSuccess_STDAMP(deltaIndex)+1;
    end

    %% Bayes optimal DAMP
    [x_hat_BODAMP,~]=BODAMP(y,A,delta,arrP,arrR,nIteration,b);
    numError_BODAMP=nnz(quantize(x_hat_BODAMP,arrR)-b);
    if numError_BODAMP==0
      arrNumSuccess_BODAMP(deltaIndex)=arrNumSuccess_BODAMP(deltaIndex)+1;
    end

    %% SOAV optimization (original)
    gamma=1;
    lambda=1.9;
    AAinv=(A*A')^(-1);
    x_hat_SOAV_orig=SOAV_DR(y,A,AAinv,arrQ_orig,arrR,lambda,gamma,nIteration);
    numError_SOAV_orig=nnz(quantize(x_hat_SOAV_orig,arrR)-b);
    if numError_SOAV_orig==0
      arrNumSuccess_SOAV_orig(deltaIndex)=arrNumSuccess_SOAV_orig(deltaIndex)+1;
    end

    %% SOAV optimization (proposed)
    x_hat_SOAV_opt=SOAV_DR(y,A,AAinv,arrQ_opt,arrR,lambda,gamma,nIteration);
    numError_SOAV_opt=nnz(quantize(x_hat_SOAV_opt,arrR)-b);
    if numError_SOAV_opt==0
      arrNumSuccess_SOAV_opt(deltaIndex)=arrNumSuccess_SOAV_opt(deltaIndex)+1;
    end

  end
end
arrSuccessRate_STDAMP=arrNumSuccess_STDAMP/nSample;
arrSuccessRate_BODAMP=arrNumSuccess_BODAMP/nSample;
arrSuccessRate_SOAV_orig=arrNumSuccess_SOAV_orig/nSample;
arrSuccessRate_SOAV_opt=arrNumSuccess_SOAV_opt/nSample;

%% Display results
close all;
figure;
h=plot(arrDelta,arrSuccessRate_STDAMP,'--^','LineWidth',1,'MarkerSize',8);
hold on;
h=plot(arrDelta,arrSuccessRate_BODAMP,'--^','LineWidth',1,'MarkerSize',8,'MarkerFaceColor',[0.8500 0.3250 0.0980]);
h=plot(arrDelta,arrSuccessRate_SOAV_orig,'-s','LineWidth',1,'MarkerSize',8,'Color',[0.4660 0.6740 0.1880]);
h=plot(arrDelta,arrSuccessRate_SOAV_opt,'-s','LineWidth',1,'MarkerSize',8,'Color',[0.3010 0.7450 0.9330],'MarkerFaceColor',[0.3010 0.7450 0.9330]);
grid on;
x_plot=[1 1];
y_plot=[0 1];
delta=0.256; % threshold for p1=0.1;
deltaThr=delta.*x_plot;
plot(deltaThr,y_plot,'-k','LineWidth',1);
xlabel('$\Delta$','Interpreter','latex');
ylabel('success rate','Interpreter','latex');
objLegend=legend('STDAMP','BODAMP','SOAV (original)','SOAV (proposed)');
objLegend.Interpreter='latex';
objLegend.Location='southeast';
objLegend.FontSize=16;
axis([0 1 0 1]);
fig=gca;
fig.XTick=0:0.1:1;
fig.FontSize=20;
fig.TickLabelInterpreter='latex';

saveas(h, 'SuccessRate.eps', 'epsc');
