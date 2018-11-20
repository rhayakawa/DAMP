function [x_hat,arrMSE]=BODAMP(y,A,delta,arrP,arrR,nIteration,b)
%% BODAMP Bayes optimal DAMP 

  [M,N]=size(A);

  % x^{0}
  x=zeros(N,1);
  % z^{0}
  z=zeros(M,1);
  % theta^{0}
  theta=sqrt(arrP*(arrR.^(2))');
  % w^{0}
  w=zeros(N,1);
  arrMSE=zeros(1,nIteration);
  arrMSE(1)=norm(x-b)^(2)/N;
  for iterationIndex=2:nIteration
    % z^{t}
    z=y-A*x+1/delta*z*mean(drv_softThr_opt(w,arrP,arrR,theta,delta));
    % theta^{t}
    theta=norm(z)/sqrt(N); % approximation
    % x^{t+1}
    w=x+A'*z;
    x=softThr_opt(w,arrP,arrR,theta,delta);
    arrMSE(iterationIndex)=norm(x-b)^(2)/N;
  end
  x_hat=x;
  
end
