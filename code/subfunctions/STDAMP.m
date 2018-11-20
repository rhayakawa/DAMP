function [x_hat,arrMSE]=STDAMP(y,A,delta,arrQ,arrP,arrR,nIteration,b)
%% STDAMP soft thresholding DAMP 

  [M,N]=size(A);

  % x^{0}
  x=zeros(N,1);
  x=ones(N,1)*(arrP*arrR');
  % z^{0}
  z=zeros(M,1);
  % theta^{0}
  theta=sqrt(arrP*(arrR.^(2))');
  % w^{0}
  w=zeros(N,1);
  % array for MSE
  arrMSE=zeros(1,nIteration);
  arrMSE(1)=norm(x-b)^(2)/N;
  for iterationIndex=2:nIteration
    % z^{t}
    z=y-A*x+1/delta*z*mean(drv_softThr(w,theta/sqrt(delta),arrQ,arrR));
    % theta^{t}
    theta=norm(z)/sqrt(N); 
    % x^{t+1}
    w=x+A'*z;
    x=softThr(w,theta/sqrt(delta),arrQ,arrR);
    arrMSE(iterationIndex)=norm(x-b)^(2)/N;
  end
  x_hat=x;
  
end
