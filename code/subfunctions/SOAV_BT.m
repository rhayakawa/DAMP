function x_hat=SOAV_BT(y,A,arrQ,arrR,alpha,gamma,nIteration)
%% SOAV SOAV optimization via Beck-Teboulle proximal gradient algorithm

  [M,N]=size(A);

  u=zeros(N,1);
  x=zeros(N,1);
  t=1;
  Ay=A'*y;
  for k=2:nIteration
    r=u-gamma^(-1)*alpha*(A'*(A*u)-Ay);
    x_former=x;
    x=softThr(r,gamma^(-1),arrQ,arrR);
    t_former=t;
    t=(1+sqrt(4*t_former^(2)+1))/2;
    u=x_former+(1+(t_former-1)/t)*(x-x_former);
  end
  x_hat=x;
  
end
