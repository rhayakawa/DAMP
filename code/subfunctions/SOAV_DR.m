function x_hat=SOAV_DR(y,A,AAinv,arrQ,arrR,lambda,gamma,nIteration)
%% SOAV SOAV optimization via Douglas-Rachford algorithm in the noise-free case

  [M,N]=size(A);

  x_ZF=A'*(AAinv*y);
  z_til=zeros(N,1);
  for k=2:nIteration
    z=z_til-A'*(AAinv*(A*z_til))+x_ZF;
    z_til=z_til+lambda*(softThr(2*z-z_til,gamma,arrQ,arrR)-z);
  end
  x_hat=z;
  
end
