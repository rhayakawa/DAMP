function value=Psi(sigma2,delta,arrP,arrR,arrU)
  
  sigma=sqrt(sigma2);
  L=length(arrR);
  T=zeros(L,L,L+1);
  for l=1:L
    T(l,:,:)=T(l,:,:)-sqrt(delta)/sigma*arrR(l);
  end
  for k=1:L
    T(:,k,:)=T(:,k,:)+sqrt(delta)/sigma*arrR(k);
  end
  for k_til=1:L+1
    T(:,:,k_til)=T(:,:,k_til)+arrU(k_til);
  end
  phi=normpdf(T);
  Phi=normcdf(T);
  
  arrPsi_l=zeros(1,L);
  for l=1:L
    for k=1:L
      arrPsi_l(l)=arrPsi_l(l)+(arrR(k)-arrR(l))^(2)*(Phi(l,k,k+1)-Phi(l,k,k));
    end
    for k=2:L
      arrPsi_l(l)=arrPsi_l(l)+sigma2/delta*(-T(l,k,k)*phi(l,k,k)+Phi(l,k,k)+2*arrU(k)*phi(l,k,k)+arrU(k)^(2)*Phi(l,k,k));
    end
    for k=1:(L-1)
      arrPsi_l(l)=arrPsi_l(l)+sigma2/delta*(T(l,k,k+1)*phi(l,k,k+1)-Phi(l,k,k+1)-2*arrU(k+1)*phi(l,k,k+1)-arrU(k+1)^(2)*Phi(l,k,k+1));
    end
  end
  
  value=arrPsi_l*arrP';
  
end
