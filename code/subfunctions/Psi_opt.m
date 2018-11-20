function value=Psi_opt(sigma2,delta,arrP,arrR)
    
  sigma=sqrt(sigma2);
  value=integral(@func,-Inf,Inf);
  
  function y=func(z)
      y=0;
      for i=1:length(arrP)
          y=y+arrP(i)*((softThr_opt(arrR(i)+sigma/sqrt(delta)*z,arrP,arrR,sigma,delta)-arrR(i)).^(2)).*normpdf(z);
      end
  end

end
