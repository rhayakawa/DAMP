function value=drv_softThr_opt(x,p,r,sigma,delta)
  X=x*ones(1,length(r));
  P=ones(length(x),1)*p;
  R=ones(length(x),1)*r;
  
  xr2_min=min(((X-R).^(2)),[],2);
  XR2_min=xr2_min*ones(1,length(r));
  E=exp(-1/2*delta/sigma^(2)*((X-R).^(2)-XR2_min));
  E(isnan(E))=1;
  p_e=sum(P.*E,2);
  p_r_e=sum(P.*R.*E,2);
  p_r2_e=sum(P.*(R.^(2)).*E,2);
  
  if sigma^(2)>1e-30
    value=delta/sigma^(2)*(p_r2_e.*p_e-p_r_e.^(2))./(p_e.^(2));
  else
    value=zeros(size(x)); % actually the value is Inf for the elements of r
  end
  value(isnan(value))=0;
end
