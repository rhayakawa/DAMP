function value=softThr_opt(x,p,r,sigma,delta)
  
  bRow=false;
  if isrow(x)
      bRow=true;
      x=x';
  end
  r=r(p~=0);
  p=p(p~=0);

  X=x*ones(1,length(r));
  P=ones(length(x),1)*p;
  R=ones(length(x),1)*r;
  
  xr2_min=min(((X-R).^(2)),[],2);
  XR2_min=xr2_min*ones(1,length(r));
  E=exp(-1/2*delta/sigma^(2)*((X-R).^(2)-XR2_min));
  E(isnan(E))=1;
  p_e=sum(P.*E,2);
  p_r_e=sum(P.*R.*E,2);
  
  value=p_r_e./p_e;
  if bRow
      value=value';
  end

end
