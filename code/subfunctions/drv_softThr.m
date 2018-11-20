function eta_drv=drv_softThr(u,c,arrQ,arrR)

  L=length(arrR);
  eta_drv=ones(size(u));
  for l=1:L
    index=(u>arrR(l)+c*arrQ(l));
    eta_drv(index)=0;
    index2=(u>arrR(l)+c*arrQ(l+1));
    eta_drv(index2)=1;
  end
  
end  