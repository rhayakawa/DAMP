function eta=softThr(u,c,arrQ,arrR)

  L=length(arrR);
  eta=u-c*arrQ(1);
  for l=1:L
    index=(u>arrR(l)+c*arrQ(l));
    eta(index)=arrR(l);
    index2=(u>arrR(l)+c*arrQ(l+1));
    eta(index2)=u(index2)-c*arrQ(l+1);
  end
  
end
