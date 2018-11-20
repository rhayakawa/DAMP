function value=quantize(x,arrR)

  value=arrR(1)*ones(size(x));
  for index=2:length(arrR)
    value(x>(arrR(index-1)+arrR(index))/2)=arrR(index);
  end
  
end