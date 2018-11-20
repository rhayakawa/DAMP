function DCTMat=myDCTMatrix(N)

  DCTMat=zeros(N,N);
  arrN=1:2:(2*N-1);
  arrK=0:(N-1);
  matN=ones(N,1)*arrN;
  matK=arrK'*ones(1,N);
  DCTMat=sqrt(2/N)*cos(pi/(2*N)*matN.*matK);
  DCTMat(1,:)=DCTMat(1,:)/sqrt(2);

end