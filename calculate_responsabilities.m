function responsabilities = calculate_responsabilities(x, y, C, d, Sigma, pi)
% CALCULATE_RESPONSABILITIES -

  M=size(d,3);
  N=size(x,3);
  respNum =zeros(size(x,2),M);
  %= calculate numerator for each mode
  for i=1:M
    w=pi(i)*sqrt(Sigma(:,:,i));
    lktemp=diag(y-[C(:,:,i)'*x+d(:,:,i)]);
    SigmaInvMat=kron(eye(N),Sigma(:,:,i)\1);
    lk=diag(-1/2*lktemp*SigmaInvMat*lktemp');

    %= logsumexp trick see: https://www.cs.bgu.ac.il/~cv202/wiki.files/logsumexp_trick.pdf
    a=max(lk)+log(w);
    likelihood_mod=exp(lk+log(w)-a);
    % likelihood=mvnpdf(y',(C(:,i)'*x+d(:,i))',Sigma(:,:,i));
    respNum(:,i)=likelihood_mod;
  end
  normalizing_constant=sum(respNum ,2);

  responsabilities =(respNum./normalizing_constant)';
end
