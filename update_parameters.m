function [C, d, pi_new] = update_parameters(x, y, responsabilities)
% UPDATE_PARAMETERS -

  M=size(responsabilities,1);
  N=size(x,2);
  N_k=sum(responsabilities,2);
  pi_new(:)=N_k/N;
  theta=[x' ones(size(x',1),1)];
  for i=1:M
      Gamma=diag(sqrt(responsabilities(i,:)));
      pseudo=(Gamma*theta)'*(Gamma'*theta);
      temp=(Gamma*theta)'*(Gamma'*y');
      % temp=(pseudo\x)*Gamma;
      phi(:,i)=pseudo\temp;
      C(:,:,i)=phi(1:end-1,i);
      d(:,:,i)=phi(end,i);
  end
end
