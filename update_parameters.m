function [Phi, pi_new] = update_parameters(X, Y, Responsabilities)
% UPDATE_PARAMETERS -
    modes=size(Responsabilities,1);
    [n N]=size(X);

    A=kron(ones(n,1),eye(n));
    B=kron(eye(N),ones(1,n));
    C=kron(repmat(eye(n),1,N),ones(n,1));
    theta=sparse([(A*X*B).*C;kron(repmat(eye(n),1,N),ones(1,1))]');

    N_k=sum(Responsabilities,2);
    pi_new(:)=N_k/N;
    for i=1:modes
        responsabilities=Responsabilities(i,:);
        resp2=cellfun(@(x) x*eye(n),mat2cell(responsabilities',ones(1,N)),'UniformOutput',0);
        Gamma=sqrt(sparse(blkdiag(resp2{:})));

        Phi(i,:)=-((Gamma*theta)\(Gamma*Y(:)));
    end

end
