function [Phi, pi_new, Sigma] = update_parameters(X, Y, OldPhi, Responsabilities)
% UPDATE_PARAMETERS -
    modes=size(Responsabilities,1);
    [n N]=size(X);

    A=kron(ones(n,1),eye(n));
    B=kron(eye(N),ones(1,n));
    C=kron(repmat(eye(n),1,N),ones(n,1));
    theta=sparse([(A*X*B).*C;kron(repmat(eye(n),1,N),ones(1,1))]');

    y_cell = cellfun(@(x) x',mat2cell(Y',ones(1,N))','UniformOutput',0);
    y_diag=sparse(blkdiag(y_cell{:}));


    N_k=sum(Responsabilities,2);
    pi_new(:)=N_k/N;
    for i=1:modes
        responsabilities=Responsabilities(i,:);
        resp2=cellfun(@(x) x*eye(n),mat2cell(responsabilities',ones(1,N)),'UniformOutput',0);
        Gamma=sqrt(sparse(blkdiag(resp2{:})));

        Phi(i,:)=-((Gamma*theta)\(Gamma*Y(:)));

        y_lin=reshape(theta*-Phi(i,:)',n,N);
        y_lin_cell = cellfun(@(x) x',mat2cell(y_lin',ones(1,N))','UniformOutput',0);
        y_lin_diag=sparse(blkdiag(y_lin_cell{:}));
        err_diag=y_diag-y_lin_diag;
        sigma_weighted=err_diag*diag(responsabilities)*err_diag.';
        sigma_num=zeros(n);
        for j=0:N-1
            sigma_num=sigma_num+sigma_weighted(n*j+1:n*j+1+(n-1),n*j+1:n*j+1+(n-1));
        end
        Sigma(:,:,i)=sigma_num/N_k(i);
    end

end
