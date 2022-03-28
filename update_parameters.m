function [Phi, pi_new, Sigma] = update_parameters(X, Y, OldPhi, Responsibilities)
% UPDATE_PARAMETERS -
    modes=size(Responsibilities,1);
    [n N]=size(X);

    A=kron(ones(n,1),eye(n));
    B=kron(eye(N),ones(1,n));
    C=kron(repmat(eye(n),1,N),ones(n,1));
    theta=sparse([(A*X*B).*C;kron(repmat(eye(n),1,N),ones(1,1))].');

    y_cell = cellfun(@(x) x.',mat2cell(Y.',ones(1,N)).','UniformOutput',0);
    y_diag=sparse(blkdiag(y_cell{:}));


    N_k=sum(Responsibilities,2);
    pi_new(:)=N_k/N;

    % NOTE: Sigma is spheric ($\Sigma=\lambda I$)
    % update formula is given by $\lambda_{k}^{r+1}=\frac{\zeta+\sum_{i=1}^{n} \tau_{i k}^{r}\left(y_{i}-\mu_{k}^{r+1}\right)^{\prime}\left(y_{i}-\mu_{k}^{r+1}\right)}{v+d \sum_{i=1}^{n} \tau_{i k}^{r}+d+2}$
    % as seen in http://dx.doi.org/10.1007/s11222-015-9561-x
    sigma0=1e-5;
    d=n;
    nu=d+2;
    zeta=2*sigma0^(1/d);

    for i=1:modes
        responsibilities=Responsibilities(i,:);
        resp2=cellfun(@(x) x*eye(n),mat2cell(responsibilities',ones(1,N)),'UniformOutput',0);
        Gamma=sqrt(sparse(blkdiag(resp2{:})));

        Phi(i,:)=-((Gamma*theta)\(Gamma*Y(:)));

        y_lin=reshape(theta*-Phi(i,:).',n,N);
        y_lin_cell = cellfun(@(x) x.',mat2cell(y_lin.',ones(1,N)).','UniformOutput',0);
        y_lin_diag=sparse(blkdiag(y_lin_cell{:}));
        err_diag=y_diag-y_lin_diag;

        % sigma_weighted=err_diag*diag(responsibilities)*err_diag.';
        % sigma_num=zeros(n);
        % for j=0:N-1
        %     sigma_num=sigma_num+sigma_weighted(n*j+1:n*j+1+(n-1),n*j+1:n*j+1+(n-1));
        % end
        % Sigma(:,:,i)=sigma_num/N_k(i);

        lambda_num=zeta+sum(diag((err_diag.'*err_diag).*diag(responsibilities)));
        sum_d2= sum(Responsibilities,2);
        lambda_denum=nu+d*sum_d2(i)+d+2;
        lambda=lambda_num/lambda_denum;
        Sigma(:,:,i)=lambda*eye(n);
    end

end
