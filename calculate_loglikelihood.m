function log_likelihood = calculate_loglikelihood(x,y,Phi,Sigma,Pi)
% CALCULATE_LOGLIKELIHOOD -
% Calculate $\sum_{n=1}^{N}\ln\{\sum_{k=1}^{K}\pi_k\mathcal{N}(x_n| \mu_k,\Sigma_k)\}$

    modes=size(Phi,1);
    [n, N]=size(x);

    respNum = zeros(size(x,2),modes);


    A=kron(ones(n,1),eye(n));
    B=kron(eye(N),ones(1,n));
    C=kron(repmat(eye(n),1,N),ones(n,1));
    x_lin=sparse([(A*x*B).*C;kron(repmat(eye(n),1,N),ones(1,1))]);
    y_cell = cellfun(@(x) x.',mat2cell(y.',ones(1,N)).','UniformOutput',0);
    y_diag=sparse(blkdiag(y_cell{:}));

    %= Calculate $\sum_{k=1}^{K}\pi_k\mathcal{N}(x_n| \mu_k,\Sigma_k)$
    for i=1:modes
        y_lin=reshape(x_lin.'*-Phi(i,:).',n,N);
        y_lin_cell = cellfun(@(x) x.',mat2cell(y_lin.',ones(1,N)).','UniformOutput',0);
        y_lin_diag=sparse(blkdiag(y_lin_cell{:}));
        err_diag=y_diag-y_lin_diag;
        SigmaInvMat=kron(sparse(eye(N)),sparse(Sigma(:,:,i)\eye(n)));

        respNum(:,i)=sparse(Pi(1,i)*(1/(2*pi)^(n/2))*(1/sqrt(det(Sigma(:,:,i))))*exp(-1/2*max(err_diag.'*SigmaInvMat*err_diag)).');
    end
    sum_pi_norm = sum(respNum,2);
    log_likelihood =sum(log(sum_pi_norm));
end
