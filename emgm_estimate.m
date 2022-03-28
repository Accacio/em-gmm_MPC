function [Phi,Responsibilities,pi_new,Sigma,loglikelihood,info] = emgm_estimate(X,Y,phi_init,sigma_init,modes,emMaxIter,maxErr)
% EMGM_NESTIMATE - ESTIMATE N DIMENSIONAL

    Pi=repmat(1/modes,1,modes);
    OldclusterSize=zeros(1,modes);
    [n, ~]=size(X);

    if(isempty(phi_init))
        % TODO(accacio): do some multidimensional magic
        % indexpts=randi([1 size(x,2)-1],1,modes);
        % dx=max(x(:,2:end)-x(:,1:end-1),[],2);
        % dx=num2cell(dx,2);

        % y_square=reshape(y,sqrt(size(y,2)),sqrt(size(y,2)));
        % [V{1:2}]=gradient(y_square,dx{:});
        % grad=cell2mat(cellfun(@(x) reshape(x,[],1),V,'UniformOutput',0))';
        % C=grad(:,indexpts);
        % x_indexed=x(:,indexpts);
        % y_indexed=y(:,indexpts);
        % d=y_indexed-sum(C.*x_indexed);
        % C=reshape(C,2,1,modes);
        % d=reshape(d,1,1,modes);
        % C=5*rand(size(x,1),1,modes);
        % d=5*rand(1,1,modes);
        Phi=20*rand(modes,n^2+n);

    else
        Phi=phi_init;
    end

    % OldPhi=zeros(size(Phi));
    OldPhi=Phi;

    if(isempty(sigma_init))
        eps=10;
        % eps=10000;
        Sigma(:,:,1:modes)=repmat(eps*eye(n),1,1,modes);
    else
        % Sigma=reshape(kron(sigma_init,eye(n)),n,n,modes)
        Sigma=sigma_init;
    end

    loglikelihood = zeros(1,emMaxIter);
    for emInd=1:emMaxIter

        Responsibilities=calculate_responsibilities(X,Y,Phi,Sigma,Pi);
        % newR=Responsibilities./max(Responsibilities)
        % newR(newR>=0.5)=1;
        % newR(newR<0.5)=0;
        % Responsibilities=newR;
        % Responsibilities(Responsibilities<0.5)=0;
        [Phi, pi_new, Sigma] = update_parameters(X, Y, Phi, Responsibilities);
        Pi=pi_new;
        loglikelihood(emInd) = calculate_loglikelihood(X,Y,Phi,Sigma,pi_new);

        % [Phi, pi_new, ~] = update_parameters(X, Y, Phi, Responsibilities);
        % Phi
        % Responsibilities
        % if (sum(sum(isnan(Phi))) | sum(sum(sum(isnan(Sigma)))) | sum(sum(isnan(pi_new))))
        %     disp('t')
        % end
        
        % [~,z_hat]=max(Responsibilities,[],1);
        % for i=1:modes
        %     z_i=find(z_hat==i);
        %     clusterSize(i)=size(z_i,2);
        %     if OldclusterSize(i)==clusterSize(i)
        %         % Sigma(:,:,i)=Sigma(:,:,i)*.9;
        %         Sigma(:,:,i)=Sigma(:,:,i)*.01;
        %     else
        %         Sigma(:,:,i)=Sigma(:,:,i)*1.;
        %     end
        % end

        % sigma_converged=1;
        % for i=1:modes
        %     sigma_converged= sigma_converged & Sigma(1,1,1)<maxErr;
        % end
        % if sigma_converged
        %     disp(['Sigma converged after ' num2str(emInd) ' iter'])
        %     info.step = emInd;
        %     return;
        % end

        if norm(OldPhi-Phi,'fro')<maxErr
        % if emInd>1 && abs(loglikelihood(emInd)-loglikelihood(emInd-1))<maxErr
            disp(['Phi Converged after ' num2str(emInd) ' iter'])
            info.step = emInd;
            return;
        end

        if emInd>1 && abs(loglikelihood(emInd)-loglikelihood(emInd-1))<maxErr
            disp(['loglikelihood Converged after ' num2str(emInd) ' iter'])
            info.step = emInd;
            return;
        end
        % OldclusterSize=clusterSize;
        OldPhi=Phi;
    end
    info.step = emInd;
    % disp(['Max iterations reached'])
end
