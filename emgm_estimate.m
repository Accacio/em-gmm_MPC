function [C,d,responsabilities,pi,Sigma] = emgm_estimate(x,y,m,n,modes,emMaxIter,maxErr)
% EMGM_ESTIMATE -

    pi=repmat(1/modes,1,modes);
    OldclusterSize=zeros(1,modes);

    if(isempty(m)||isempty(n))
        indexpts=randi([1 size(x,2)-1],1,modes);
        dx=max(x(:,2:end)-x(:,1:end-1),[],2);
        dx=num2cell(dx,2);

        y_square=reshape(y,sqrt(size(y,2)),sqrt(size(y,2)));
        [V{1:2}]=gradient(y_square,dx{:});
        grad=cell2mat(cellfun(@(x) reshape(x,[],1),V,'UniformOutput',0))';
        C=grad(:,indexpts);
        x_indexed=x(:,indexpts);
        y_indexed=y(:,indexpts);
        d=y_indexed-sum(C.*x_indexed);
        C=reshape(C,2,1,modes);
        d=reshape(d,1,1,modes);
        % C=5*rand(size(x,1),1,modes);
        % d=5*rand(1,1,modes);
    else
        C=m;
        d=n;
    end


    eps=1000;
    Sigma(:,:,1:modes)=repmat(eps*eye(1),1,1,modes);

    for emInd=1:emMaxIter

        responsabilities=calculate_responsabilities(x,y,C,d,Sigma,pi);

        [C, d, pi] = update_parameters(x, y, responsabilities);
        [~,z_hat]=max(responsabilities,[],1);
        for i=1:modes
            z_i=find(z_hat==i);
            clusterSize(i)=size(z_i,2);
            if OldclusterSize(i)==clusterSize(i)
                Sigma(:,:,i)=Sigma(:,:,i)*.9;
            else
                Sigma(:,:,i)=Sigma(:,:,i)*1.;
            end
        end
        if sum(Sigma < maxErr)==modes
            break;
        end
        OldclusterSize=clusterSize;
    end

end
