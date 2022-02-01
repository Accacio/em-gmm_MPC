%% === GENERATION ===
%% Clean variables
close all
clear
warning off
%%%
%% Define systems

% Optimization settings
options = optimset('Display', 'off');
% options = [];

%= Number of systems
M=1;

Cair_mean=8;
Cwalls_mean=5;
Roaia_mean=5;
Riwia_mean=2.5;
Rowoa_mean=1;

Cair  = repmat(Cair_mean,1,M)+(-.5+rand(1,M));
Cwalls = repmat(Cwalls_mean,1,M)+(-.5+rand(1,M));
Roaia = repmat(Roaia_mean,1,M)+(-.5+rand(1,M));
Riwia = repmat(Riwia_mean,1,M)+(-.5+rand(1,M));
Rowoa = repmat(Rowoa_mean,1,M)+(-.5+rand(1,M));

%= Define systems using 3R2C
for i=1:M
    csys(:,:,1,i)=model_3R2C(Roaia(i),Riwia(i),Rowoa(i),Cwalls(i),Cair(i));
end
ni=size(csys.B(:,:,1,1),2);

clear C* R*

Te=.25;
dsys=c2d(csys,Te);
%%%
%% Set functions for MPC

% Prediction horizon
n=2;

%= Gains Q and R for $\sum_{j=1}^n \|v\|^2_{Q}+\|u\|^2_{R}$
for i=1:M
    Q(:,:,i)=eye(n*size(dsys(:,:,1,i).C,1)); % no x no
    R(:,:,i)=eye(n*size(dsys(:,:,1,i).B,2)); % nc x nc
end

paren = @(x, varargin) x(varargin{:}); %
                                       % Apply index in created elements
curly = @(x, varargin) x{varargin{:}}; %

%= Output prediction matrices, such $\vec{Y}=\mathcal{M}\vec{x}_k+\mathcal{C}\vec{U}$
Mmat_fun=@(sys,n) ...
         cell2mat( ...
    (repmat(mat2cell(sys.C,1),1,n) .* ...
     (repmat(mat2cell(sys.A,size(sys.A,2)),1,n).^num2cell(1:n)) ...
    )'  ...
                 );
Cmat_fun=@(sys, n) ...
         cell2mat( ...
    paren( ...
    (repmat(mat2cell(sys.C, 1), 1, n+1) .* ...
     (horzcat( ...
    zeros(size(sys.A)), ...
    repmat(mat2cell(sys.A, size(sys.A,2)),1,n).^num2cell(1:n)...
             )) .* ...
     repmat(mat2cell(sys.B, size(sys.B,1), size(sys.B,2)), 1, n+1)) ...
    , tril(toeplitz(1:n))+1));

%= H and f, such $\frac{1}{2}\vec{U}^TH\vec{U}+f^T\vec{U}$
H_fun=@(Cmat,Q,R) round(Cmat'*Q*Cmat+R*eye(size(Q)),10);
f_fun=@(Cmat,Mmat,Q,xt,Wt) Cmat'*Q*(Mmat*xt-Wt);


%= Prediction matrices for the systems
for i=1:M
    Mmat(:,:,i)=Mmat_fun(dsys(:,:,1,i),n);
    Cmat(:,:,i)=Cmat_fun(dsys(:,:,1,i),n);
    H(:,:,i)=H_fun(Cmat(:,:,i),Q(:,:,i),R(:,:,i));
end
clear -regexp [^f].*_fun % Delete all functions but f_fun

X0(:,1) = [21 3.2]';
X0(:,2) = [20. 6.]';
Wt(:,1) = [23]';
Wt(:,2) = [21]';

for i=1:M
    f(:,:,i)=f_fun(Cmat(:,:,i),Mmat(:,:,i),Q(:,:,1),X0(:,i),Wt(:,i));
end

% TODO(accacio): Use yalmip
% umin(1:2)=-inf;
umin(1:2)=0;
umax(1:2)=inf;
%%%
%% Get Lambdas

%= Generate n-dimensional rectangular grid for combinations
% see https://accacio.gitlab.io/blog/matlab_combinations/
% values=-4:.5:4;
% values=-10:1:10;
% values=-1:.1:1;
values=linspace(0,.001,n^3);
% values=.2:.1:1.2;
% values=max(umin(1),values);

[ v{1:n} ]=ndgrid(values);
theta(:,:) =cell2mat(cellfun(@(x) reshape(x,[],1),v,'UniformOutput',0))';

for i=1:M
    for cur_theta=1:size(theta,2)
        % QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0)
        [u(:,i) ,J(:,i),~,~,l] = quadprog(H(:,:,i), f(:,:,i), ...
                                          eye(ni*n), theta(:,cur_theta), ...
                                          [], [], ...
                                          umin(:,i)*ones(ni*n,1), ...  % Lower Bound
                                          umax(:,i)*ones(ni*n,1), ...  % Upper Bound
                                          [], options);
        % [u(:,i) ,J(:,i),~,~,l] = quadprog(H(:,:,i), f(:,:,i), ...
        %                                   eye(ni*n), min(max(theta(:,cur_theta),umin(:,i)),umax(:,i)), ...
        %                                   [], [], ...
        %                                   umin(:,i)*ones(ni*n,1), ...  % Lower Bound
        %                                   umax(:,i)*ones(ni*n,1), ...  % Upper Bound
        %                                   [], options);
        lambda(:,cur_theta,i)=l.ineqlin;
    end
end

%%%
%% Plots

if(size(theta,1)==2)
for i=1:M
    figure
    for j=1:size(lambda,1)
        sgtitle([' System ' num2str(i) ' normal behavior'],'interpreter','latex')
        subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),j)
        scatter3(theta(1,:),theta(2,:),lambda(j,:,i),10,'filled');
        view(135,30)
        title(['$\lambda_{' num2str(j) '}$'],'interpreter','latex')
        xlabel('$\theta_1$','interpreter','latex')
        ylabel('$\theta_2$','interpreter','latex')
    end
end
end

% Get only lambdas different from zero
% for i=1:M
% suberror=1e-4;
% colIdx=find(sum(lambda(:,:,i)>suberror)==size(lambda(:,:,i),1));
% lambda_est=paren(lambda(:,:,i),':',colIdx);
% theta_est=paren(theta,':',colIdx);

% figure
% for j=1:size(lambda,1)
%     sgtitle([' System ' num2str(i) ' normal behavior no $\lambda=0$'],'interpreter','latex')
%     subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),j)
%     scatter3(theta_est(1,:),theta_est(2,:),lambda_est(j,:));
% view(135,30)
%     title(['$\lambda_{' num2str(j) '}$ '],'interpreter','latex')
%     xlabel('$\theta_1$','interpreter','latex')
%     ylabel('$\theta_2$','interpreter','latex')
% end
% end


%%%
%% Simulate
% u=ones(1,n);
% [y,t,x]=lsim(dsys(:,:,1,2),[u(:,1) 0],0:Te:Te); % Update systems
%%%
%= Cheating influence on lambda
% let's suppose $\tilde{\vec{\lambda}}=T\vec{\lambda}$

T = (10*rand([size(lambda,1) size(lambda,1) M]));
% T = 20*eye(n)
% T = diag(fix(20*rand(1,n)));
T = (T+T')/2;
% T=T+2*eye(size(lambda,1))
for i=1:M
    lambda_tilde(:,:,i) = T*lambda(:,:,i);
end

if(size(theta,1)==2)
for i=1:M
    figure
    for j=1:size(lambda,1)
        sgtitle([' System ' num2str(i) ' cheating'],'interpreter','latex')
        subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),j)
        scatter3(theta(1,:),theta(2,:),lambda_tilde(j,:,i),10,'filled');
        view(135,30)
        title(['$\tilde{\lambda}_{' num2str(j) '}$'],'interpreter','latex')
        xlabel('$\theta_1$','interpreter','latex')
        ylabel('$\theta_2$','interpreter','latex')
    end
end
end

%% === ESTIMATION ===
PI=pi;
emMaxIter=500;
modes=n^2;
maxErr=1e-8;
x=theta;

%= Estimate normal
for component=1:n;
    y=lambda(component,:);
    [C(:,:,:,component),d(:,:,:,component),responsabilities(:,:,component),~, ~] = emgm_estimate (x,y,[],[],modes,emMaxIter,maxErr);
    C_estimated(:,:,component) = -reshape(C(:,:,:,component),size(C,1),size(C,3));
    d_estimated(:,:,component) = -reshape(d(:,:,:,component),size(d,1),size(d,3));
end
return
for component=1:n
    index_of_zero=find(sum(theta==zeros(size(theta)))==n);
    % index_of_zero=1;
    [~,z_hat_zero]=max(responsabilities(:,index_of_zero,component));
    H_est(component,:,i)=C_estimated(:,z_hat_zero,component)';
    f_est(component,:,i)=d_estimated(:,z_hat_zero,component);
end
H_est(:,:,i)
H(:,:,1)
f_est(:,:,i)
f(:,:,1)

%= Estimate cheating
for component=1:n;
    y=lambda_tilde(component,:);

    % Cinit= reshape(repmat(paren(T*H_est,1,':')',1,modes),2,1,modes)+rand(2,1,modes);
    % dinit = reshape(repmat(paren(T*f_est,1,':')',1,modes),1,1,modes)+rand(1,1,modes);
    % Cinit= reshape(repmat(paren(T*H_est,1,':')',1,modes),2,1,modes);
    % dinit = reshape(repmat(paren(T*f_est,1,':')',1,modes),1,1,modes);

    % Cinit=C(:,:,:,component);
    % dinit=d(:,:,:,component);
    % [C_cheat(:,:,:,component),d_cheat(:,:,:,component),responsabilities_cheat(:,:,component),~, ~] = emgm_estimate (x,y,Cinit,dinit,modes,emMaxIter,maxErr);
    [C_cheat(:,:,:,component),d_cheat(:,:,:,component),responsabilities_cheat(:,:,component),~, ~] = emgm_estimate (x,y,[],[],modes,emMaxIter,maxErr);
    C_cheat_estimated(:,:,component) = -reshape(C_cheat(:,:,:,component),size(C_cheat,1),size(C_cheat,3));
    d_cheat_estimated(:,:,component) = -reshape(d_cheat(:,:,:,component),size(d_cheat,1),size(C_cheat,3));
end
% TODO(accacio): should not necessarily be zero, can be any theta < -P\s
for component=1:n
    index_of_zero=find(sum(theta==zeros(size(theta)))==n);
    % index_of_zero=1;
    [~,z_hat_zero]=max(responsabilities(:,index_of_zero,component));
    H_cheat_est(component,:,i)=C_cheat_estimated(:,z_hat_zero,component)';
    f_cheat_est(component,:,i)=d_cheat_estimated(:,z_hat_zero,component);
end

H_cheat_est(:,:,i)
T*H(:,:,1)
f_cheat_est(:,:,i)
T*f(:,:,1)


% C_estimated
% C_cheat_estimated

invT_est=H*inv(H_cheat_est)
invT=inv(T)

%%
rgb=@(x,y,z) [x, y,z]/255;
colors={ rgb(84, 177, 159),
         rgb(217, 108, 25),
         rgb(211, 101, 159),
         rgb(128, 63, 189),
       };

% Plot normal Behavior

if(size(theta,1)==2)
figure
for component=1:n
    y=lambda(component,:);
    sgtitle('EM-GM Using MPC Data (Normal Behavior)','interpreter','latex')
    subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),component)
    plot_responsibles(x, y, responsabilities(:,:,component), C, d, colors);
    view(135,30)
    title(['$\lambda_{' num2str(component) '}$' ],'interpreter','latex')
    xlabel('$\theta_1$','interpreter','latex')
    ylabel('$\theta_2$','interpreter','latex')
end
end

% Plot cheating
if(size(theta,1)==2)
figure
for component=1:n
    y=lambda_tilde(component,:);
    sgtitle('EM-GM Using MPC Data (Cheating)','interpreter','latex');
    subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),component)
    plot_responsibles(x, y, responsabilities_cheat(:,:,component), C_cheat, d_cheat, colors);
    view(135,30)
    title([' $\lambda_{' num2str(component) '}$' ],'interpreter','latex')
    xlabel('$\theta_1$','interpreter','latex')
    ylabel('$\theta_2$','interpreter','latex')
end
end
