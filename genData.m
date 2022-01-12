%% GENERATION ===
%% Clean variables
close all
clear
%%%
%% Define systems

% Optimization settings
options = optimset('Display', 'off');
% options = [];

% Number of systems
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

% Define systems using 3R2C
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

% Gains Q and R for $\sum_{j=1}^n \|v\|^2_{Q}+\|u\|^2_{R}$
for i=1:M
    Q(:,:,i)=eye(n*size(dsys(:,:,1,i).C,1)); % no x no
    R(:,:,i)=eye(n*size(dsys(:,:,1,i).B,2)); % nc x nc
end

paren = @(x, varargin) x(varargin{:}); %
                                       % Apply index in created elements
curly = @(x, varargin) x{varargin{:}}; %

% Output prediction matrices, such $\vec{Y}=\mathcal{M}\vec{x}_k+\mathcal{C}\vec{U}$
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

% H and f, such $\frac{1}{2}\vec{U}^TH\vec{U}+f^T\vec{U}$
H_fun=@(Cmat,Q,R) round(Cmat'*Q*Cmat+R*eye(size(Q)),10);
f_fun=@(Cmat,Mmat,Q,xt,Wt) Cmat'*Q*(Mmat*xt-Wt);


% Prediction matrices for the systems
for i=1:M
    Mmat(:,:,i)=Mmat_fun(dsys(:,:,1,i),n);
    Cmat(:,:,i)=Cmat_fun(dsys(:,:,1,i),n);
    H(:,:,i)=H_fun(Cmat(:,:,i),Q(:,:,1),R(:,:,1));
end
clear -regexp [^f].*_fun % Delete all functions but f_fun

X0(:,1) = [19 3.2]';
X0(:,2) = [20. 6.]';
Wt(:,1) = [21]';
Wt(:,2) = [21]';

for i=1:M
    f(:,:,i)=f_fun(Cmat(:,:,i),Mmat(:,:,i),Q(:,:,1),X0(:,i),Wt(:,i));
end

% TODO(accacio): Use yalmip
umin(1:2)=-inf;
umax(1:2)=inf;
%%%
%% Get Lambdas

% Generate n-dimensional rectangular grid for combinations
% see https://accacio.gitlab.io/blog/matlab_combinations/
% values=0:.05:4;
values=-10:1:10;
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
        lambda(:,cur_theta,i)=l.ineqlin;
    end
end

%%%
%% Plot

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
%% Cheating influence on lambda
% let's suppose $\tilde{\vec{\lambda}}=T\vec{\lambda}$

T = ceil(rand([size(lambda,1) size(lambda,1) M])*10);
for i=1:M
    lambda_tilde(:,:,i) = T(:,:,i)*lambda(:,:,i);
end

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


%% ESTIMATION ===
PI=pi;
emMaxIter=500;
modes=n^2;
maxErr=1e-4;
x=theta;

% estimate normal
for component=1:n;
    y=lambda(component,:);
    [C(:,:,:,component),d(:,:,:,component),responsabilities(:,:,component),~, ~] = emgm_estimate(x,y,modes,emMaxIter,eps);
    C_estimated(:,:,component) = sort(reshape(C(:,:,:,component),size(C,1),size(C,3)));
end

% estimate cheating
for component=1:n;
    y=lambda_tilde(component,:);
    [C_cheat(:,:,:,component),d_cheat(:,:,:,component),responsabilities_cheat(:,:,component),~, ~] = emgm_estimate(x,y,modes,emMaxIter,eps);
    C_cheat_estimated(:,:,component) = sort(reshape(C_cheat(:,:,:,component),size(C_cheat,1),size(C_cheat,3)));
end

C_estimated
C_cheat_estimated

%%
rgb=@(x,y,z) [x, y,z]/255;
colors={ rgb(84, 177, 159),
         rgb(217, 108, 25),
         rgb(211, 101, 159),
         rgb(128, 63, 189),
       };

% Plot normal Behavior
figure
for component=1:n
    y=lambda(component,:);
    sgtitle(['EM-GM Using MPC Data (Normal Behavior)'],'interpreter','latex')
    subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),component)
    plot_responsibles(x, y, responsabilities(:,:,component), C, d, colors);
    view(135,30)
    title([' $\lambda_{' num2str(component) '}$' ],'interpreter','latex')
    xlabel('$\theta_1$','interpreter','latex')
    ylabel('$\theta_2$','interpreter','latex')
end

% Plot cheating
figure
for component=1:n
    y=lambda_tilde(component,:);
    sgtitle(['EM-GM Using MPC Data (Cheating)'],'interpreter','latex')
    subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),component)
    plot_responsibles(x, y, responsabilities_cheat(:,:,component), C_cheat, d_cheat, colors);
    view(135,30)
    title([' $\lambda_{' num2str(component) '}$' ],'interpreter','latex')
    xlabel('$\theta_1$','interpreter','latex')
    ylabel('$\theta_2$','interpreter','latex')
end
