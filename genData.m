%% === GENERATION ===
%% Clean variables
tic
close all
clear
warning off

doplots=1; %= do plots?
%%%
%% Define systems

% Optimization settings
options = optimset('Display', 'off');
% options = [];
% rand('seed',2);

%= Number of systems
M=1;

Cair_mean=8;
Cwalls_mean=5;
Roaia_mean=5;
Riwia_mean=2.5;
Rowoa_mean=1.1;

Cair  = repmat(Cair_mean,1,M)+(-.5+rand(1,M));
Cwalls = repmat(Cwalls_mean,1,M)+(-.5+rand(1,M));
Roaia = repmat(Roaia_mean,1,M)+(-.5+rand(1,M));
Riwia = repmat(Riwia_mean,1,M)+(-.5+rand(1,M));
Rowoa = repmat(Rowoa_mean,1,M)+(-.5+rand(1,M)); %#ok

%= Define systems using 3R2C
for i=M:-1:1 % make it backward to "preallocate"
    csys(:,:,1,i)=model_3R2C(Roaia(i),Riwia(i),Rowoa(i),Cwalls(i),Cair(i));
end
ni=size(csys.B(:,:,1,1),2);

clear C* R*

Te=.25;
dsys=c2d(csys,Te);
%%%
%% Set functions for MPC


n=2; %= Prediction horizon
% n=2; %~ 0.06s
% n=3; %~ 0.09s
% n=4; %~ 0.18s
% n=5; %~ 0.4s
% n=6; %~ 1.3s
% n=7; %~ 6.2s
% n=8; %~ 39.8s

%= Gains Q and R for $\sum_{j=1}^n \|v\|^2_{Q}+\|u\|^2_{R}$
for i=M:-1:1 % make it backward to "preallocate"
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
H_fun=@(Cmat,Q,R) round(Cmat.'*Q*Cmat+R*eye(size(Q)),10);
f_fun=@(Cmat,Mmat,Q,xt,Wt) Cmat.'*Q*(Mmat*xt-Wt);


%= Prediction matrices for the systems
for i=M:-1:1
    Mmat(:,:,i)=Mmat_fun(dsys(:,:,1,i),n);
    Cmat(:,:,i)=Cmat_fun(dsys(:,:,1,i),n);
    H(:,:,i)=H_fun(Cmat(:,:,i),Q(:,:,i),R(:,:,i));
end
clear -regexp [^f].*_fun % Delete all functions but f_fun

%= Initial state and reference
X0(:,1) = [21 3.2].';
X0(:,2) = [20. 6.].';
Wt(:,1) = [20].'; %#ok
Wt(:,2) = [21].'; %#ok

for i=M:-1:1
    f(:,:,i)=f_fun(Cmat(:,:,i),Mmat(:,:,i),Q(:,:,1),X0(:,i),Wt(:,i));
end

% TODO(accacio): Use yalmip
umin(1:2)=-inf;
% umin(1:2)=0;
umax(1:2)=inf;

Gamma=eye(ni);
Gamma_bar=kron(eye(n),Gamma);

%%%
%% Get Lambdas

%= Generate n-dimensional rectangular grid for combinations
% see https://accacio.gitlab.io/blog/matlab_combinations/
% values={linspace(0,.1,2)};
values={linspace(0,.1,10)};
% values={linspace(0,1,20)}; % both active
% values={linspace(0,1,20), linspace(3,5,20),}; % 1 active
% values={linspace(4,10,20),linspace(0,1,20)}; % 2 active
% values={linspace(0,3,20), linspace(10,30,20),};

% [ v{1:n} ]=ndgrid(values,2*values);
[ v{1:n} ]=ndgrid(values{:});
theta(:,:) =cell2mat(cellfun(@(x) reshape(x,[],1),v,'UniformOutput',0))';

lambda=zeros(n,size(theta,2),M);
u=zeros(n,M);
J=zeros(n,M);
for i=1:M
    for cur_theta=1:size(theta,2)
        % QUADPROG(H,f,A,b,Aeq,beq,LB,UB,X0)
        [u(:,i) ,J(:,i),~,~,l] = quadprog(H(:,:,i), f(:,:,i), ...
                                          Gamma_bar, theta(:,cur_theta), ...
                                          [], [], ...
                                          umin(:,i)*ones(ni*n,1), ...  % Lower Bound
                                          umax(:,i)*ones(ni*n,1), ...  % Upper Bound
                                          [], options);
        lambda(:,cur_theta,i)=l.ineqlin;
    end
end
%%%
%% Simulate
% u=ones(1,n);
% [y,t,x]=lsim(dsys(:,:,1,2),[u(:,1) 0],0:Te:Te); % Update systems
%%%
%= Selfish behavior influence on lambda
% let's suppose $\tilde{\vec{\lambda}}=T\vec{\lambda}$

T = (10*rand([size(lambda,1) size(lambda,1) M]));
% T = repmat(20*(eye(n)),1,1,2)
% T = diag(fix(20*rand(1,n)));
% T=T+2*eye(size(lambda,1))
for i=M:-1:1
    T(:,:,i) = (T(:,:,i)+T(:,:,i).')/2;
    lambda_tilde(:,:,i) = T(:,:,i)*lambda(:,:,i);
end


%% === ESTIMATION ===
modes=2^n;
PI=repmat(1/modes,1,modes);
emMaxIter=200;
maxErr=1e-8;
X=theta;

%= Initialize estimation
P_1=[1/(Gamma_bar(1,:)*inv(H(:,:,1))*Gamma_bar(1,:).') 0; 0 0].';
P_2=[0 0; 0 1/(Gamma_bar(2,:)*inv(H(:,:,1))*Gamma_bar(2,:).')].';
P_complet=inv(Gamma_bar*inv(H(:,:,1))*Gamma_bar).';
s_complet=P_complet*Gamma_bar*inv(H(:,:,1))*f;
s_1=[inv(Gamma_bar(1,:)*inv(H(:,:,1))*Gamma_bar(1,:).')*Gamma_bar(1,:)*inv(H(:,:,1))*f; 0];
s_2=[0 ;inv(Gamma_bar(2,:)*inv(H(:,:,1))*Gamma_bar(2,:).')*Gamma_bar(2,:)*inv(H(:,:,1))*f];
Phi_init_orig=[P_complet(:).' s_complet.';
          P_1(:).' s_1.';
          P_2(:).' s_2.';
          zeros(1,n*n+n);
         ];
Phi_init=Phi_init_orig+1.*rand(size(Phi_init_orig));
% Phi_init=2*rand(modes,n^2+n);

%= Estimate normal behavior
Y=lambda(:,:,1);
% [Phi,Responsibilities,~, ~] = emgm_Nestimate (X,Y,[],modes,emMaxIter,maxErr);
[Phi,Responsibilities,pi_new, Sigma] = emgm_estimate (X,Y,Phi_init,modes,emMaxIter,maxErr);
Phi_init
Phi_init_orig
Phi

P_1_tilde=(T(:,:,1)*P_1).';
P_2_tilde=(T(:,:,1)*P_2).';
P_complet_tilde=(T(:,:,1)*P_complet).';
s_complet_tilde=T(:,:,1)*s_complet;
s_1_tilde=T(:,:,1)*[inv(Gamma_bar(1,:)*inv(H(:,:,1))*Gamma_bar(1,:).')*Gamma_bar(1,:)*inv(H(:,:,1))*f; 0];
s_2_tilde=T(:,:,1)*[0 ;inv(Gamma_bar(2,:)*inv(H(:,:,1))*Gamma_bar(2,:).')*Gamma_bar(2,:)*inv(H(:,:,1))*f];
Phi_init_orig_tilde=[P_complet_tilde(:).' s_complet_tilde.';
          P_1_tilde(:).' s_1_tilde.';
          P_2_tilde(:).' s_2_tilde.';
          zeros(1,n*n+n);
         ];
Phi_init_tilde=Phi_init_orig_tilde+1.*rand(size(Phi_init_orig_tilde));
% Phi_init_tilde=2*rand(modes,n^2+n);

%= Estimate selfish behavior
Y=lambda_tilde(:,:,1);
[Phi_tilde,Responsibilities_tilde,~, Sigma_tilde] = emgm_estimate (X,Y,Phi_init_tilde,modes,emMaxIter,maxErr);
Phi_init_tilde
Phi_init_orig_tilde
Phi_tilde

%
% % NOTE(accacio): only if values have zero
% index_of_zero=find(sum(theta==zeros(size(theta)))==n);
% index_of_zero=1;
% [~, z_hat_zero]=max(Responsibilities(:,index_of_zero)); %#ok
% zero_params=Phi(z_hat_zero,:);

%
% NOTE(accacio): only if values have zero
% index_of_zero=find(sum(theta==zeros(size(theta)))==n);
% index_of_zero=1;

% [~, z_hat_zero_tilde]=max(Responsibilities_tilde(:,index_of_zero));
% zero_params_tilde=Phi_tilde(z_hat_zero_tilde,:);
% H_est_tilde=reshape(zero_params_tilde(1:n^2),n,n).';
% display(H_est_tilde);
% H_tilde=T(:,:,1)*H(:,:,1);
% display(H_tilde);
% f_est=zero_params_tilde((n^2+1):end).';
% display(f_est);
% display(T(:,:,1)*f(:,:,1));

% invT_est=H(:,:,1)/(H_est_tilde);
% invT=inv(T(:,:,1));
% display(invT_est);
% display(invT);
toc

% syms('P%d%d',[2 2])
% invP=subexpr(inv(P))
% Delta=det(P)
% active_none=sym(zeros(2));
% active_1=sym(zeros(2)); active=[1]; inactive=setdiff(1:2,active);
% active_1(active,active)=subs(simplify(eval(adjoint(invP(active,active)).'./det(P(inactive,inactive))*det(P))),Delta,'Delta')
% active_2=sym(zeros(2)); active=[2]; inactive=setdiff(1:2,active);
% active_2(active,active)=subs(simplify(eval(adjoint(invP(active,active)).'./det(P(inactive,inactive))*det(P))),Delta,'Delta')

% P=inv(Gamma_bar*inv(H(:,:,1))*Gamma_bar.')
% invP=inv(P)
% active_1=zeros(n); active=[1]; inactive=setdiff(1:n,active);
% adjoint(invP(active,active)).'./det(P(inactive,inactive))*det(P)
% H_est
% active_1(active,active)=adjoint(invP(active,active)).'./det(P(inactive,inactive))*det(P)
% active_2=zeros(n); active=[2]; inactive=setdiff(1:n,active);
% active_2(active,active)=adjoint(invP(active,active)).'./det(P(inactive,inactive))*det(P)
% -det(P)/P(2,2,1)

%%
rgb=@(x,y,z) [x, y,z]/255;
colors={ rgb( 84, 177, 159),  ...
         rgb(217, 108,  25), ...
         rgb(211, 101, 159), ...
         rgb(128,  63, 189), ...
       };
% #54B19F #D96C19 #D3659F #803FBD

%% Plots

if(doplots==1 && size(theta,1)==2 )
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


if(doplots==1 && size(theta,1)==2)
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


%= Plot normal Behavior
if(doplots ==1 && size(theta,1)==2)
figure
for component=1:n
    y=lambda(component,:);
    sgtitle('EM-GM Using MPC Data (Normal Behavior)','interpreter','latex')
    subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),component)
    plot_responsibles(X, y, Responsibilities, colors);
    view(135,30)
    title(['$\lambda_{' num2str(component) '}$' ],'interpreter','latex')
    xlabel('$\theta_1$','interpreter','latex')
    ylabel('$\theta_2$','interpreter','latex')
end
end

%= Plot cheating
if(doplots ==1 && size(theta,1)==2)
figure
for component=1:n
    y=lambda_tilde(component,:);
    sgtitle('EM-GM Using MPC Data (Cheating)','interpreter','latex');
    subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),component)
    plot_responsibles(X, y, Responsibilities_tilde, colors);
    view(135,30)
    title([' $\lambda_{' num2str(component) '}$' ],'interpreter','latex')
    xlabel('$\theta_1$','interpreter','latex')
    ylabel('$\theta_2$','interpreter','latex')
end
end
