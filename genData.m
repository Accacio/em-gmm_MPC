%% === GENERATION ===
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
Wt(:,1) = [25].'; %#ok
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
values={linspace(-0,7,2*(n^2+n))};
% values={linspace(-0,.1,n^2+n)}; % both active
% values={lispace(0,1,20), linspace(3,5,20),}; % 1 active
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
