%% Clean variables
clear
close all
warning off

doplots=1; %= do plots?

genData

%% === ESTIMATION ===
modes=2^n;
emMaxIter=200;
n_small_em=10;
repeat_small_em=10;
maxErr=1e-8;
emMaxIter_pre=2;
maxErr_pre=1e-8;

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
[~,indices]=sort(sum(Phi_init_orig.^2,2));
Phi_init_orig=Phi_init_orig(flip(indices),:);

%% = Estimate normal behavior
X=theta;
Y=lambda(:,:,1);
tic
Phi_pre=zeros(modes,n^2+n,n_small_em);
for small_em=1:repeat_small_em

    % run multiple small runs (=n_msall_em= runs with =emMaxIter_pre= maximum iterations )
    for i=1:n_small_em
        Phi_init=1*rand(modes,n^2+n);
        [Phi_pre(:,:,i),Responsibilities,pi_new, Sigma_pre(:,:,:,i),loglikelihood,info] = emgm_estimate (X,Y,Phi_init,[],modes,emMaxIter_pre,maxErr_pre);
        loglikelihoods(:,i)=loglikelihood(info.step);
    end

    % save phi with best log-likelihood
    [~, indices] = sort(loglikelihoods);
    Phi_init_chosen(:,:,small_em)=Phi_pre(:,:,indices(end));
    Sigma_init_chosen(:,:,:,small_em)=Sigma_pre(:,:,:,indices(end));
    loglikelihoods_em_small(:,small_em)=loglikelihoods(indices(end));
end

[~, indices] = sort(loglikelihoods_em_small);
Phi_init_chosen_em=Phi_init_chosen(:,:,indices(end));
Sigma_init_chosen_em=Sigma_init_chosen(:,:,:,indices(end));

[Phi,Responsibilities,pi_new, Sigma,loglikelihood,info] = emgm_estimate (X,Y,Phi_init_chosen_em,Sigma_init_chosen_em,modes,emMaxIter,maxErr);

[~,indices]=sort(sum(Phi.^2,2));
Phi=Phi(flip(indices),:);
Phi_init_chosen
Phi_init_orig
Phi
norm(Phi-Phi_init_orig,'fro')
toc
% return

P_1_tilde=(T(:,:,1)*P_1).';
P_2_tilde=(T(:,:,1)*P_2).';
P_complet_tilde=(T(:,:,1)*P_complet).';
s_complet_tilde=T(:,:,1)*s_complet;
s_1_tilde=T(:,:,1)*s_1;
s_2_tilde=T(:,:,1)*s_2;
Phi_init_orig_tilde=[P_complet_tilde(:).' s_complet_tilde.';
                     P_1_tilde(:).' s_1_tilde.';
                     P_2_tilde(:).' s_2_tilde.';
                     zeros(1,n*n+n);
                    ];

[~,indices]=sort(sum(Phi_init_orig_tilde.^2,2));
Phi_init_orig_tilde=Phi_init_orig_tilde(flip(indices),:);

%% = Estimate selfish behavior
X=theta;
Y=lambda_tilde(:,:,1);
tic
Phi_tilde_pre=zeros(modes,n^2+n,n_small_em);
for small_em=1:repeat_small_em

    % run multiple small runs (=n_msall_em= runs with =emMaxIter_pre= maximum iterations )
    for i=1:n_small_em
        Phi_init_tilde=1*rand(modes,n^2+n);
        [Phi_tilde_pre(:,:,i),Responsibilities_tilde,~,Sigma_tilde(:,:,:,i),loglikelihood_tilde,info_tilde] = emgm_estimate (X,Y,Phi_init_tilde,[],modes,emMaxIter_pre,maxErr_pre);
        loglikelihoods_tilde_pre(:,i)=loglikelihood_tilde(info_tilde.step);
    end

    % save phi with best log-likelihood
    [~, indices] = sort(loglikelihoods_tilde_pre);
    Phi_init_chosen_tilde(:,:,small_em)=Phi_tilde_pre(:,:,indices(end));
    Sigma_init_chosen_tilde(:,:,:,small_em)=Sigma_tilde(:,:,:,indices(end));

    loglikelihoods_em_small_tilde(:,small_em)=loglikelihoods_tilde_pre(indices(end));
end

[~, indices] = sort(loglikelihoods_em_small_tilde);
Phi_init_chosen_em_tilde=Phi_init_chosen_tilde(:,:,indices(end));
Sigma_init_chosen_em_tilde=Sigma_init_chosen_tilde(:,:,:,indices(end));

[Phi_tilde,Responsibilities_tilde,~,Sigma_tilde,loglikelihood_tilde,info_tilde] = emgm_estimate (X,Y,Phi_init_chosen_em_tilde,Sigma_init_chosen_em_tilde,modes,emMaxIter,maxErr);

[~,indices]=sort(sum(Phi_tilde.^2,2));
Phi_tilde=Phi_tilde(flip(indices),:);
Phi_init_chosen_em_tilde
Phi_init_orig_tilde
Phi_tilde
norm(Phi_tilde-Phi_init_orig_tilde,'fro')
toc

% NOTE(accacio): only if values have zero
% index_of_zero=find(sum(theta==zeros(size(theta)))==n);
% index_of_zero=1;
% [~, z_hat_zero]=max(Responsibilities(:,index_of_zero)); %#ok
% zero_params=Phi(z_hat_zero,:);

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

%%
rgb=@(x,y,z) [x, y,z]/255;
colors={ rgb( 84, 177, 159),  ...
         rgb(217, 108,  25), ...
         rgb(211, 101, 159), ...
         rgb(128,  63, 189), ...
       };
% #54B19F #D96C19 #D3659F #803FBD

%% Plots
% if(doplots==1 && size(theta,1)==2 )
% for i=1:M
%     figure
%     for j=1:size(lambda,1)
%         sgtitle([' System ' num2str(i) ' normal behavior'],'interpreter','latex')
%         subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),j)
%         scatter3(theta(1,:),theta(2,:),lambda(j,:,i),10,'filled');
%         view(135,30)
%         title(['$\lambda_{' num2str(j) '}$'],'interpreter','latex')
%         xlabel('$\theta_1$','interpreter','latex')
%         ylabel('$\theta_2$','interpreter','latex')
%     end
% end
% end

% if(doplots==1 && size(theta,1)==2)
% for i=1:M
%     figure
%     for j=1:size(lambda,1)
%         sgtitle([' System ' num2str(i) ' cheating'],'interpreter','latex')
%         subplot(round(sqrt(2)),round(sqrt(2))+1*(round(sqrt(2))<=floor(sqrt(2))),j)
%         scatter3(theta(1,:),theta(2,:),lambda_tilde(j,:,i),10,'filled');
%         view(135,30)
%         title(['$\tilde{\lambda}_{' num2str(j) '}$'],'interpreter','latex')
%         xlabel('$\theta_1$','interpreter','latex')
%         ylabel('$\theta_2$','interpreter','latex')
%     end
% end
% end

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

figure
plot(1:info.step,loglikelihood(1:info.step))
title(['Log-likelihood (Normal Behavior)' ],'interpreter','latex')
figure
plot(1:info_tilde.step,loglikelihood_tilde(1:info_tilde.step))
title(['Log-likelihood (Cheating)' ],'interpreter','latex')
