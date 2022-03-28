clear
close all
warning off

doplots=1; %= do plots?

genData

%% === ESTIMATION ===
modes=2^n;
PI=repmat(1/modes,1,modes);
emMaxIter=200;
maxErr=1e-18;
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

%= Estimate normal behavior
tic
Y=lambda(:,:,1);
[Phi,Responsibilities,pi_new, Sigma,loglikelihood,info] = emgm_estimate (X,Y,Phi_init,[],modes,emMaxIter,maxErr);

% Phi_init
Phi_init_orig
Phi
norm(Phi-Phi_init_orig,'fro')
toc

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
Phi_init_tilde=Phi_init_orig_tilde+1.*rand(size(Phi_init_orig_tilde));

%= Estimate selfish behavior
tic
Y=lambda_tilde(:,:,1);
[Phi_tilde,Responsibilities_tilde,~,Sigma_tilde,loglikelihood_tilde,info_tilde] = emgm_estimate (X,Y,Phi_init_tilde,[],modes,emMaxIter,maxErr);
% Phi_init_tilde
Phi_init_orig_tilde
Phi_tilde
norm(Phi_tilde-Phi_init_orig_tilde,'fro')
toc

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
