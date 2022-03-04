function plot_responsibles(x, y, responsabilities, colors)
% PLOT_RESPONSIBLES -
  [~,z_hat]=max(responsabilities,[],1);
  M=size(responsabilities,1);
  hold on

  % for i=1:M
  %     x_m=C(:,:,i).*x;
  %     [X1 ,X2 ] = meshgrid(x(1,:),x(2,:));
  %     [Xm1 ,Xm2 ] = meshgrid(x_m(1,:),x_m(2,:));
  %     Z=Xm1+Xm2+d(:,:,i);

  %     CO(:,:,1) = colors{i}(1)*ones(size(x,2)); % red
  %     CO(:,:,2) = colors{i}(2)*ones(size(x,2)); % green
  %     CO(:,:,3) = colors{i}(3)*ones(size(x,2)); % blue
  %     surf(X1,X2,Z,CO,'LineStyle',':')
  % end



  for i=1:M
    z_i=find(z_hat==i);
    scatter3(x(1,z_i),x(2,z_i),y(z_i),10,colors{i},'filled')
  end
  hold off
end
