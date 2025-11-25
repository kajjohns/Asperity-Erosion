
plot_profiles = false;

folder_name = 'Cascadia_outputs';

eval(['load ./' folder_name '/logrho.txt'])
figure
plot(logrho)


%discard burn-in samples
discard = 100;

eval(['load ./' folder_name '/M_radii.txt'])
eval(['load ./' folder_name '/M_ring_tau.txt'])
eval(['load ./' folder_name '/dhat.txt'])
eval(['load ./' folder_name '/locked_index.txt'])
eval(['load ./' folder_name '/creep_rates.txt'])



%toss out burn-in samples
M_radii(1:discard,:) = [];
M_ring_tau(1:discard,:) = [];
dhat(1:discard,:) = [];
locked_index(1:discard,:) = [];
creep_rates(1:discard,:) = [];


%construct ring_tau
clear  Ring_Taus
for k=1:size(M_radii,1)

    radii = M_radii(k,:)';
    ring_taus = M_ring_tau(k,:)';

    %forward model calculation
    i_locked = get_locked_radii_indices(ptsDZ,radii,centroids_rot,patch_stuff.dip_faces,minrad);
    
    %accumulated stress
    Ring_Taus(:,k) = get_ring_tau(ring_taus,ptsDZ,D,radii,centroids_rot,i_locked,patch_stuff.dip_faces,minrad);

end




c1 = colormap(flipud(slanCM('inferno')));
c2 = colormap(flipud(slanCM('Blues')));
mycmap = [c2(1:2:end,:); c1(1:2:end,:)];


str = mean(strikes)*pi/180; %negative sign to rotate back
R = [cos(str) -sin(str); sin(str) cos(str)];
centroids_rot = (R*centroids(:,1:2)')';
centroids_rot(:,3) = centroids(:,3);
centroids_rot(:,2) = centroids_rot(:,2)-min(centroids_rot(:,2));



nd_ll = local2llh(nd(:,1:2)',origin)';

load cmap


figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),mean(locked_index),'edgecolor','none'); 
colormap(flipud(pink))
colorbar     
set(gca,'fontsize',15)
title('Probability of Locking')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([40 50])
   


figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),mean(Ring_Taus,2),'edgecolor','none'); 
colormap(flipud(pink))
colorbar     
set(gca,'fontsize',15)
title('Mean Negative Stressing Rate')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([39 51])
   

figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),mean(creep_rates),'edgecolor','none'); 
colormap(turbo)
colorbar    
set(gca,'fontsize',15)
title('Mean Creep Rate')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([39 51])



figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),1-mean(creep_rates,1)'./(srate*1000),'edgecolor','none'); 
colormap(cmap)
colorbar
set(gca,'fontsize',15)
title('Coupling Ratio')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([40 50])
caxis([-1 1])
colormap(mycmap)

addpath Cascadia_obs/SSEs
load bartlow_SSE_contours.mat
load bartlow_tremor_contours.mat


contour(tremor_X,tremor_Y,tremor_conts,[0 60],'color','g','linewidth',2)
contour(SSE_X,SSE_Y,SSE_rate,[0 0.005],'color','m','linewidth',2)


figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),(srate*1000)-mean(creep_rates)','edgecolor','none'); 
colormap(cmap)
colorbar
set(gca,'fontsize',15)
title('Slip Deficit Rate')
view(2)
hold on
plot_coast(bbox)
xlim([-129 -120])
ylim([39 51])







%plot fit to data
Nbase =size(baselines,1);
dhat_base = dhat(:,1:Nbase);
dhat_u = dhat(:,1+Nbase:end);

sites_ll = local2llh(xysites',origin)';
vert_ll = local2llh(xysites_vert',origin)';

BaseEnds_ll = [local2llh(BaseEnds(:,1:2)',origin)' local2llh(BaseEnds(:,3:4)',origin)'];


figure; 
subplot(121)
hold on;  
axis equal
for k=1:size(BaseEnds,1)
    cline([BaseEnds_ll(k,1) BaseEnds_ll(k,3)],[BaseEnds_ll(k,2) BaseEnds_ll(k,4)],[Vbase(k)/L(k) Vbase(k)/L(k)])
end
title('Baseline elongation rates')
colormap(cmap)
colorbar
caxis([-0.1 0.1])
plot_coast(bbox)
xlim([-129 -120])
ylim([40 50])



subplot(122)
hold on;  
axis equal
v = mean(dhat,1);
for k=1:size(BaseEnds,1)
    cline([BaseEnds_ll(k,1) BaseEnds_ll(k,3)],[BaseEnds_ll(k,2) BaseEnds_ll(k,4)],[v(k)/L(k) v(k)/L(k)])
end
title('Baseline elongation rates')
colormap(cmap)
colorbar
caxis([-0.1 0.1])
plot_coast(bbox)
xlim([-127 -117])
ylim([40 50])



% 
% figure; 
% subplot(131)
% hold on;  
% axis equal
% for k=1:size(BaseEnds,1)
%     cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k)/L(k) Vbase(k)/L(k)])
% end
% title('Baseline elongation rates')
% plot_coast_xy(bbox,origin,'k')
% load cmap
% colormap(cmap)
% colorbar
% caxis(0.5*[-max(abs(Vbase./L)) max(abs(Vbase./L))])
% xlim([-200 600])
% ylim([-600 600])
% 
% set(gca,'fontsize',15)
% 
% subplot(132)
% hold on;  
% axis equal
% v = mean(dhat,1);
% for k=1:size(BaseEnds,1)
%     cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[v(k)/L(k) v(k)/L(k)])
% end
% title('Baseline elongation rates')
% plot_coast_xy(bbox,origin,'k')
% load cmap
% colormap(cmap)
% colorbar
% caxis(0.5*[-max(abs(Vbase./L)) max(abs(Vbase./L))])
% xlim([-200 600])
% ylim([-600 600])
% 
% set(gca,'fontsize',15)
% 
% 
% subplot(133)
% hold on;  
% axis equal
% v = mean(dhat,1);
% for k=1:size(BaseEnds,1)
%     cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[(Vbase(k)-v(k))/L(k) (Vbase(k)-v(k))/L(k)])
% end
% title('Residual')
% plot_coast_xy(bbox,origin,'k')
% load cmap
% colormap(cmap)
% colorbar
% caxis(0.5*[-max(abs(Vbase./L)) max(abs(Vbase./L))])
% xlim([-200 600])
% ylim([-600 600])
% 
% set(gca,'fontsize',15)
% 
% set(gcf,'renderer','Painters')
% 


%plot vertical
figure; 
subplot(121)
hold on
scatter(vert_ll(:,1),vert_ll(:,2),50,Vu,'fill')
axis equal    
title('Observed Vertical')
plot_coast_xy(bbox,origin,'k')
colormap(cmap)
colorbar
grid on
set(gca,'fontsize',15)
xlim([-200 600])
ylim([-600 600])
caxis([-5 5])

vu = mean(dhat_u,1);
subplot(122)
hold on
scatter(xysites_vert(:,1),xysites_vert(:,2),50,vu,'fill')
axis equal    
title('Predicted Vertical')
plot_coast_xy(bbox,origin,'k')
colormap(cmap)
colorbar
grid on
set(gca,'fontsize',15)
xlim([-200 600])
ylim([-600 600])
caxis([-5 5])





%plot vertical
figure; 
subplot(131)
hold on
scatter(xysites_vert(:,1),xysites_vert(:,2),50,Vu,'fill')
axis equal    
title('Observed Vertical')
plot_coast_xy(bbox,origin,'k')
colormap(cmap)
colorbar
grid on
set(gca,'fontsize',15)
xlim([-200 600])
ylim([-600 600])
caxis([-5 5])

vu = mean(dhat(:,length(Vbase)+1:end));
subplot(132)
hold on
scatter(xysites_vert(:,1),xysites_vert(:,2),50,vu,'fill')
axis equal    
title('Predicted Vertical')
plot_coast_xy(bbox,origin,'k')
colormap(cmap)
colorbar
grid on
set(gca,'fontsize',15)
xlim([-200 600])
ylim([-600 600])
caxis([-5 5])


subplot(133)
hold on
scatter(xysites_vert(:,1),xysites_vert(:,2),50,Vu-vu','fill')
axis equal    
title('Residual')
plot_coast_xy(bbox,origin,'k')
colormap(cmap)
colorbar
grid on
set(gca,'fontsize',15)
xlim([-200 600])
ylim([-600 600])
caxis([-5 5])


if plot_profiles

%plot profiles
ypos = -400:50:400;

cr = mean(creep_rates);

for k=1:length(ypos)

    ind = abs(centroids(:,2)-ypos(k))<5;
    cr_prof = cr(ind);
    xpos_prof = centroids(ind,1);
    [x,i]=sort(xpos_prof);

    %smooth profile a bit
    crp = cr_prof(i);
    crp = medfilt1(crp,6);

    ltrate = rates(ind);
    ltrate =  ltrate(i);

    figure(100+k)
    plot(x-x(1),crp)
    hold on
    plot(x-x(1),ltrate,'r')
    grid on

    xlabel('distance from trench (km)')
    ylabel('creep rate/slip rate (mm/yr)')

end

end
