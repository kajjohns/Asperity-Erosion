

%folder name with output files
folder_name = 'Hikurangi_outputs';

eval(['load ./' folder_name '/logrho.txt'])
figure
plot(logrho)


%discard burn-in samples
discard = 1;

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





%prepare slow slip zones
addpath ./Hikurangi_slow_slip_zones
zone1 =csvread('./Hikurangi_slow_slip_zones/zone1.csv',1,0);
zone2 =csvread('./Hikurangi_slow_slip_zones/zone2.csv',1,0);
zone3 =csvread('./Hikurangi_slow_slip_zones/zone3.csv',1,0);
zone4 =csvread('./Hikurangi_slow_slip_zones/zone4.csv',1,0);
zone1(:,2) = -zone1(:,2);
zone2(:,2) = -zone2(:,2);
zone3(:,2) = -zone3(:,2);
zone4(:,2) = -zone4(:,2);





%construct ring_tau

for k=1:size(locked_index,1)

    radii = M_radii(k,:)';
    ring_taus = M_ring_tau(k,:)';

    %forward model calculation
    %i_locked = get_locked_radii_indices(ptsDZ,radii,centroids_rot,patch_stuff.dip_faces,minrad);
    
    %accumulated stress
    %Ring_Taus(:,k) = get_ring_tau(ring_taus,ptsDZ,D,radii,centroids_rot,i_locked,patch_stuff.dip_faces,minrad);

    [i_locked, Ring_Taus(:,k) ] = get_locked_and_ring_tau(...
    ptsDZ, radii, ring_taus, D, centroids_rot, minrad, ...
    cos_t, sin_t, tan_dips, KDTree);


end


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
ylim([-43 -37])
xlim([172.5 180])
   

plot([zone1(:,1);zone1(1,1)],[zone1(:,2);zone1(1,2)],'b','linewidth',3)
plot([zone2(:,1);zone2(1,1)],[zone2(:,2);zone2(1,2)],'b','linewidth',3)
plot([zone2(:,1);zone2(1,1)],[zone2(:,2);zone2(1,2)],'b','linewidth',3)
plot([zone3(:,1);zone3(1,1)],[zone3(:,2);zone3(1,2)],'b','linewidth',3)
plot([zone4(:,1);zone4(1,1)],[zone4(:,2);zone4(1,2)],'b','linewidth',3)


figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),mean(creep_rates),'edgecolor','none'); 
colormap(turbo)
colorbar    
set(gca,'fontsize',15)
title('Mean Creep Rate')
view(2)
hold on
plot_coast(bbox)
ylim([-43 -37])
xlim([172.5 180])


plot([zone1(:,1);zone1(1,1)],[zone1(:,2);zone1(1,2)],'b','linewidth',3)
plot([zone2(:,1);zone2(1,1)],[zone2(:,2);zone2(1,2)],'b','linewidth',3)
plot([zone2(:,1);zone2(1,1)],[zone2(:,2);zone2(1,2)],'b','linewidth',3)
plot([zone3(:,1);zone3(1,1)],[zone3(:,2);zone3(1,2)],'b','linewidth',3)
plot([zone4(:,1);zone4(1,1)],[zone4(:,2);zone4(1,2)],'b','linewidth',3)

coupling_ratio = 1-mean(creep_rates)'./(srate*1000);
figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),coupling_ratio,'edgecolor','none'); 
colormap(cmap)
colorbar
set(gca,'fontsize',15)
title('Coupling Ratio')
view(2)
hold on
plot_coast(bbox)
ylim([-43 -37])
xlim([172.5 180])
caxis([-1 1])

plot([zone1(:,1);zone1(1,1)],[zone1(:,2);zone1(1,2)],'b','linewidth',3)
plot([zone2(:,1);zone2(1,1)],[zone2(:,2);zone2(1,2)],'b','linewidth',3)
plot([zone2(:,1);zone2(1,1)],[zone2(:,2);zone2(1,2)],'b','linewidth',3)
plot([zone3(:,1);zone3(1,1)],[zone3(:,2);zone3(1,2)],'b','linewidth',3)
plot([zone4(:,1);zone4(1,1)],[zone4(:,2);zone4(1,2)],'b','linewidth',3)

slip_deficit_rate = (srate*1000)-mean(creep_rates)';
figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),slip_deficit_rate,'edgecolor','none'); 
colormap(cmap)
colorbar
set(gca,'fontsize',15)
title('Slip Deficit Rate')
view(2)
hold on
plot_coast(bbox)
ylim([-43 -37])
xlim([172.5 180])


plot([zone1(:,1);zone1(1,1)],[zone1(:,2);zone1(1,2)],'b','linewidth',3)
plot([zone2(:,1);zone2(1,1)],[zone2(:,2);zone2(1,2)],'b','linewidth',3)
plot([zone2(:,1);zone2(1,1)],[zone2(:,2);zone2(1,2)],'b','linewidth',3)
plot([zone3(:,1);zone3(1,1)],[zone3(:,2);zone3(1,2)],'b','linewidth',3)
plot([zone4(:,1);zone4(1,1)],[zone4(:,2);zone4(1,2)],'b','linewidth',3)



figure
h=trisurf(el,nd_ll(:,1),nd_ll(:,2),nd(:,3),mean(Ring_Taus,2),'edgecolor','none'); 
colormap(flipud(pink))
colorbar
set(gca,'fontsize',15)
title('Mean Ring Stress')
view(2)
hold on
plot_coast(bbox)
ylim([-43 -37])
xlim([172.5 180])



%plot fit to data
Nbase =size(baselines,1);
dhat_base = dhat(:,1:Nbase);
dhat_u = dhat(:,1+Nbase:end);

sites_ll = local2llh(xysites',origin)';
vert_ll = local2llh(vert_xy',origin)';


%variance reduction
%ignore data in volcanic zone
BaseEnds_ll = [local2llh(BaseEnds(:,1:2)',origin)' local2llh(BaseEnds(:,3:4)',origin)'];
ind1 = BaseEnds_ll(:,1)>175.6 & BaseEnds_ll(:,1)<176.8 & BaseEnds_ll(:,2)>-38.93 & BaseEnds_ll(:,2)<-37.7;
ind2 = BaseEnds_ll(:,1)>175.3 & BaseEnds_ll(:,1)<175.9 & BaseEnds_ll(:,2)>-39.5 & BaseEnds_ll(:,2)<-38.9;
ind = ind1 | ind2;
Vbase_hat = mean(dhat_base)';
var_reduction_hz = 1 - norm(Vbase(~ind)-Vbase_hat(~ind))^2/norm(Vbase(~ind))^2
R=Vbase(~ind)./Sigbase(~ind)-Vbase_hat(~ind)./Sigbase(~ind);
chi2_hz = R'*R/length(R)

ind = vert_ll(:,1)>175.6 & vert_ll(:,1)<176.8 & vert_ll(:,2)>-38.93 & vert_ll(:,2)<-37.7;
Vu_hat = mean(dhat_u)';
var_reduction_vert = 1 - norm(Vu(~ind)-Vu_hat(~ind))^2/norm(Vu(~ind))^2
R=Vu(~ind)./Sigu(~ind)-Vu_hat(~ind)./Sigu(~ind);
chi2_vert = R'*R/length(R)



figure; 
subplot(121)
hold on;  
axis equal
for k=1:size(BaseEnds,1)
    cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k)/L(k) Vbase(k)/L(k)])
end
title('Baseline elongation rates')
plot_coast_xy(bbox,origin,'k')
load cmap
colormap(cmap)
colorbar
caxis([-max(abs(Vbase./L)) max(abs(Vbase./L))])
xlim([-500 300])
ylim([-300 300])


subplot(122)
hold on;  
axis equal
v = mean(dhat);
for k=1:size(BaseEnds,1)
    cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[v(k)/L(k) v(k)/L(k)])
end
title('Baseline elongation rates')
plot_coast_xy(bbox,origin,'k')
load cmap
colormap(cmap)
colorbar
caxis([-max(abs(Vbase./L)) max(abs(Vbase./L))])
xlim([-500 300])
ylim([-300 300])



figure; 
subplot(131)
hold on;  
axis equal
for k=1:size(BaseEnds,1)
    cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[Vbase(k)/L(k) Vbase(k)/L(k)])
end
title('Baseline elongation rates')
plot_coast_xy(bbox,origin,'k')
load cmap
colormap(cmap)
colorbar
caxis(0.5*[-max(abs(Vbase./L)) max(abs(Vbase./L))])
xlim([-500 300])
ylim([-300 300])


set(gca,'fontsize',15)

subplot(132)
hold on;  
axis equal
v = mean(dhat);
for k=1:size(BaseEnds,1)
    cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[v(k)/L(k) v(k)/L(k)])
end
title('Baseline elongation rates')
plot_coast_xy(bbox,origin,'k')
load cmap
colormap(cmap)
colorbar
caxis(0.5*[-max(abs(Vbase./L)) max(abs(Vbase./L))])
xlim([-500 300])
ylim([-300 300])


set(gca,'fontsize',15)


subplot(133)
hold on;  
axis equal
v = mean(dhat);
for k=1:size(BaseEnds,1)
    cline([BaseEnds(k,1) BaseEnds(k,3)],[BaseEnds(k,2) BaseEnds(k,4)],[(Vbase(k)-v(k))/L(k) (Vbase(k)-v(k))/L(k)])
end
title('Residual')
plot_coast_xy(bbox,origin,'k')
load cmap
colormap(cmap)
colorbar
caxis(0.5*[-max(abs(Vbase./L)) max(abs(Vbase./L))])
xlim([-500 300])
ylim([-300 300])


set(gca,'fontsize',15)

set(gcf,'renderer','Painters')




% % %plot vertical
vert_ll = local2llh(vert_xy',origin)';

figure; 
subplot(121)
hold on
scatter(vert_ll(:,1),vert_ll(:,2),50,Vu,'fill')
title('Observed Vertical')
plot_coast(bbox)
ylim([-43 -37])
xlim([171 179])
colormap(cmap)
colorbar
grid on
set(gca,'fontsize',15)
caxis([-5 5])

subplot(122)
hold on
scatter(vert_ll(:,1),vert_ll(:,2),50,mean(dhat_u)','fill')
title('Predicted Vertical')
colormap(cmap)
colorbar
grid on
set(gca,'fontsize',15)
caxis([-5 5])
plot_coast(bbox)
ylim([-43 -37])
xlim([171 179])

set(gcf,'Position',[39   243   952   387])



%% compute moment rates
A = patch_stuff.area_faces*10^6;  %convert km^2 to m^2
sr = (srate)-mean(creep_rates)'/1000; % m/yr
Mo_total = 30e9*sum(sr.*A)
Mw_500_total = (2/3) * log10(500*Mo_total*1e7) - 10.7

locking_prob = mean(locked_index)';
Mo_locking = 30e9*sum(sr.*A.*locking_prob)
Mw_500_locking = (2/3) * log10(500*Mo_locking*1e7) - 10.7

%make histograms of moments
sr = repmat(srate,1,size(creep_rates,1))-creep_rates'/1000; % m/yr
Mo_total = 30e9*sum(sr.*repmat(A,1,size(creep_rates,1)));

Mo_locking = 30e9*sum(sr.*repmat(A,1,size(creep_rates,1)).*locked_index');
Mw_500_locking = (2/3) * log10(500*Mo_locking*1e7) - 10.7



figure
histogram(Mo_total)
hold on
histogram(Mo_locking)

xlabel('Moment rate or Mo [N·m/yr]')
ylabel('Counts')
legend('total','locked areas')
set(gca,'FontSize',15)

% ---- Add second x-axis on top for Mw ----
T = 500;
% Mo is moment *rate* in N·m/yr; convert to 500-yr moment (N·m) then to Mw.
% Using Hanks & Kanamori with dyne-cm: Mw = (2/3)(log10(M0[dyncm]) - 10.7)
Mo2Mw = @(Mo) (2/3) * log10(T * Mo * 1e7) - 10.7;

ax1 = gca;

% Make sure Mo axis only draws ticks at the bottom and no top box line
ax1.XAxisLocation = 'bottom';
ax1.Box = 'off';
ax1.Layer = 'top';    % draw ticks above the bars
drawnow

% Create a transparent axes on top
ax2 = axes('Position', ax1.Position, ...
           'Color', 'none', ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right');

% Keep y limits synced and hide y-ticks on the top axes
linkaxes([ax1 ax2], 'y');
ax2.YTick = [];

% Set Mw limits to match Mo limits via the transform
ax2.XLim = Mo2Mw(ax1.XLim);

% Label + styling
xlabel(ax2, 'Equivalent 500-year Mw')
set(ax2,'FontSize',15,'XColor','k')
