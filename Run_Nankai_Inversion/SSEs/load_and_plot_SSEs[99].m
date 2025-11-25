clear; close all;

load Nankai_LSSEs.mat
load Nankai_SSSEs.mat

figure
hold on

%% PLOT S-SSEs
% plots Kano and Kato S-SSEs
for i = 1:length(kano_count_vals)
    start_ind = kano_label_ind(i)+1;
    end_ind = kano_label_ind(i)+kano_count_vals(i)-1;
    if kano_count_vals(i)>5
        fill([kano_sse_cont(start_ind:end_ind,1); kano_sse_cont(start_ind,1)],[kano_sse_cont(start_ind:end_ind,2); kano_sse_cont(start_ind,2)],'m','EdgeColor','none')
    end
end

% plots individual S-SSE offshore Kii Peninsula
scatter(136.5,33.11,52*50,'MarkerEdgeColor','m','MarkerFaceColor','m','LineWidth',3)


%% PLOT L-SSEs
% plots L-SSEs from Takagi et al.
for i = 1:length(lsse_count_vals)
    start_ind = lsse_label_ind(i)+1;
    end_ind = lsse_label_ind(i)+lsse_count_vals(i)-1;
    plot([lsse_cont(start_ind:end_ind,1); lsse_cont(start_ind,1)],[lsse_cont(start_ind:end_ind,2); lsse_cont(start_ind,2)],'Color',[255,233,53]./256,'LineWidth',3)
end

% individual circles representing further east L-SSEs
for i = 1:3
    scatter(LSSEs(i,1),LSSEs(i,2),LSSEs(i,3)*50,'MarkerEdgeColor',[255,233,53]./256,'LineWidth',3)
end


set(gca,'fontsize',24)
set(gcf,'Units','inches','Position',[0.1944,3.625,13.5417,10.3611]);
