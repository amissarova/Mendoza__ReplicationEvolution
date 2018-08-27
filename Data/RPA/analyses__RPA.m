
%% get RPA data
addpath(genpath('~/Develop/matlab'));
cd ~/Develop/Mendoza__ReplicationEvolution/RPA

%%
D = dataset('file' , 'RPA.tab');

%%
clrs1 = lines(6); clrs2 = summer(12); clrs_set = [clrs1(4,:) ; clrs2(4,:)];
legend_titles = {'no RPA foci' , 'RPA foci'};
figure('units','centimeters','position',[5 5 7 10]); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(D.MITOSIS , D.Category , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); 
N = length(h);
for j=1:N
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_set(j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
legend('location' , 'nw');
ylim([0 40]);
set(gca , 'Xtick' , '');
%print('-dpng' , 'mitosis' , '-r300');

clrs1 = lines(6); clrs2 = summer(12); clrs_set = [clrs1(4,:) ; clrs2(4,:)];
legend_titles = {'no RPA foci' , 'RPA foci'};
figure('units','centimeters','position',[5 5 7 10]); hold on; grid on; set(gca , 'FontSize' , 10);
h1 = boxplot(D.ChrBridge , D.Category , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); 
N = length(h);
for j=1:N
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_set(j,:) ,'FaceAlpha',.5 , 'Display' , legend_titles{j});
end
legend('location' , 'nw');
ylim([0 25]);
set(gca , 'Xtick' , '');
%print('-dpng' , 'chr_bridge' , '-r300');