%% all figures regarding subtelomeric underreplication
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

%%
D = DS;
unq_dist_to_the_end = [0.2 : 0.2 : 765];
median_unrep_cdc20 = NaN(length(unq_dist_to_the_end) , 1);
median_unrep_dbf2 = NaN(length(unq_dist_to_the_end) , 1);
for I = 1:length(unq_dist_to_the_end)
    idx = find(D.dist_to_the_end_kb >= unq_dist_to_the_end(I)-0.1 & D.dist_to_the_end_kb <= unq_dist_to_the_end(I)+0.1);
    median_unrep_cdc20(I) = nanmedian(D.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
    median_unrep_dbf2(I) = nanmedian(D .percent_unreplicated_not_trimmed_dbf2_smooth(idx));
end
DS.percent_underreplicated_cdc20_not_trimmed_DM_dist = NaN(length(DS) , 1);
DS.percent_underreplicated_dbf2_not_trimmed_DM_dist = NaN(length(DS) , 1);
for I = 1:length(DS)
    vq = interp1(unq_dist_to_the_end , median_unrep_cdc20 , DS.dist_to_the_end_kb(I) , 'spline');
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(I) = DS.percent_unreplicated_not_trimmed_cdc20(I) - vq;
    vq = interp1(unq_dist_to_the_end, median_unrep_dbf2 , DS.dist_to_the_end_kb(I) ,  'spline');
    DS.percent_underreplicated_dbf2_not_trimmed_DM_dist(I) = DS.percent_unreplicated_not_trimmed_dbf2(I) - vq;
end
save('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat','DS');
%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__features.mat' , 'DS');
D = DS;
idx = find(strcmp(D.TYPE , 'ARS'));
D = D(idx , :);
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat' , 'DS');
DS.dist_to_ARS = NaN(length(DS) , 1);
for I = 1:length(DS)
    idx = find(D.chr_num == DS.chr_num(I));
    data = D.middle_point_kb(idx);
    DS.dist_to_ARS(I) = nanmin(abs(DS.middle_point_kb(I) - data));
end
save('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat','DS');

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

figure; 
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
%DS = DS(idx , :);
%
subplot(3,2,1); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Distance to the end');
ylim([-50 100]);
%ylabel('Degree Underreplication');

subplot(3,2,2); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_the_end_kb(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) );
[C,p] = corrcoef(DS.dist_to_the_end_kb(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Distance to the end');
ylim([-50 50]);
%ylabel('Degree Underreplication DM');

subplot(3,2,3); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.Trep_spline(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.Trep_spline(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Replication timing, minutes after release from alpha-factor arrest');
ylim([-50 100]);
%ylabel('Degree Underreplication');
xlim([0 65]);

subplot(3,2,4); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.Trep_spline(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) );
[C,p] = corrcoef(DS.Trep_spline(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Replication timing');
%ylabel('Degree Underreplication DM');
xlim([0 65]);
ylim([-100 50]);

subplot(3,2,5); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_ARS(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.dist_to_ARS(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
ylim([-50 100]);
xlabel('Distance to ARS');
%ylabel('Degree Underreplication');
xlim([0 45]);

subplot(3,2,6); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_ARS(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
[C,p] = corrcoef(DS.dist_to_ARS(idx) , DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Distance to ARS');
%ylabel('Degree Underreplication DM');
xlim([0 45]);
ylim([-100 50]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig6A');

%%
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');

figure; 
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
%DS = DS(idx , :);
%
subplot(3,1,1); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.dist_to_the_end_kb(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('KB from the telomere');
ylim([-50 100]);
colorbar;
%ylabel('Degree Underreplication');

subplot(3,1,2); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.Trep_spline(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.Trep_spline(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
xlabel('Replication timing, minutes after release from alpha-factor arrest');
ylim([-50 100]);
%ylabel('Degree Underreplication');
xlim([0 65]);
colorbar;

subplot(3,1,3); hold on; grid on; set(gca , 'FontSize' , 12);
dscatter(DS.dist_to_ARS(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx) );
[C,p] = corrcoef(DS.dist_to_ARS(idx) , DS.percent_unreplicated_not_trimmed_cdc20_smooth(idx));
title(sprintf('R = %.2f, P = %.2f.' , C(1,2) , p(1,2)));
ylim([-50 100]);
xlabel('Distance to the closest ARS, KB');
%ylabel('Degree Underreplication');
xlim([0 45]);
colorbar;

print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig6A');


%% RT
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp.mat');
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
DS = DS(idx , :);
figure; hold on; grid on; set(gca , 'FontSize' , 12);
data = DS.Trep_spline;
% trep_grid = [0 quantile(data,0.9) ...
%     quantile(data,0.95) quantile(data,0.99) quantile(data,0.995) quantile(data,1)];
trep_grid = [0 40 45 50 70];
data = NaN(length(DS) , length(trep_grid)-1);
clrs = winter(4*length(trep_grid)-4);
for I = 1:length(trep_grid)-1
    idx = find(DS.Trep_spline > trep_grid(I) & DS.Trep_spline < trep_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    scatter(repmat(I,length(idx),1)+randn(length(idx),1)*0.03-0.015 , ...
       data(1:length(idx) , I) , 10 , clrs(4*I-2,:) , 'filled');
end
h = boxplot(data , 'color' , [.2 .2 .2] , 'notch' , 'on' , 'symbol' , '');
set(h , 'LineWidth' , 1.5);
ylim([-35 40]);
set(gca , 'Xtick' , []);
set(gca , 'Xtick' , [1:length(trep_grid)-1] , 'XtickLabel' , ...
    {'<= 40 min' , '40-45 min' , '45-50 min' , '> 50 min'} );
xlabel('Replication timing, min after \alpha factor release');
ylabel('underreplication DM');
set(gcf, 'PaperPosition' , [0 0 12 12]);
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig6B');


%% legend for 6B
figure; hold on; grid on; set(gca , 'FontSize' , 20);
clrs = winter(4*length(trep_grid)-4);
for I = 1:length(trep_grid)-1
    idx = find(DS.Trep_spline > trep_grid(I) & DS.Trep_spline < trep_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    scatter(repmat(I,length(idx),1)+randn(length(idx),1)*0.03-0.015 , ...
       data(1:length(idx) , I) , 10 , clrs(4*I-2,:) , 'filled');
end
legend('<= 0.9' , '0.9-0.95' , '0.95-0.99' , '0.99-0.995' , '> 0.995',...
    'location' , 'best');

%% dist to ARS
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
DS = DS(idx , :);
figure; hold on; grid on; set(gca , 'FontSize' , 16);
data = DS.dist_to_ARS;
ARS_grid = [0 ...
    quantile(data,0.999) quantile(data,1)];
data = NaN(length(DS) , length(ARS_grid)-1);
clrs = spring(4*length(ARS_grid)-4);
for I = 1:length(ARS_grid)-1
    idx = find(DS.dist_to_ARS > ARS_grid(I) & DS.dist_to_ARS < ARS_grid(I+1) );
    data(1:length(idx) , I) = DS.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx);
    scatter(repmat(I,length(idx),1)+randn(length(idx),1)*0.03-0.015 , ...
       data(1:length(idx) , I) , 10 , clrs(4*I-2,:) , 'filled');
end
boxplot(data , 'color' , [.2 .2 .2] , 'notch' , 'on' , 'symbol' , '');
ylim([-30 30]);
set(gca , 'Xtick' , [1:length(ARS_grid)-1] , 'XtickLabel' , ...
    {'<99.9%' , '>99.9%'});
xlabel('Quantile of distance to ARS');
ylabel('Degree underreplication DM');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/Manuscript/Fig6C');



%%
idx = find(~isnan(DS.percent_unreplicated_not_trimmed_cdc20_smooth) & ...
    DS.percent_unreplicated_not_trimmed_cdc20_smooth < 10000 & ...
	DS.percent_unreplicated_not_trimmed_cdc20_smooth > -10000 & ...
    ~isnan(DS.percent_underreplicated_cdc20_not_trimmed_DM_dist) & ...
    DS.percent_underreplicated_cdc20_not_trimmed_DM_dist < 10000 & ...
	DS.percent_underreplicated_cdc20_not_trimmed_DM_dist > -10000  );
D = DS(idx , :);
figure; 
clrs1 = parula(12);
clrs2 = bone(12);

subplot(2,1,1);hold on; grid on; set(gca , 'FontSize' , 12);
idx = find(D.Trep_spline >= quantile(D.Trep_spline , 0.99));
[y,x] = ksdensity(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , [-30:2:40]);
plot(x,y/sum(y) , 'LineWidth' , 3 , 'Display' , 'Top 1% by replication time' , 'color' , clrs1(1,:));
idx = find(D.Trep_spline < quantile(D.Trep_spline , 0.99));
[y,x] =  ksdensity(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , [-30:2:40] );
plot(x,y/sum(y) , 'LineWidth' , 3 , 'Display' , 'Bottom 99% by replication time','color' , clrs1(4,:));
xlabel('Degree underreplication DM');
legend('location' , 'best');
ylabel('PDF');
[~,p] = ttest2(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(D.Trep_spline >= quantile(D.Trep_spline , 0.99)),...
D.percent_underreplicated_cdc20_not_trimmed_DM_dist(D.Trep_spline < quantile(D.Trep_spline , 0.99)) );
title(sprintf('P = %.2f.' , p));


subplot(2,1,2);hold on; grid on; set(gca , 'FontSize' , 12);
idx = find(D.Trep_spline >= quantile(D.dist_to_ARS , 0.99));
[y,x] = ksdensity(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , [-30:2:40]);
plot(x,y/sum(y) , 'LineWidth' , 3 , 'Display' , 'Top 1% distance to ARS','color' , clrs2(2,:));
idx = find(D.Trep_spline < quantile(D.dist_to_ARS , 0.99));
[y,x] =  ksdensity(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(idx) , [-30:2:40] );
plot(x,y/sum(y) , 'LineWidth' , 3 , 'Display' , 'Bottom 99% distance to ARS','color' , clrs2(7,:));
xlabel('Degree underreplication DM');
legend('location' , 'best');
ylabel('PDF');
[~,p] = ttest2(D.percent_underreplicated_cdc20_not_trimmed_DM_dist(D.dist_to_ARS >= quantile(D.dist_to_ARS , 0.99)),...
D.percent_underreplicated_cdc20_not_trimmed_DM_dist(D.dist_to_ARS < quantile(D.dist_to_ARS , 0.99)) );
title(sprintf('P = %.2f.' , p));








