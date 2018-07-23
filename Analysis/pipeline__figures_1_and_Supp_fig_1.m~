%% pipeline__figures_1_and_Supp_fig_1
addpath(genpath('~/Develop/matlab'));

%% 1. Boxplot -- Chromatin bridges duration in psf2 and psf2rad9
% option -- boxplot for each mutant/replicate
D = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/1 - clean data chromatin bridge by replicates.txt');

clrs1 = winter(12); clrs_psf2 = [clrs1(1,:) ; clrs1(4,:) ; clrs1(6,:) ; clrs1(8,:)];
clrs2 = spring(12); clrs_psf2rad9 = [clrs2(5,:) ; clrs2(7,:)];

data = double(D(:,:));
figure('units','centimeters','position',[5 5 12 10]);
hold on; grid on;
h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
N = length(h);
for j=1:4
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_psf2(j,:) ,'FaceAlpha',.7 );
end
for j=5:6
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_psf2rad9(j-4,:) ,'FaceAlpha',.7 );
end
set(gca , 'Xtick' , []);
ylim([0 60]);
ylabel('Duration of chromatin bridge, minutes');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/boxplot__chromatin_bridges_psf2_psf2rad9__by_reps' , '-r300');


% option -- boxplot for pooled data
D = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/1 - clean data chromatin bridge by replicates.txt');

clrs1 = winter(12); clrs_psf2 = [clrs1(6,:) ];
clrs2 = spring(12); clrs_psf2rad9 = [clrs2(6,:)];

psf2 = [D.psf2_1 ; D.psf2_2 ; D.psf2_3 ; D.psf2_4]; psf2 = psf2(~isnan(psf2));
psf2rad9 = [D.psf2rad9_1 ; D.psf2rad9_2]; psf2rad9 = psf2rad9(~isnan(psf2rad9));

N = max([length(psf2) length(psf2rad9)]);
data = NaN(N,2);
data(1:length(psf2) , 1) = psf2; data(1:length(psf2rad9) , 2) = psf2rad9;
figure('units','centimeters','position',[5 5 12 10]);
hold on; grid on;
h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
N = length(h);
for j=1:1
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_psf2(j,:) ,'FaceAlpha',.7 );
end
for j=2:2
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_psf2rad9(j-1,:) ,'FaceAlpha',.7 );
end
set(gca , 'Xtick' , []);
ylim([0 60]);
ylabel('Duration of chromatin bridge, minutes');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/boxplot__chromatin_bridges_psf2_psf2rad9__pooled' , '-r300');


%% 2. Boxplot -- Nuclear in psf2 and psf2rad9
D = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/1 - clean data nuclear division by replicates.txt');

psf2 = [D.psf2_1 ; D.psf2_2 ; D.psf2_3 ; D.psf2_4]; psf2 = psf2(~isnan(psf2));
psf2rad9 = [D.psf2rad9_1 ; D.psf2rad9_2]; psf2rad9 = psf2rad9(~isnan(psf2rad9));

data = double(D(:,:));

clrs1 = winter(12); clrs_psf2 = [clrs1(1,:) ; clrs1(4,:) ; clrs1(6,:) ; clrs1(8,:)];
clrs2 = spring(12); clrs_psf2rad9 = [clrs2(5,:) ; clrs2(7,:)];

figure; hold on; grid on;
h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , '');
set(h1 , 'LineWidth' , 1.5);
h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
N = length(h);
for j=1:4
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_psf2(j,:) ,'FaceAlpha',.7 );
end
for j=5:6
    patch(get(h(N-j+1),'XData'),get(h(N-j+1),'YData'), clrs_psf2rad9(j-4,:) ,'FaceAlpha',.7 );
end
set(gca , 'Xtick' , []);
ylim([0 301]);
ylabel('Time until nuclear division, minutes');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/boxplot__nuclear_division_psf2_psf2rad9' , '-r300');

%% 3. eCDFs -- Nuclear in psf2 and psf2rad9
% option -- ecdf for each mutant/replicate

D = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/1 - clean data nuclear division by replicates.txt');

clrs1 = winter(12); clrs_psf2 = [clrs1(1,:) ; clrs1(4,:) ; clrs1(6,:) ; clrs1(8,:)];
clrs2 = spring(12); clrs_psf2rad9 = [clrs2(5,:) ; clrs2(7,:)];

psf2 = cell(4,1); psf2{1} = D.psf2_1; psf2{2} = D.psf2_2; psf2{3} = D.psf2_3; psf2{4} = D.psf2_4;
psf2rad9 = cell(2,1); psf2rad9{1} = D.psf2rad9_1; psf2rad9{2} = D.psf2rad9_2;
figure('units','centimeters','position',[5 5 12 10]);
hold on; grid on;
for I = 1:4
    data = psf2{I}; data = data(~isnan(data));
    [y,x] = ecdf(data);
    plot(x,y , 'LineWidth' , 2.5 , 'color' , clrs_psf2(I,:) , 'Display' , ...
        strcat('{\Delta}\it{psf2} -' ,  sprintf('rep. %d',I) ));
end
for I = 1:2
    data = psf2rad9{I}; data = data(~isnan(data));
    [y,x] = ecdf(data);
    plot(x,y , 'LineWidth' , 2.5 , 'color' , clrs_psf2rad9(I,:), 'Display' , ...
        strcat('{\Delta}\it{psf2}{\Delta}\it{rad9} -' ,  sprintf('rep. %d',I) ));
end
xlim([0 250]);
legend('location' , 'se');
ylabel('eCDF');
xlabel('Time before N.D., minutes');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/ecdf__nuclear_division_psf2_psf2rad9__by_reps' , '-r300');

% option -- ecdf for pooled data
D = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/1 - clean data nuclear division by replicates.txt');

clrs1 = winter(12); clrs_psf2 = [clrs1(6,:) ];
clrs2 = spring(12); clrs_psf2rad9 = [clrs2(6,:)];

psf2 = [D.psf2_1 ; D.psf2_2 ; D.psf2_3 ; D.psf2_4]; psf2 = psf2(~isnan(psf2));
psf2rad9 = [D.psf2rad9_1 ; D.psf2rad9_2]; psf2rad9 = psf2rad9(~isnan(psf2rad9));

figure('units','centimeters','position',[5 5 12 10]);
hold on; grid on;
    
data = psf2; data = data(~isnan(data));
[y,x] = ecdf(data);
plot(x,y , 'LineWidth' , 2.5 , 'color' , clrs_psf2(1,:) , 'Display' , ...
        '{\Delta}\it{psf2}' );

data = psf2rad9; data = data(~isnan(data));
[y,x] = ecdf(data);
plot(x,y , 'LineWidth' , 2.5 , 'color' , clrs_psf2rad9(1,:), 'Display' , ...
        '{\Delta}\it{psf2}{\Delta}\it{rad9}' );

xlim([0 250]);
legend('location' , 'se');
ylabel('eCDF');
xlabel('Time before N.D., minutes');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/ecdf__nuclear_division_psf2_psf2rad9__pooled' , '-r300');

%%
D_ND = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ND.tab');
D_ChrBridge = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ChrBridge.tab');

N = 9;
% 
pairs_to_plot = [3 4; 4 6; 1 2; 7 8; 9 10 ; 5 6 ; 3 11 ; 3 12; 3 13 ;  ];
%sets_to_plot = [3 4 1 8];

clrs1 = spring(12); clrs2 = winter(12); clrs3 = parula(12); clrs4 = summer(12);
clrs_to_plot = cell(N,1);
clrs_to_plot{1} = [.6 .6 .6; clrs2(6,:)];
clrs_to_plot{2} = [clrs2(6,:); clrs1(5,:)];
clrs_to_plot{3} = [clrs2(2,:); clrs1(7,:)];
clrs_to_plot{4} = [.6 .6 .6; clrs2(9,:)];
clrs_to_plot{5} = [.6 .6 .6; clrs2(9,:)];
clrs_to_plot{6} = [clrs4(2,:); clrs1(5,:)];
clrs_to_plot{7} = [.6 .6 .6; clrs3(4,:)];
clrs_to_plot{8} = [.6 .6 .6; clrs2(9,:)];
clrs_to_plot{9} = [.6 .6 .6; clrs3(7,:)];

for I = 1:N
    figure('units','centimeters','position',[5 5 5 5]); hold on; grid on;
    current_pair_to_plot = pairs_to_plot(I,:);
    
    temp_clrs_to_plot = clrs_to_plot{I};
    data = double(D_ND(:,current_pair_to_plot(1))); data = data(~isnan(data)); 
    [y,x] = ecdf(data);
    if I ~= 6
        plot(x,y*100,'LineWidth',2.5,'color', temp_clrs_to_plot(1,:));
    else
        plot(x,y*100,'LineWidth',4.5,'color', temp_clrs_to_plot(1,:) , 'LineStyle' , '-');
    end
    
    data = double(D_ND(:,current_pair_to_plot(2))); data = data(~isnan(data)); 
    [y,x] = ecdf(data);
    plot(x,y*100,'LineWidth',2.5,'color', temp_clrs_to_plot(2,:));
    if I < 4 | I > 5
        xlim([0 250]);
    else
        xlim([0 120]);
    end
    save_name = sprintf('1c__eCDFs__ND__vars_%d_%d' , current_pair_to_plot(1) , current_pair_to_plot(2));
    print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r300');
end

%
for I = 1:N
    figure('units','centimeters','position',[5 5 5 5]); hold on; grid on;
    current_pair_to_plot = pairs_to_plot(I,:);
    
    temp_clrs_to_plot = clrs_to_plot{I};
    data = double(D_ChrBridge(:,current_pair_to_plot)); 
	h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , ''); set(h1 , 'LineWidth' , 1.5);
	h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
	for j=1:2
        patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), temp_clrs_to_plot(j,:) ,'FaceAlpha',.8 );
    end
	set(gca , 'Xtick' , []);    
	ylim([0 50]);
    save_name = sprintf('1c__boxplots__ChrBridge__vars_%d_%d' , current_pair_to_plot(1) , current_pair_to_plot(2));
    print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r300');
end


%%
figure('units','centimeters','position',[5 5 25 5]);
for I = 1:N
    current_pair_to_plot = pairs_to_plot(I,:);
    subplot(1,N,I); hold on; grid on;
    if I <=3
        temp_clrs_to_plot = clrs_to_plot{I};
        data = double(D_ChrBridge(:,current_pair_to_plot)); 
        h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , ''); set(h1 , 'LineWidth' , 1.5);
        h = findobj(gca,'Tag','Box'); set(h1 , 'LineWidth' , 1.5);
        for j=1:2
            patch(get(h(2-j+1),'XData'),get(h(2-j+1),'YData'), temp_clrs_to_plot(j,:) ,'FaceAlpha',.8 );
        end
        set(gca , 'Xtick' , []);    
        ylim([0 45]);
    end
end
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/1c__boxplots__ChrBridges' , '-r300');

%% legend

clrs1 = spring(12); clrs2 = winter(12); 
clrs_to_plot = cell(5,1);
clrs_to_plot{1} = [.6 .6 .6; clrs2(6,:)];
clrs_to_plot{2} = [clrs2(6,:); clrs1(5,:)];
clrs_to_plot{3} = [clrs2(2,:); clrs1(7,:)];
clrs_to_plot{4} = [.6 .6 .6; clrs2(9,:)];
clrs_to_plot{5} = [.6 .6 .6; clrs2(9,:)];

clrs_set = [.6 .6 .6; clrs2(6,:) ; clrs2(2,:) ; clrs2(9,:) ; clrs1(5,:) ; clrs1(7,:)];
legend_titles = {'Untreated' , '+HU' , '\it{psf2ts}' , '\it{pol2-12}' , ...
    '\it{rad9}{\Delta}+HU' , '\it{rad9}{\Delta}\it{psf2ts}'};
figure; hold on; set(gca , 'FontSize' , 20);
for I = 1:6
    plot([1 1 ] , [2 2 ] , 'LineWidth' , 4 , 'color' , clrs_set(I,:) , 'Display' , legend_titles{I})
end
legend('location' , 'best');
print('-dpng' , '~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/1c__legend' , '-r300');
    
%% MAIN fig    
    
D_ND = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ND.tab');
D_ChrBridge = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ChrBridge.tab');

N = 4;
% 
sets_to_plot = cell(N,1);
sets_to_plot{1} = [3 4]; sets_to_plot{2} = [14 1 12]; sets_to_plot{3} = [5 6]; sets_to_plot{4} = [9 10]; 

clrs1 = spring(12); clrs2 = winter(12); clrs3 = parula(12); clrs4 = summer(12); clrs5 = hot(12); clrs6 = lines(6);
clrs_to_plot = cell(N,1);

clrs_to_plot{1} = [.4 .4 .4; clrs6(4,:)];
clrs_to_plot{2} = [.4 .4 .4; clrs4(4,:); clrs2(7,:)];
clrs_to_plot{3} = [.4 .4 .4; clrs1(7,:)];
clrs_to_plot{4} = [.4 .4 .4; clrs2(7,:)];

for I = 1:N
    figure('units','centimeters','position',[5 5 5 5]);
    hold on; grid on;
    current_pair_to_plot = sets_to_plot{I};
    temp_clrs_to_plot = clrs_to_plot{I};
    
    for J = 1:length(current_pair_to_plot)
        data = double(D_ND(:,current_pair_to_plot(J))); data = data(~isnan(data)); 
        [y,x] = ecdf(data);
        if J == 1
            plot(x,y*100,'LineWidth',3.5,'color', temp_clrs_to_plot(J,:));
        else
            plot(x,y*100,'LineWidth',2.5,'color', temp_clrs_to_plot(J,:));
        end
        if I <= 3
            xlim([0 250]);
        else
            xlim([0 120]);
        end
    end
    
    save_name = sprintf('1c__eCDFs__ND__main__subpanel_%d' , I);
    print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r300');
end

%
for I = 1:3
    figure('units','centimeters','position',[5 5 5 5]); hold on; grid on;
    current_pair_to_plot = sets_to_plot{I};
    temp_clrs_to_plot = clrs_to_plot{I};
    
    data = double(D_ChrBridge(:,current_pair_to_plot)); 
	h1 = boxplot(data , 'color' , [.2 .2 .2] , 'symbol' , ''); set(h1 , 'LineWidth' , 1.3);
	h = findobj(gca,'Tag','Box'); 
    K = length(h);
	for j=1:K
        patch(get(h(K-j+1),'XData'),get(h(K-j+1),'YData'), temp_clrs_to_plot(j,:) ,'FaceAlpha',.8 );
    end
	set(gca , 'Xtick' , []);  
    if I == 1
        ylim([0 75]);
    elseif I == 2
        ylim([0 75]);
    else
        ylim([0 75]);
    end
    save_name = sprintf('1c__boxplots__ChrBridge__main__subpanel_%d' , I);
    print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r300');
end

%% Supp fig
D_ND = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ND.tab');
D_ChrBridge = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ChrBridge.tab');

N = 3;
% 
sets_to_plot = cell(N,1);
sets_to_plot{1} = [3 4 1 11 12 13]; sets_to_plot{2} = [1 2 4 6]; sets_to_plot{3} = [7 8]; 

clrs1 = spring(12); clrs2 = winter(12); clrs3 = parula(12); clrs4 = summer(12); clrs5 = hot(12); clrs6 = lines(6);
clrs_to_plot = cell(N,1);

clrs_to_plot{1} = [.4 .4 .4; clrs6(4,:) ; clrs4(4,:) ; clrs5(3,:) ; clrs2(7,:) ; clrs3(10,:)];
clrs_to_plot{2} = [clrs4(4,:); clrs1(3,:) ; clrs6(4,:) ; clrs1(7,:)];
clrs_to_plot{3} = [.4 .4 .4; clrs2(7,:)];

figure('units','centimeters','position',[5 5 15 5]);
for I = 1:N
    subplot(1,N,I); hold on; grid on;
    current_pair_to_plot = sets_to_plot{I};
    temp_clrs_to_plot = clrs_to_plot{I};
    
    for J = 1:length(current_pair_to_plot)
        data = double(D_ND(:,current_pair_to_plot(J))); data = data(~isnan(data)); 
        [y,x] = ecdf(data);
        if I ~= 2 & J == 1
            plot(x,y*100,'LineWidth',3.5,'color', temp_clrs_to_plot(J,:));
        else
            plot(x,y*100,'LineWidth',2.5,'color', temp_clrs_to_plot(J,:));
        end
        if I < 3
            xlim([0 250]);
        else
            xlim([0 120]);
        end
    end
    
    save_name = '1c__eCDFs__ND__Supp';
    print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r300');
end
D_ND = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ND.tab');
D_ChrBridge = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ChrBridge.tab');

