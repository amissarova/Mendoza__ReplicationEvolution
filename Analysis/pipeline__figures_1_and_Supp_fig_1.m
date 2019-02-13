%% pipeline__figures_1_and_Supp_fig_1
addpath(genpath('~/Develop/matlab'));

%% MAIN fig 1C  
    
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
        if current_pair_to_plot(J) < 7 | current_pair_to_plot(J) > 10
            idx_max = find(data == 300);
            for Z = 1:length(idx_max)
                data(idx_max(Z)) = 10000000;
            end
        else
            idx_max = find(data == 200);
            for Z = 1:length(idx_max)
                data(idx_max(Z)) = 10000000;
            end
        end
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
    print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r600');
end
%%
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
        ylim([0 37]);
    elseif I == 2
        ylim([0 75]);
    else
        ylim([0 75]);
    end
    save_name = sprintf('1c__boxplots__ChrBridge__main__subpanel_%d' , I);
    print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r600');
end

%% Supp fig 1C-1

D_ND = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ND.tab');
D_ChrBridge = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ChrBridge.tab');

N = 3;
% 
sets_to_plot = cell(N,1);
sets_to_plot{1} = [3 4]; sets_to_plot{2} = [5 6]; sets_to_plot{3} = [1 2]; 

clrs1 = spring(12); clrs2 = winter(12); clrs3 = parula(12); clrs4 = summer(12); clrs5 = hot(12); clrs6 = lines(6);
clrs_to_plot = cell(N,1);

clrs_to_plot{1} = [.4 .4 .4; clrs6(4,:) ];
clrs_to_plot{2} = [clrs1(6,:); clrs5(3,:)];
clrs_to_plot{3} = [clrs4(4,:); clrs3(4,:)];

figure('units','centimeters','position',[5 5 15 5]);
for I = 1:N
    subplot(1,N,I); hold on; grid on;
    current_pair_to_plot = sets_to_plot{I};
    temp_clrs_to_plot = clrs_to_plot{I};
    
    for J = 1:length(current_pair_to_plot)
        data = double(D_ND(:,current_pair_to_plot(J))); data = data(~isnan(data)); 
        if current_pair_to_plot(J) < 7 | current_pair_to_plot(J) > 10
            idx_max = find(data == 300);
            for Z = 1:length(idx_max)
                data(idx_max(Z)) = 10000000;
            end
        else
            idx_max = find(data == 200);
            for Z = 1:length(idx_max)
                data(idx_max(Z)) = 10000000;
            end
        end
        [y,x] = ecdf(data);
        if I~=2 | J ~= 2
            plot(x,y*100,'LineWidth',2.5,'color', temp_clrs_to_plot(J,:));       
        else
            plot(x,y*100,'LineWidth',2.5,'color', temp_clrs_to_plot(J,:) , 'LineStyle' , '--');     
        end
        xlim([0 250]);
    end
end
save_name = '1c__eCDFs__ND__Supp_1';
print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r300');

%%
figure('units','centimeters','position',[5 5 15 5]);
for I = 1:N
    subplot(1,N,I); hold on; grid on;
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
        ylim([0 45]);
	elseif I == 2
        ylim([0 45]);
    else
        ylim([0 45]);
	end
end
save_name = '1c__boxplots__ChrBridge__Supp_1';
print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r600');

%% Supp fig 1C-2

D_ND = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ND.tab');
D_ChrBridge = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/AM_stat_fig1C_ChrBridge.tab');

N = 2;
% 
sets_to_plot = cell(N,1);
sets_to_plot{1} = [14 1 11 12 13]; sets_to_plot{2} = [7 8]; 

clrs1 = spring(12); clrs2 = winter(12); clrs3 = parula(12); clrs4 = summer(12); clrs5 = hot(12); clrs6 = lines(6);
clrs_to_plot = cell(N,1);

clrs_to_plot{1} = [.4 .4 .4; clrs4(4,:) ; clrs5(3,:) ; clrs2(7,:) ; clrs3(10,:)];
clrs_to_plot{2} = [.4 .4 .4; clrs2(7,:)];

figure('units','centimeters','position',[5 5 12 5]);
for I = 1:N
    subplot(1,N,I); hold on; grid on;
    current_pair_to_plot = sets_to_plot{I};
    temp_clrs_to_plot = clrs_to_plot{I};
    
    for J = 1:length(current_pair_to_plot)
        data = double(D_ND(:,current_pair_to_plot(J))); data = data(~isnan(data));
        if current_pair_to_plot(J) < 7 | current_pair_to_plot(J) > 10
            idx_max = find(data == 300);
            for Z = 1:length(idx_max)
                data(idx_max(Z)) = 10000000;
            end
        else
            idx_max = find(data == 200);
            for Z = 1:length(idx_max)
                data(idx_max(Z)) = 10000000;
            end
        end
            
        [y,x] = ecdf(data);
            plot(x,y*100,'LineWidth',2.5,'color', temp_clrs_to_plot(J,:));    
        if I == 1
            xlim([0 250]);
        else
            xlim([0 120]);
        end
    end
end
save_name = '1c__eCDFs__ND__Supp_2';
print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r300');

%%
figure('units','centimeters','position',[5 5 12 5]);
for I = 1:N
    subplot(1,N,I); hold on; grid on;
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
        ylim([0 65]);
	elseif I == 2
        ylim([0 20]);
	end
end
save_name = '1c__boxplots__ChrBridge__Supp_2';
print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r600');

%% Supp fig 1C-3
D_ND = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/1335_MM.txt');
figure('units','centimeters','position',[5 5 6 5]); hold on; grid on;
clrs1 = lines(6);
data = D_ND.Spc42_noHU;
[y,x] = ecdf(data);
plot(x,y*100,'LineWidth',2.5,'color', [.4 .4 .4]);

data = D_ND.Spc42_HU;
[y,x] = ecdf(data);
plot(x,y*100,'LineWidth',2.5,'color', clrs1(4,:));  

xlim([0 120]);
save_name = '1c__eCDFs__ND__Supp_3_Spc42';
print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r600');

%
D_ND = dataset('file' , '~/Develop/Mendoza__ReplicationEvolution/figure1_data_analysis/5066_MM.txt');
figure('units','centimeters','position',[5 5 6 5]); hold on; grid on;
clrs1 = lines(6);
data = D_ND.Htb2_noHU;
[y,x] = ecdf(data);
plot(x,y*100,'LineWidth',2.5,'color', [.4 .4 .4]);

data = D_ND.Htb2_HU_0308;
[y,x] = ecdf(data);
plot(x,y*100,'LineWidth',2.5,'color', clrs1(4,:));  

xlim([0 120]);
save_name = '1c__eCDFs__ND__Supp_3_Htb2';
print('-dpng' , strcat('~/Develop/Mendoza__ReplicationEvolution/Figures/fig1/', save_name ) , '-r600');

%%





