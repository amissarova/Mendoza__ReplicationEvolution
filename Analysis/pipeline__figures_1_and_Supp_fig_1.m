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





