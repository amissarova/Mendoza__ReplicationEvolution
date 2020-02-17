%% Correlation of chromosome arm length vs lenghth-of-underrep
% Lucas Febuary 2020
%  replicating : Supplementary Figure 16. Chromosome arms with longer under-replicated
% 
% use bootstrapping to calculate robustness of correlation


%% load data

FIGNAME = '~/Downloads/FigS16_bootstrapping_to_calculate_robustness_of_correlation' ; 

cd('~/Develop/Mendoza__ReplicationEvolution/');
load('Data/MutantChr200_Michi.mat');
load('Data/DS_stat__features.mat');

% get chromosome arm length by taking the max value in each chromosome arm
DS = grpstats( DS , {'chr_num' 'arm'} , 'max' ,'DataVars' , {'middle_point_kb' 'end_point'}); 
DS = DS( ~strcmp(DS.arm,'center'),:);
DS.arm_length_kb = DS.max_middle_point_kb ;
DS.arm_length_nt = DS.max_end_point ; 

% convert to tables
DS = dataset2table(DS);
G = dataset2table(G);

%% join the tables and keep only the variables that we need
Q = outerjoin( G , DS( : , {'chr_num' , 'arm_length_kb' 'arm' })  , 'Key' ,'chr_num' ,'MergeKeys',true) ;
Q = Q( : , {'Replicate' 'MutantID' 'chr_num' 'length_chr' ...
    'length_underreplicated_left' 'length_underreplicated_right' ...
    'arm_length_kb' 'arm'});

%  %underrep is repeated for each replicate. removes replicates
Q = grpstats(Q , {'MutantID' 'chr_num' 'arm' 'arm_length_kb'} , 'mean' ,'DataVars' , ...
    { 'length_underreplicated_left' 'length_underreplicated_right' });

% only want CDCD20 data
Q = Q( strcmp(Q.MutantID,'cdc20') ,:);

% convert to a tall table
Q = stack(Q , {'mean_length_underreplicated_left' 'mean_length_underreplicated_right'}) ; 
Q.cl = cellstr(Q.mean_length_underreplicated_left_mean_length_underreplicated__1);

% only chr arms that match 
idxL = (strcmp(Q.arm,'left')&regexpcmp(Q.cl,'_left$')) ;
idxR = (strcmp(Q.arm,'right')&regexpcmp(Q.cl,'_right$')) ;
Q = Q(idxL|idxR,:);

% to correctly calcualte right size, subtract the left
for I = 1:height(Q)
    if strcmp(Q.arm{I},'right')
        Q.arm_length_kb(I) = Q.arm_length_kb(I) - Q.arm_length_kb( strcmp(Q.arm,'left') & Q.chr_num==Q.chr_num(I)) ; 
    end
end

%% scatter plot and correlation to make sure I have the same results as Alsu

figure; 
plot(Q.arm_length_kb,Q.mean_length_underreplicated_left_mean_length_underreplicated_ri,'ok')
xlabel('chr arm length')
ylabel('% underrep')
title(corr( Q.arm_length_kb,Q.mean_length_underreplicated_left_mean_length_underreplicated_ri))

%% how does correlation change as we vary the threshold of minimum chromosome arm length ?

figure; 
hold on; 
for I = 1:300
    idx = Q.arm_length_kb > I & Q.mean_length_underreplicated_left_mean_length_underreplicated_ri < 999  ; 
    X = Q.mean_length_underreplicated_left_mean_length_underreplicated_ri( idx ) ; 
    Y = Q.arm_length_kb( idx ) ; 
    [c,p] = corr(X,Y);
    if p<0.05
        plot(I,c,'ok','MarkerFaceColor',[.7 .7 .7]);
    else
        plot(I,c,'.k');
    end
end
ylabel('Correlation')
xlabel('Threshold (min chr arm length, kb)')

%% figures : histograms and cdf of bootstrapps

X = Q.mean_length_underreplicated_left_mean_length_underreplicated_ri ; 
Y = Q.arm_length_kb ; 
[bootstat,bootsam] = bootstrp(5e4,@corr,X,Y);

fh = figure('units','centimeters','position',[5 5 10 7]) ;
yyaxis left ; 
histogram( bootstat , 50)
xlim([-1 0])
ylabel('# of bootstraps')

yyaxis right; 
[f,x] = ecdf(bootstat); 
plot(x,f,'LineWidth',2)
set(gca,'xtick',-1:0.1:1)
set(gca,'ytick',0:0.1:1)
ax = gca  ; 
ylabel('Fraction of bootstraps')
title('all arms ; 95% have corr > 0.4')
ax.YGrid = 'on' ; 
ax.XGrid = 'on' ; 
xlim([-1 0])
print('-dpng' , [FIGNAME '_all'] ,'-r300');
close ; 

clear 'X' 'Y' 'bootstat' ; 
%
idx = Q.arm_length_kb > 100 & Q.arm_length_kb < 800 & Q.mean_length_underreplicated_left_mean_length_underreplicated_ri < 65 ;

X = Q.mean_length_underreplicated_left_mean_length_underreplicated_ri( idx ) ; 
Y = Q.arm_length_kb( idx ) ; 
[bootstat,bootsam] = bootstrp(5e4,@corr,X,Y);

fh = figure('units','centimeters','position',[5 5 10 7]) ;
yyaxis left ; 
histogram( bootstat , 50)
xlim([-1 0])
ylabel('# of bootstraps')

yyaxis right; 
[f,x] = ecdf(bootstat); 
plot(x,f,'LineWidth',2)
set(gca,'xtick',-1:0.1:1)
set(gca,'ytick',0:0.1:1)
ax = gca  ; 
ylabel('Fraction of bootstraps')
ax.YGrid = 'on' ; 
ax.XGrid = 'on' ; 
xlim([-1 0])
title('arm length > 100kb & < 800kb , & underrep < 65%')
print('-dpng' , [FIGNAME '_subset'] ,'-r300');
close ; 