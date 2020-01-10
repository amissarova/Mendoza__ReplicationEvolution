
%% Origin firing timing from Alvino 2007
load ~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat
ORI = readtable('~/Google Drive File Stream/My Drive/CareyLab/ExternalData/Alvino07/alvino_suppl_data3.xls');
ORI.Tp10_ = strcmp(ORI.Tp10_,'Y');
ORI.Tp12_5_ = strcmp(ORI.Tp12_5_,'Y');
ORI.Tp17_5_ = strcmp(ORI.Tp17_5_,'Y');
ORI.Tp25_ = strcmp(ORI.Tp25_,'Y');
ORI.Tp15_ = strcmp(ORI.Tp15_,'Y');

% remove rDNA and regions deleted in our yeast strain
% N_begin = 450.663;
% N_end = 491;
% chr XII 
%  deleted region: chr III, 148-152kb



%% load data from Muller14
T = readtable('~/Google Drive File Stream/My Drive/CareyLab/ExternalData/Muller14/coverage.tab','FileType','text','Delimiter','\t');
T = T( regexpcmp(T.chr,'^chr'),:) ;

% SRR926346' # S
% SRR926345' # G2
% SRR926343' # exponential
% SRR926344' # stationary


%% calc normalized CN
G = grpstats( T ,'chr',@modefit,'DataVars',{'SRR926343' 'SRR926344' 'SRR926345' 'SRR926346'});
T.exponential = NaN(height(T),1);
T.stationary = NaN(height(T),1);
T.G2 = NaN(height(T),1);
T.S = NaN(height(T),1);
for I = 1:height(T)
    T.exponential(I) = T.SRR926343(I) / G.modefit_SRR926343( strcmp(G.chr,T.chr{I})) ; 
    T.stationary(I) = T.SRR926344(I) / G.modefit_SRR926344( strcmp(G.chr,T.chr{I})) ; 
    T.G2(I) = T.SRR926345(I) / G.modefit_SRR926345( strcmp(G.chr,T.chr{I})) ; 
    T.S(I) = T.SRR926346(I) / G.modefit_SRR926346( strcmp(G.chr,T.chr{I})) ; 
end


%%
figname = '~/Downloads/Muller14.eps' ;
delete(figname);
for chrI = 1:16
    figure; 
    hold on ; 
    h=gcf;
    set(h,'PaperOrientation','landscape');
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition', [0 0 1 0.4]);

    ws = 20 ;
    chr = ['chr' num2roman(chrI)];
    idx =  strcmp(T.chr,chr);
    plot( T.start(idx)./1000 , smooth(T.S(idx),ws),'.' ,'DisplayName','S')
    plot( T.start(idx)./1000 , smooth(T.exponential(idx),ws),'.' ,'DisplayName','exponential')
    plot( T.start(idx)./1000 , smooth(T.stationary(idx),ws),'.' ,'DisplayName','stationary')
    plot( T.start(idx)./1000 , smooth(T.G2(idx),ws),'.' ,'DisplayName','G2')
    legend('location','best')
    axis tight;
    ylim([0.5 2])
    set(gca,'xtick',0:25:max(xlim))
    title(chr);
    xlabel('Position along chr (kb)');
    ylabel('Copy Number');
    
    X = ORI.Coordinate_kb_(ORI.Tp25_ & ORI.Chromosome == chrI) ;
    ph = plot( X , 1 ,'sk','MarkerFaceColor',[.7 .7 .7],'MarkerSize',10);
    set(ph,'HandleVisibility','off')
    X = ORI.Coordinate_kb_(ORI.Tp10_ & ORI.Chromosome == chrI) ;
    print('-dpsc2',figname,'-append');
    close;
end
%ph = plot( X , 1 ,'ok','MarkerFaceColor','b') ; 
%set(ph,'HandleVisibility','off')

%% Ori CN in S phase arrested cells
%%
T.midpt_kb = round((T.start + 100)./1000);
G = grpstats( T , {'chr' 'midpt_kb'} , {'mean' 'median'} , 'DataVars' , {'S' 'exponential' 'stationary' 'G2'});
G.chr_num = cellfun( @roman2num ,regexprep(G.chr,'chr','')) ; 
G = G( ~isnan(G.chr_num),:);
Y = G.median_G2 ; % which CN data to use? 
mat = table2array(ORI(:,3:end)) ; 
clrs = parula(5); 
fh = figure('units','centimeters','position',[5 5 12 8]);
hold on ;
CN = NaN( height(ORI) , 2);
t10 = NaN(0);
t12 = NaN(0);
t15 = NaN(0);
t17 = NaN(0);
t25 = NaN(0);
for I = 1:height(ORI)
    chr = ORI.Chromosome(I);
    posKB = round( ORI.Coordinate_kb_(I)) ;
    idx_in_G =  G.midpt_kb==posKB & G.chr_num == chr ; 
    CN_this_pos = Y( idx_in_G );
    CN_all_chr  = median(Y( G.chr_num == chr ));
 %  fprintf('%d\t%d\t%0.0f\t%0.04f\t%0.04f\n' , I , chr , posKB , mean(CN_G1_all_chr) , mean(CN_G1_this_pos));
  %  if ORI.Tp10_(I)
        xloc = find(mat(I,:),1) + random('normal',0,0.1,1);
        plot( xloc , (CN_this_pos) ,'ok' ,'MarkerFaceColor',clrs(find(mat(I,:),1),:) );
   % elseif ORI.Tp25_(I)
   %     plot( mean(CN_G1_all_chr) , mean(CN_G1_this_pos) ,'ok'); % 
 %   end
    if find(mat(I,:),1) == 1
        t10 = vertcat(t10,CN_this_pos); 
    elseif find(mat(I,:),1) == 2
        t12 = vertcat(t12,CN_this_pos); 
    elseif find(mat(I,:),1) == 3
        t15 = vertcat(t15,CN_this_pos);
    elseif find(mat(I,:),1) == 4
        t17 = vertcat(t17,CN_this_pos); 
    elseif find(mat(I,:),1) == 5
        t25 = vertcat(t25,CN_this_pos);    
    end
end
ylabel('Copy number in G2 phase cells')
xlabel('Origin firing timing (min after alpha factor release)')
set(gca,'xtick',1:5)
xlim([0.5 5.5])
set(gca,'xticklabel',{'10''' '12.5''' '15''' '17.5''' '25'''})
line( xlim , [1 1] , 'LineStyle','--','LineWidth',2,'Color',[.7 .7 .7])
ylim([0.75 2])

line( [0.75 1.25] , [median(t10) median(t10)] ,'Color','r' , 'LineWidth',3)
line( [1.75 2.25] , [median(t12) median(t12)] ,'Color','r' , 'LineWidth',3)
line( [2.75 3.25] , [median(t15) median(t15)] ,'Color','r' , 'LineWidth',3)
line( [3.75 4.25] , [median(t17) median(t17)] ,'Color','r' , 'LineWidth',3)
line( [4.75 5.25] , [median(t25) median(t25)] ,'Color','r' , 'LineWidth',3)
%print('-dpng','~/Downloads/CopyNumber_Origins_Alvino07_G1.png','-r600')
%close ;

