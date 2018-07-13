
%% Analysis of early-replicating regions (defined by HU block): are do they have higher copy number (not CNR)? 
% If so, we have some early S-phase cells.  (data are in CareyLab/ExternalData/Data/Alvino07 )
% LBC July 2018



%% Origin firing timing from Alvino 2007
load ~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat
ORI = readtable('~/Google Drive File Stream/My Drive/CareyLab/ExternalData/Alvino07/alvino_suppl_data3.xls');
ORI.Tp10_ = strcmp(ORI.Tp10_,'Y');
ORI.Tp12_5_ = strcmp(ORI.Tp12_5_,'Y');
ORI.Tp17_5_ = strcmp(ORI.Tp17_5_,'Y');
ORI.Tp25_ = strcmp(ORI.Tp25_,'Y');
ORI.Tp15_ = strcmp(ORI.Tp15_,'Y');

%%
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
    posKB = ORI.Coordinate_kb_(I);
    idx_in_G =  G.chr_num == chr ; 
    idx_in_chr = find(round(G.start_point_kb{ find(idx_in_G,1)}) == round(posKB)) ;
    CN_G1_this_pos = cell2mat( cellfun( @(X)median(X(idx_in_chr)) , G.CN_G1(idx_in_G) ,'UniformOutput',false) );
    CN_G1_all_chr  = cell2mat( cellfun( @(X)median(X) , G.CN_G1(idx_in_G) ,'UniformOutput',false) );
 %  fprintf('%d\t%d\t%0.0f\t%0.04f\t%0.04f\n' , I , chr , posKB , mean(CN_G1_all_chr) , mean(CN_G1_this_pos));
  %  if ORI.Tp10_(I)
        xloc = find(mat(I,:),1) + random('normal',0,0.1,1);
        plot( xloc , (CN_G1_this_pos) ,'ok' ,'MarkerFaceColor',clrs(find(mat(I,:),1),:) );
   % elseif ORI.Tp25_(I)
   %     plot( mean(CN_G1_all_chr) , mean(CN_G1_this_pos) ,'ok'); % 
 %   end
    CN(I,:) = [mean(CN_G1_all_chr) , mean(CN_G1_this_pos)  ];
    if find(mat(I,:),1) == 1
        t10 = vertcat(t10,CN_G1_this_pos); 
    elseif find(mat(I,:),1) == 2
        t12 = vertcat(t12,CN_G1_this_pos); 
    elseif find(mat(I,:),1) == 3
        t15 = vertcat(t15,CN_G1_this_pos);
    elseif find(mat(I,:),1) == 4
        t17 = vertcat(t17,CN_G1_this_pos); 
    elseif find(mat(I,:),1) == 5
        t25 = vertcat(t25,CN_G1_this_pos);    
    end
end
ylabel('Copy number in G1 arrested cells')
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
print('-dpng','~/Downloads/CopyNumber_Origins_Alvino07_G1.png','-r600')
close ;

%
% in M phase arrested cells
fh = figure('units','centimeters','position',[5 5 12 8]);
hold on ;
CN = NaN( height(ORI) , 2);
for I = 1:height(ORI)
    chr = ORI.Chromosome(I);
    posKB = ORI.Coordinate_kb_(I);
    idx_in_G =  G.chr_num == chr ; 
    idx_in_chr = find(round(G.start_point_kb{ find(idx_in_G,1)}) == round(posKB)) ;
    CN_G1_this_pos = cell2mat( cellfun( @(X)median(X(idx_in_chr)) , G.CN_M(idx_in_G) ,'UniformOutput',false) );
    CN_G1_all_chr  = cell2mat( cellfun( @(X)median(X) , G.CN_M(idx_in_G) ,'UniformOutput',false) );
 %  fprintf('%d\t%d\t%0.0f\t%0.04f\t%0.04f\n' , I , chr , posKB , mean(CN_G1_all_chr) , mean(CN_G1_this_pos));
  %  if ORI.Tp10_(I)
        xloc = find(mat(I,:),1) + random('normal',0,0.1,1);
        plot( xloc , (CN_G1_this_pos) ,'ok' ,'MarkerFaceColor',clrs(find(mat(I,:),1),:) );
   % elseif ORI.Tp25_(I)
   %     plot( mean(CN_G1_all_chr) , mean(CN_G1_this_pos) ,'ok'); % 
 %   end
    CN(I,:) = [mean(CN_G1_all_chr) , mean(CN_G1_this_pos)  ];
    if find(mat(I,:),1) == 1
        t10 = vertcat(t10,CN_G1_this_pos); 
    elseif find(mat(I,:),1) == 2
        t12 = vertcat(t12,CN_G1_this_pos); 
    elseif find(mat(I,:),1) == 3
        t15 = vertcat(t15,CN_G1_this_pos);
    elseif find(mat(I,:),1) == 4
        t17 = vertcat(t17,CN_G1_this_pos); 
    elseif find(mat(I,:),1) == 5
        t25 = vertcat(t25,CN_G1_this_pos);    
    end
end
ylabel('Copy number in M arrested cells')
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
print('-dpng','~/Downloads/CopyNumber_Origins_Alvino07_M.png','-r600')
close ;

%% boxplot binning by  Trep. 
% load data
load ~/Develop/Mendoza__ReplicationEvolution/Data/MutantChr200_Michi.mat
load('~/Develop/Mendoza__ReplicationEvolution/Data/DS_stat__200bp_new.mat')

for I = 1:height(G)
    trep = DS.Trep_spline( DS.chr_num == G.chr_num(I))
    boxplot( G.CN_G1{I} , round(trep))
    
end


