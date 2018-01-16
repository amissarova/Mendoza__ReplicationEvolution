%% load DEE data
DATASAVEDIR = '~/Develop/Mendoza__ReplicationEvolution/Data/' ; 
FIGSAVEDIR  = [DATASAVEDIR '../Figures/'] ; 
load([DATASAVEDIR 'DEE.mat'] );
%% relationship between expression & % of the time some signal is lost
% confirm that auxotrophic markers are frequently 'lost'
%

X = DEE.expr_l2_nonzero ;   %  expression when reads is not zero
Y = DEE.fraction_experiments_zero_reads * 100+0.1 ;   % vs fraction of experiments w/zero reads

[xData, yData] = prepareCurveData( X , Y );
% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

Ypred = feval(fitresult ,X );
DEE.Los_minus_Losspred   =  Y - Ypred ;  

figure; hold on ;
plot(X,Y,'.k','DisplayName','all genes')
ylabel('% experiments w/zero read counts')
xlabel('log2( median read counts )')
set(gca,'xtick',0:100)
grid on
set(gca,'yscale','log')
set(gca,'ytick',     [0.1 0.2 0.5 1 2 5 10 20 50 100])
set(gca,'yticklabel',[0   0.2 0.5 1 2 5 10 20 50 100])

genes = {'TRP1' 'LEU2' 'URA3' 'HIS3'  'MET17' 'LYS2' 'ADE2'};
for I = 1:numel(genes)
    idx = strcmp(DEE.GENE_SGD,genes{I});
    plot( X(idx) , Y(idx),'.','DisplayName',genes{I},'MarkerSize',25)
    fprintf('%0.02f %0.02f\t%s\n' , X(idx) , Y(idx),genes{I});
end

%idx = (T.dist_to_end < 1e4);
%plot( X(idx) , Y(idx),'.c','DisplayName','<10kb','MarkerSize',10)

legend('location','sw')

xl = linspace(min(X),8,1e4);
plot( xl , feval(fitresult,xl) , '-r','DisplayName','fit')
ylim([0.1 60])
xlim([0.01 12])

% display residuals -- should be high for auxotrophic markers
DEE(ismember(DEE.GENE_T,genes) ,{'GENE_T','Los_minus_Losspred'});

%% save a minimal .mat file for joining w/other datasets
DEE = DEE(  : , {'ORF' 'GENE_T' 'fraction_experiments_zero_reads' 'expr_l2'   'Los_minus_Losspred'  'Chr' 'Start' 'End' 'dist_to_end' }) ; 
save( [DATASAVEDIR '/DEE_minimal.mat' ]  , 'DEE'); 
%% boxplot loss corrected for expression
idx = DEE.dist_to_end < 50*1000  ;
G = round( DEE.dist_to_end(idx) ./ 5000).*5000 ;
ug = unique(G) ;
figure;  hold on;
bh = boxplot( DEE.Los_minus_Losspred(idx) , G  ,'notch','on','plotstyle','compact','Symbol',''...
    ,'Positions',1:numel(ug) , 'Labels',ug/1000 ,'LabelOrientation','horizontal')
ylabel('gene loss residual (by expr)')
xlabel('distance from chr end (kb)')
line(xlim,[0 0],'Color',[.7 .7 .7])
ylim([-10 10])
for I = 1:numel(bh)
    set(bh(I),'Color',[.7 .7 .7])
end

%% loss not corrected for expression
idx = DEE.expr_l2 > log2(10) % 10 reads ;
G = round( DEE.dist_to_end ./ 5000).*5000 ;
ug = unique(G(idx)) ;
fh = figure('units','centimeters','position',[5 5 7 5 ]);
hold on;
Ypred = feval(fitresult ,X );
bh = boxplot( Y(idx)  , G(idx)  ,'notch','off','plotstyle','compact','Symbol',''...
    ,'Positions',1:numel(ug) , 'Labels',ug/1000 ,'LabelOrientation','horizontal' , 'Whisker',0 ) ; 
ylabel('% genes w/0 expression') ; 
xlabel('distance from chr end (kb)') ; 
line(xlim,[0 0],'Color',[.7 .7 .7])  ;

ylim([0 15]) ; 
xlim([0.5 10]) ; 
 
for I = 1:numel(bh)
    set(bh(I),'Color',[.7 .7 .7]) ; 
end
for I = 1:(numel(ug)-1)
    [ ~,p] = ttest2(  Y(idx & G==ug(I)) , Y(idx & DEE.dist_to_end>50000) , p );
    fprintf('%d\t%d\t%d\t%0.02f\t%0.02f\t' , I , ug(I) , sum(idx & G==ug(I)) , mean( Y(idx & G==ug(I))) , p );
    if p < 0.0001
%        fprintf('Y1\n');
        text( I-0.2 , max(ylim) , '***')
    elseif p < 0.001
%        fprintf('Y2\n');
        text( I-0.2 , max(ylim) , '**')
    else
%        fprintf('NO\n');
    end
end

%%
THRESHOLD_PRCTILE_EXPRESSION_TO_CONSIDER = 0 ; 
DEE.ChrRightArm = regexpcmp(DEE.ORF,'Y.R') ; 
DEE = sortrows(DEE , { 'Chr' 'ChrRightArm' 'dist_to_end'} ,'descend');
idx_genes_to_check = find(DEE.expr_l2 >= prctile(DEE.expr_l2, THRESHOLD_PRCTILE_EXPRESSION_TO_CONSIDER )) ; 
exprmat = table2array( DEE( : , regexpcmp(DEE.Properties.VariableNames , '^.RR\d')));

c = NaN(height(DEE),2);
for I = 1:height(DEE)
    if ismember(I,idx_genes_to_check)
    idx_exp_zero_expr = find(exprmat(I,:) == 0) ; % zero expression in these experiments
    idx_rest_of_chr =  find( DEE.Chr == DEE.Chr(I) & DEE.ChrRightArm == DEE.ChrRightArm(I) & DEE.dist_to_end < DEE.dist_to_end(I)) ; 
    rest_of_chr_has_zero_expr = all( exprmat( idx_rest_of_chr , idx_exp_zero_expr) == 0 , 1) ; 
    c(I,1) = sum(rest_of_chr_has_zero_expr) ; 
    c(I,2) = numel(idx_rest_of_chr) ; 
    c(I,3) = numel(idx_exp_zero_expr) ; 
    c(I,4) = DEE.dist_to_end(I) ;
    end
end


G = round( DEE.dist_to_end ./ 5000) .* 5 ;
idx = c(:,4) < 50000 ;

fh = figure('units','centimeters','position',[5 5 7 25 ]);


subplot(3,1,1)
boxplot( 100*(c(idx,1) ./ c(idx,3)) , c(idx,2) ) ; % by # of genes
xlim([0 10])
ylim([0 100])
xlabel('# genes from chr end')
ylabel('Conditional % time rest of arm is deleted')

subplot(3,1,2)
boxplot( 100*(c(idx,1) ./ c(idx,3)) , G(idx) ) ; % by dist from end
xlim([0 10])
ylim([0 100])
xlabel('Kb from chr end')
ylabel('Conditional % time rest of arm is deleted')


subplot(3,1,3)
boxplot( 100*(c(idx,1) ./ ncols(exprmat)) ,  G(idx) ) ; % by dist from end & of all experiments
xlim([0 10])

ylabel('% time rest of arm is deleted')
xlabel('Kb from chr end')

%  scatter plots intead of boxplots
idx = c(:,4) < 50000  ; 
fh = figure('units','centimeters','position',[5 5 7 7 ]);
subplot(2,1,1)
plot( DEE.dist_to_end(idx)./1000 ,  100*(c(idx,1) ./ ncols(exprmat)) ,'ok','MarkerFaceColor',[.7 .7 .7]); %absolute
ylabel('absolute')
set(gca,'xscale','log')
set(gca,'xtick',[.5 1 2 5 10 25 50])
xlim([.1 50])
set(gca,'ytick',0:10:100)


subplot(2,1,2)
plot( DEE.dist_to_end(idx)./1000 ,  100*(c(idx,1) ./ c(idx,3)) ,'ok','MarkerFaceColor',[.7 .7 .7]);  %conditional 
ylabel('conditional')
set(gca,'xscale','log')
xlabel('kb to the chromosome end')
set(gca,'xtick',[.5 1 2 5 10 25 50])
set(gca,'ytick',0:20:100)
%ylabel({'% of strains with the' 'entire arm deleted'})
legend('off')   
axis tight;
xlim([.1 50])

