% %%
% M-M:
% 9, 10, 20, 29, 30
% M-G1:
% 1, 2, 5, 6, 21, 25
% M-G1 noBrdU:
% 18, 27, 32
% G1-HU:
% 19, 26, 31
% %%

%% load data
fn = dir('*tab');
fnG1 = {'1_S1_mapped.coverage_features.tab' '2_S2_mapped.coverage_features.tab'...
    '5_S3_mapped.coverage_features.tab'  '6_S4_mapped.coverage_features.tab' ...
    '21_S10_mapped.coverage_features.tab' '25_S11_mapped.coverage_features.tab'};
fnM = {'9_S5_mapped.coverage_features.tab' '10_S6_mapped.coverage_features.tab'...
    '20_S9_mapped.coverage_features.tab' '29_S14_mapped.coverage_features.tab' ...
    '30_S15_mapped.coverage_features.tab'};
fn = [fnG1 fnM];
clrs = [ repmat('r',1,numel(fnG1)) repmat('k',1,numel(fnM))] ;
lgt = [ repmat('G',1,numel(fnG1)) repmat('M',1,numel(fnM))] ;

T = table();
vn = {'chr' 'midpt' 'Exp' 'cnts_div_modefit' 'cnts_dmf_smooth' 'cnts_dmed_smooth'  'fn' } ;
for fnI = 1:numel(fn)
    Q = readtable(fn{fnI},'FileType','text');
    Q.Properties.VariableNames{1} = 'chr' ; 
    Q.Properties.VariableNames{4} = 'counts' ; 
    Q = Q(regexpcmp(Q.chr,'^chr'),:);
    Q.chr = categorical(Q.chr);
    Q.fn = repmat(fn(fnI),height(Q),1);
    Q.Exp = repmat({lgt(fnI)},height(Q),1);
    
    % mode / median fit
    G = grpstats(Q,'chr',{@modefit 'median'},'DataVars','counts') 
    for chrI = 1:height(G)
        idx = (Q.chr == G.chr(chrI)) ; 
        Q.cnts_div_modefit(idx) = Q.counts(idx) ./ G.modefit_counts(chrI) ; 
        Q.cnts_div_median(idx) = Q.counts(idx) ./ G.median_counts(chrI) ; 
    end
    Q.cnts_dmf_smooth = movmedian(Q.cnts_div_modefit , 5);
    Q.cnts_dmed_smooth = movmedian(Q.cnts_div_median , 5);
    Q.midpt = Q.Var2 + 100 ; 
    if G.modefit_counts(1) < 20 | G.median_counts(1) < 20
        warning( [ 'removing ' fn{fnI} ] )
    else
        T = vertcat( T , Q(:,vn) );
    end
end
T.Exp = categorical(T.Exp);
T.fn = categorical(T.fn);
T.midpt_kb = round(T.midpt ./ 1000) ; 



%% calculate
G = grpstats( T , {'chr' 'midpt_kb' 'Exp'}  , 'median' , 'DataVars' , {'cnts_dmf_smooth' 'cnts_div_modefit'} )  ; 

G1 = G( G.Exp=='G' , :);
M = G( G.Exp=='M' , :);
G1.Ratio_mf = G1.median_cnts_dmf_smooth ./ M.median_cnts_dmf_smooth  ;
G1.Ratio_med = G1.median_cnts_div_modefit./ M.median_cnts_div_modefit  ;

%%
figure; 
boxplot(G1.Ratio_mf(G1.midpt_kb < 25) , G1.midpt_kb(G1.midpt_kb < 25) ) ; grid on ; xlabel('dist from left tel')
ylabel('Ratio (modefit)')

figure; 
boxplot(G1.Ratio_med(G1.midpt_kb < 25) , G1.midpt_kb(G1.midpt_kb < 25) ) ; grid on ; xlabel('dist from left tel')
ylabel('Ratio (median)')
