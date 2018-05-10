
%%
cd ~/Develop/Mendoza__ReplicationEvolution/CLB2
addpath(genpath('~/Develop/matlab'));
%% Read data from tubes

replicates = {'FCS__20180420','FCS__20180427' , 'FCS__20180505' , 'FCS__20180506'};
wells = struct();
for Z = 1:length(replicates)
    wells1 = FACSReadPlateDir( replicates{Z} ,'fcsprefix','');
    wells1 = arrayfun(@FACSAssignWellID,wells1);
	% common metadata
    [wells1.ExpDate] = deal(20180420) ;
    % extract data from filenames
    for I = 1:length(wells1)
        try
            metadata  = regexp( wells1(I).filename , '_' ,'split');
            wells1(I).ReplicateNum = Z;
            if strcmp(metadata{2} , 'pre')
                wells1(I).TimePoint = -1;
                wells1(I).TimeH = -3;
                wells1(I).Condition = 'Pregrowth';
                wells1(I).Strain = str2num(char(metadata{3}));
            else
                wells1(I).TimePoint = str2num(char(metadata{2}));
                if Z == 1
                    wells1(I).TimeH = 0.75*wells1(I).TimePoint;
                    wells1(I).TimeMin = 45*wells1(I).TimePoint;
                    wells1(I).Condition = char(strrep( char(metadata{4}) , '.fcs' , ''));
                    wells1(I).Strain = str2num(char(metadata{3}));
                else
                    wells1(I).TimeH = 0.33334*wells1(I).TimePoint;
                    wells1(I).TimeMin = 20*wells1(I).TimePoint;
                    wells1(I).Condition = char(strrep( char(metadata{4}) , '.fcs' , ''));
                    wells1(I).Strain = str2num(char(metadata{3}));
                end
            end
        catch
        end
    end
    if isempty(fieldnames(wells))
        wells = wells1;
    else
        wells = [wells , wells1];
    end
end
wells = arrayfun(@FACSflagCellDensity , wells);   
%%
% Gate by cells size 
Gates_Size = FACSDrawGateForEach(wells,'group','PlateName','max_total_cells',10000,'xvar','FSC_A','yvar','SSC_A');
wells = FACSApplyGates(wells,'PlateName',Gates_Size);
wells = arrayfun(@FACSRemoveFilteredCells,wells);

Gates_Size = FACSDrawGateForEach(wells,'group','PlateName','max_total_cells',10000,'xvar','FSC_A','yvar','FSC_H');
wells = FACSApplyGates(wells,'PlateName',Gates_Size);
wells = arrayfun(@FACSRemoveFilteredCells,wells);

%% 
wells = arrayfun(@(x)FACSGetStats(x,'flunames',{'FITC_A' 'FSC_A' 'SSC_A'}),wells);
DS = struct2ds(wells);
DS.log_data_GFP = cell(length(DS) , 1);
% trimming data only for SytoxGreen positive cells (since for cells that lost the staining we can't say to which cell cycle they belong)
DS.trimmed_data_log = cell(length(DS) , 1);
DS.trimmed_data = cell(length(DS) , 1);
DS.TimeMin = NaN(length(DS) , 1);
for I = 1:length(DS)
    data = DS.data_FITC_A{I};
    data = nanmax(data , 2);
    DS.trimmed_data{I} = data(data > 2^13 & data < 2^15);
    data = log2(data);
    DS.log_data_GFP{I} = data;
    DS.trimmed_data_log{I} = data(data > 13 & data < 15);
    if DS.ReplicateNum(I) == 1
        DS.TimeMin(I) = DS.TimePoint(I)*45;
    else
        DS.TimeMin(I) = DS.TimePoint(I)*20;
    end
end
save('~/Develop/Mendoza__ReplicationEvolution/Data/DS_CLB2.mat' , 'DS');

