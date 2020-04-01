%%
Animal = "Alfa"; Set_Path;
expftr = contains(ExpRecord.expControlFN,"generate") & ...
        contains(ExpRecord.Exp_collection, "SUHash");
row_idx = find(expftr);
Expi = 1;
for i = 1:length(row_idx)
    if isempty(ExpRecord.stimuli{row_idx(i)}) || ...
        contains(ExpRecord.comments{row_idx(i)}, 'abort')
        ExpRecord(row_idx(i),:)
        keyboard
        continue
    end
    ExpRecord.Expi(row_idx(i)) = Expi; 
    Expi = Expi + 1;
    stim_path = ExpRecord.stimuli{row_idx(i)};
    if ~ contains(stim_path,'N:\') && contains(stim_path(1:7),'Stimuli')
        ExpRecord.stimuli{row_idx(i)} = fullfile('N:\', stim_path);
    end
end
%%
writetable(ExpRecord, 'S:\Exp_Record_Alfa.xlsx')