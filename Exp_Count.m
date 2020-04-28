
%%
Animal = "Both"; Set_Path;
Beto_msk = find(ExpRecord.Exp_collection == "ReducDimen_Evol" & ...
    contains(ExpRecord.expControlFN, "generate") & ...
    contains(ExpRecord.expControlFN, "Beto"));
channels = unique(ExpRecord.pref_chan(Beto_msk))';
IT_ch = channels(channels<=32);
V4_ch = channels(channels>=49);
V1_ch = channels(channels>=33 & channels<=48);
fprintf("Reduce Dimension Evolutions %s\n","Beto")
fprintf("\nIT channels (%d): ",length(IT_ch));fprintf("%d ",IT_ch)
fprintf("\nV4 channels (%d): ",length(V4_ch));fprintf("%d ",V4_ch)
fprintf("\nV1 channels (%d): ",length(V1_ch));fprintf("%d ",V1_ch)
fprintf("\n")
%%
Animal = "Both"; Set_Path;
Beto_msk = find(ExpRecord.Exp_collection == "Manifold" & ...
    contains(ExpRecord.expControlFN, "generate") & ...
    contains(ExpRecord.expControlFN, "Beto"));
unique(ExpRecord.pref_chan(Beto_msk))'
%%
Animal = "Both"; Set_Path;
Alfa_msk = find(ExpRecord.Exp_collection == "Manifold" & ...
    contains(ExpRecord.expControlFN, "generate") & ...
    contains(ExpRecord.expControlFN, "Alfa"));
channels = unique(ExpRecord.pref_chan(Alfa_msk))';
IT_ch = channels(channels<=32);
V4_ch = channels(channels>=49);
V1_ch = channels(channels>=33 & channels<=48);
fprintf("Manifold Experiments %s\n","Alfa")
fprintf("\nIT channels (%d): ",length(IT_ch));fprintf("%d ",IT_ch)
fprintf("\nV4 channels (%d): ",length(V4_ch));fprintf("%d ",V4_ch)
fprintf("\nV1 channels (%d): ",length(V1_ch));fprintf("%d ",V1_ch)
fprintf("\n")