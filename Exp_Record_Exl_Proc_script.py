import pandas as pd
import numpy as np
df = pd.read_excel(r"Alfa_Exp_Record_draft.xlsx")
df2 = pd.read_excel(r"ExpSpecTable_Augment.xlsx")
row_num = df.shape[0]
#%%
# df.drop(df.isna().all(axis=1))
df = df.dropna(axis=0, how='all') # drop lines full with nan
df = df.reset_index(drop=True)  # (index=range(df.shape[0]))  # make the index contiguous!
#%%
msk = ~ df.ephysFN.isna()  # 500 rows
# df.ephysFN[msk].str.contains("Alfa")  # 436 rows containing Alfa
ExpEphysNames = df.ephysFN[df.ephysFN.str.contains("Alfa")==True]
RowidEphs = ExpEphysNames.index
ExpBhv2Names = df.expControlFN[df.expControlFN.str.contains("ALfa|Alfa")==True]
RowidBhv = ExpEphysNames.index
assert RowidEphs is RowidBhv
#%%
df.comments.fillna(value="", inplace=True)
df.comments = df.comments.astype("str")
df.stimuli.fillna(value="", inplace=True)
df.stimuli = df.stimuli.astype("str")
#%%
for Expi, rowi in enumerate(RowidEphs):
    if Expi != len(RowidEphs) - 1:
        nextrow = RowidEphs[Expi + 1]
    else:
        nextrow = row_num
    print("\nExp %d\t %s\t %s"%( Expi, df.ephysFN[rowi], df.expControlFN[rowi]))
    print(df.comments[rowi:nextrow].str.cat(sep="\n"))
#%%
stimuli_miss_cnt = 0
df_sort = df[df.ephysFN.str.contains("Alfa")==True]
df_sort = df_sort.reset_index(drop=True)
for Expi, rowi in enumerate(RowidEphs):
    if Expi != len(RowidEphs) - 1:
        nextrow = RowidEphs[Expi + 1]
    else:
        nextrow = row_num

    df_sort.comments[Expi] = df.comments[rowi:nextrow].str.cat(sep="\n")
    df_sort.ephysFN[Expi] = df.ephysFN[rowi]
    df_sort.expControlFN[Expi] = df.expControlFN[rowi]
    if "Stimuli" in df.stimuli[rowi]:
        df_sort.stimuli[Expi] = df.stimuli[rowi]
    else:
        df_sort.stimuli[Expi] = ""
        stimuli_miss_cnt += 1
        # print out info for further examination
        print("\nExp %d\t %s\t %s" % (Expi, df.ephysFN[rowi], df.expControlFN[rowi]))
        print(df.stimuli[rowi:nextrow].str.cat(sep=""))
        if ("Abort" in df_sort.comments[Expi]) or ("abort" in df_sort.comments[Expi]):
            print("Do aborted! No worry.")
print(stimuli_miss_cnt, "stimuli missing")
#%%
df_sort.to_excel("Alfa_Exp_Record_sort.xlsx")