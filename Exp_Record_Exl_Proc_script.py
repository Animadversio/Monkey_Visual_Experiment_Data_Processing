import pandas as pd
import numpy as np
df = pd.read_excel(r"Alfa_Exp_Record_draft.xlsx")
df2 = pd.read_excel(r"ExpSpecTable_Augment.xlsx")
row_num = df.shape[0]
#%%
# df.drop(df.isna().all(axis=1))
df = df.dropna(axis=0, how='all')
df = df.reindex(index=range(df.shape[0]))
#%%
msk = ~ df.ephysFN.isna() # 500 rows
df.ephysFN[msk].str.contains("Alfa")  # 436 rows containing Alfa
ExpEphysNames = df.ephysFN[df.ephysFN.str.contains("Alfa")==True]
RowidEphs = ExpEphysNames.index
ExpBhv2Names = df.expControlFN[df.expControlFN.str.contains("ALfa|Alfa")==True]
RowidBhv = ExpEphysNames.index
#%%
df.comments.fillna(value="", inplace=True)
df.comments = df.comments.astype("str")
#%%
for Expi, rowi in enumerate(RowidEphs):
    if Expi != len(RowidEphs) - 1:
        nextrow = RowidEphs[Expi + 1]
    else:
        nextrow = row_num
    print("\nExp %d\t %s\t %s"%( Expi, df.ephysFN[rowi], df.expControlFN[rowi]))
    print(df.comments[rowi:nextrow].str.cat(sep="\n"))