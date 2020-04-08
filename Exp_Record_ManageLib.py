import pandas as pd
import numpy as np
import os
if os.environ['COMPUTERNAME'] == 'DESKTOP-9DDE2RH':  # PonceLab-Desktop 3
    tmp_input = r"D:\ExpRecord_tmp.xlsx"
    tmp_output = r"D:\ExpRecord_out.xlsx"
    df_paths = [r"S:\Exp_Record_Alfa.xlsx", r"S:\ExpSpecTable_Augment.xlsx"]
elif os.environ['COMPUTERNAME'] == 'DESKTOP-MENSD6S':  # Home_WorkStation
    tmp_input = "E:\\Monkey_Data\\ExpRecord_tmp.xlsx"
    tmp_output = "E:\\Monkey_Data\\ExpRecord_out.xlsx"
    df_paths = ["E:\\Monkey_Data\\Exp_Record_Alfa.xlsx", "E:\\Monkey_Data\\ExpSpecTable_Augment.xlsx"]


def process_concat_cells(df, out_excel, Animal):
    """Process the raw form excel copied from onenote to well formed excel
    Filter the array using `Animal` label""
    """
    if isinstance(df,str):
        df = pd.read_excel(df)
    df = df.dropna(axis=0, how='all') # drop lines full with nan
    df = df.reset_index(drop=True)  # (index=range(df.shape[0]))  # make the index contiguous!
    row_num = df.shape[0]
    #%%
    msk = ~ df.ephysFN.isna()  # 500 rows
    # df.ephysFN[msk].str.contains("Alfa")  # 436 rows containing Alfa
    if Animal is "Alfa":
        search_str = "Alfa|ALfa"
    elif Animal is "Beto":
        search_str = "Beto"
    elif Animal is "Both":
        search_str = "Beto|Alfa|ALfa"
    else:
        search_str = "Beto|Alfa|ALfa"
    ExpEphysNames = df.ephysFN[df.ephysFN.str.contains(search_str)==True]
    RowidEphs = ExpEphysNames.index
    ExpBhv2Names = df.expControlFN[df.expControlFN.str.contains(search_str)==True]
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
    # 
    stimuli_miss_cnt = 0  # Count how many stimuli entries are missed
    df_sort = df[df.ephysFN.str.contains(search_str)==True]
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
    df_sort.to_excel(out_excel,index=False, engine='xlsxwriter')
    # Use the 'xlsxwriter' engine will avoid some Illegal Character Error in openpyxl
    # https://cooperluan.github.io/python/2015/01/08/pandas-daochu-excel-zhong-de-luan-ma-wen-ti/
    return df_sort

def concat_table(df_old, df_new, addexplabel=None, out_path=None):
    """Obsolete, use the function below instead"""
    if isinstance(df_old,str):
        out_path = df_old
        df_old = pd.read_excel(df_old)
    if isinstance(df_new,str):
        df_new = pd.read_excel(df_new)
    # Check if the experiments in the new datatable has been recorded
    break_flag = False
    for name in df_new.expControlFN:
        if (df_old.expControlFN==name).any():
            print("%s  has been recorded in the excel index %d, please check"%(name, (df_old.expControlFN==name).nonzero()[0][0]))
            break_flag = True
    if break_flag:
        raise ValueError
    if addexplabel is not None:
        df_new.Exp_collection[:] = addexplabel
    df_cat = pd.concat([df_old,df_new], axis=0, ignore_index=True)
    df_cat.to_excel(out_path,index=False)
    return df_cat

def sort_merge_table(df_sort, addexplabel=None):
    Animal_strs = ["Alfa", "Beto"]
    if isinstance(df_sort,str):
        df_sort = pd.read_excel(df_sort)
    # loop through animal name and sort corresponding exp to the collection
    for animal, out_path in zip(Animal_strs, df_paths):
        print("Sort out exp for %s, adding "%animal)
        df_old = pd.read_excel(out_path) # load the old exp collection
        id_col = []
        for idx in df_sort.index:
            name = df_sort.expControlFN[idx]
            if name is np.nan:
                print("%s Empty bhv file entry encountered" % df_sort.ephysFN[idx])
                continue
            if animal in name:
                if (df_old.expControlFN==name).any():
                    print("%s  has been recorded in the excel index %d, please check. Skipping."%(name, (df_old.expControlFN==name).nonzero()[0][0]))
                else:
                    id_col.append(idx)
        if len(id_col) == 0:
            print("\nNo new experiments to add! Continue.")
            continue
        df_ftr = df_sort.iloc[id_col].copy()
        print(df_ftr.expControlFN)
        if addexplabel is not None:
            df_ftr.Exp_collection[:] = addexplabel
        df_cat = pd.concat([df_old, df_ftr], axis=0, ignore_index=True)
        df_cat.to_excel(out_path,index=False, engine='xlsxwriter') # write to the excel of all old experiments 
    return

if __name__ == '__main__':
    os.startfile(tmp_input)
    Animal = input("Which animal to parse from %s?" % (tmp_input) )
    if len(Animal) == 0:
        Animal = "Both"#"Beto" # "Alfa" "ALfa"
    df_sort = process_concat_cells(tmp_input, tmp_output, Animal=Animal)
    try:
        os.startfile(tmp_output)
    except:
        print("Open %s failed"%tmp_output)
        pass
    Label = input("Add Exp labels to the new Exps?")
    if len(Label) == 0:
        Label = None # "ReducDimen_Evol"
    sort_merge_table(tmp_output, addexplabel=Label) # df_sort