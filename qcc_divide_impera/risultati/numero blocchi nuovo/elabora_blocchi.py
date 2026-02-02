import pandas as pd
df=pd.read_csv("dim_blocco_sw_1.csv")
sw=df.groupby(["instance_name","num_chunks"]).num_swaps.mean().reset_index()
sw0=sw.loc[sw.groupby("instance_name")["num_swaps"].idxmin()]
tu=pd.read_csv("../../../revlib/tutti2.csv")
tab_sw=sw0.merge(tu)
tab_sw.sort_values("size",inplace=True)
tab_sw['dim_chunk']=tab_sw.depth / tab_sw.num_chunks
#tab=pd.pivot_table(res, values="num_swaps", index=["instance_name","depth"], columns="num_chunks", aggfunc="mean")
tab_sw.to_csv("opt_dim_sw.csv")
#print("Numero di chunk per gli swap")
#print(tab_sw)

df=pd.read_csv("dim_blocco_de_1.csv")
de=df.groupby(["instance_name","num_chunks"]).depth_out.mean().reset_index()
de0=de.loc[de.groupby("instance_name")["depth_out"].idxmin()]
tab_de=de0.merge(tu)
tab_de.sort_values("size",inplace=True)
tab_de['dim_chunk']=tab_de.depth / tab_de.num_chunks
#tab=pd.pivot_table(res, values="depth_out", index=["instance_name","depth"], columns="num_chunks", aggfunc="mean")
#print("Numero di chunk per la depth")
#print(tab_de)
tab_de.to_csv("opt_dim_de.csv")
