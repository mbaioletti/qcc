import pandas as pd
div=pd.read_csv("esp_swap_bis.csv")
#div_swaps=df.groupby("instance_name").agg(div_swaps=('num_swaps','min'),div_depth=('depth_out','min')).reset_index()
div_swaps=div.groupby(['instance_name','timeout']).agg(div_swaps=('num_swaps','min')).reset_index()
sab=pd.read_csv("/home/marco/OneDrive/ricerca/qcc/cpp/risultati/results_sabre.csv")
sab['num_swaps']=sab['best_size']-sab['g_in']
sab_swaps=sab.groupby(['instance_name','timeout']).agg(sab_swaps=('num_swaps','min'), size=('g_in','min')).reset_index()
confronto=sab_swaps.merge(div_swaps)
confronto=confronto.sort_values("size")
for t in [10,20,30,60]:
    c=confronto[confronto.timeout==t]
    nvitt=len(c[c.div_swaps< c.sab_swaps])
    npare=len(c[c.div_swaps==c.sab_swaps])
    nscon=len(c[c.div_swaps> c.sab_swaps])
    print(t,nvitt,npare,nscon)
