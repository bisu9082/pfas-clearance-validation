import sys, pandas as pd, numpy as np
d=pd.read_csv('cr/ecotox_pfas_measured_recurated.csv').rename(columns={'parent':'canonical_smiles'})
src=open('scripts/oecd_classify_v2.py').read().split('d = pd.read_csv(SRC)')[0]
ns={'__name__':'oc'}; exec(compile(src,'oc','exec'),ns); classify=ns['classify']
out=pd.concat([d.reset_index(drop=True),pd.DataFrame([classify(s) for s in d.canonical_smiles])],axis=1)
out.to_csv('cr/oecd_103set.csv',index=False)
t=out[out.oecd_pfas].groupby('oecd_class').agg(n=('logTox','size'),
   n_acute=('logTox_median',lambda s:s.notna().sum()), median_min=('logTox','median'),
   median_chain=('pf_chain','median'))
print(t.sort_values('n',ascending=False).to_string(float_format=lambda x:f'{x:6.2f}'))
print(f'\nsum {t.n.sum()} / OECD PFAS {out.oecd_pfas.sum()} / total {len(out)}')
print('\nshort-chain matches reassigned away from 400:')
for _,r in out[(out.pf_chain<4)&(out.oecd_pfas)].iterrows():
    print(f'   chain={int(r.pf_chain)} -> {r.oecd_class:<36} {str(r["name"])[:46]}')
