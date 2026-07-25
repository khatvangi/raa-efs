import pandas as pd
import os

families = ['PF00047', 'PF00069', 'PF00089-RAA', 'PF01145', 'PF12796', 'PF00505']
results = []

for fam in families:
    for csv in [f'{fam}/usalign_results.csv', f'{fam}/usalign_EFS_C_results.csv']:
        if os.path.exists(csv):
            df = pd.read_csv(csv)
            if 'TM_score' in df.columns or 'TM-score' in df.columns:
                tm_col = 'TM_score' if 'TM_score' in df.columns else 'TM-score'
                n = len(df)
                pass_05 = (df[tm_col] >= 0.5).sum()
                pass_07 = (df[tm_col] >= 0.7).sum()
                mean_tm = df[tm_col].mean()
                results.append({
                    'family': fam,
                    'file': os.path.basename(csv),
                    'n_structures': n,
                    'mean_TM': round(mean_tm, 3),
                    'pass_0.5': f"{pass_05}/{n} ({100*pass_05/n:.1f}%)",
                    'pass_0.7': f"{pass_07}/{n} ({100*pass_07/n:.1f}%)"
                })

df_out = pd.DataFrame(results)
print(df_out.to_string(index=False))
df_out.to_csv('RAA_all_families_summary.csv', index=False)
print("\nSaved to RAA_all_families_summary.csv")
