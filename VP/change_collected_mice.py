import pickle
import numpy as np
import pandas as pd

df = pd.read_csv('collected.csv')

host_strain_correction = {'BALB/CByJ': 'BALB/C',
						 'C57BL/6 MX1++': 'C57BL/6',
						 'DBA.2': 'DBA/2J',
						 'DBA/2': 'DBA/2J',
						 'C57BL/6 IFNAR-/-': np.nan,
						 'ICR': np.nan,
						 'SJL/JOrlCrl': np.nan}

df['Host_strain'].replace(host_strain_correction, inplace=True)
df.dropna(inplace=True)

df.to_csv('changed.csv', sep=',')