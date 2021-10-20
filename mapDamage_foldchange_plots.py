#!/usr/bin/env python
# coding: utf-8

# In[1]:


ROOT = './leeHom_trimmed_seqtk_filtered_bwa_aligned_mapDamage_results/'


# In[2]:


from glob import glob
import pandas as  pd
import matplotlib.pyplot as plt


# In[3]:


threePrime_files = glob('{}/ancient_tenanus/by_bioSample/*/3pGtoA_freq.txt'.format(ROOT))
fivePrime_files  = glob('{}/ancient_tenanus/by_bioSample/*/5pCtoT_freq.txt'.format(ROOT))


# In[4]:


threePrime_data = [ pd.read_csv(fname,delimiter='\t', index_col='pos').rename(columns={'3pG>A':fname.split('/')[-2]}) for fname in threePrime_files]
fivePrime_data  = [ pd.read_csv(fname,delimiter='\t', index_col='pos').rename(columns={'5pC>T':fname.split('/')[-2]}) for fname in fivePrime_files]


# In[5]:


threePrime_table = pd.concat(threePrime_data, axis=1).transpose().rename(columns={col:'3pG{}'.format(col) for col in range(1,26)})
fivePrime_table  = pd.concat(fivePrime_data, axis=1).transpose().rename(columns={col:'5pC{}'.format(col) for col in range(1,26)})


# In[6]:


mapDamage_table = pd.concat([threePrime_table, fivePrime_table], axis=1)


# In[7]:


sorted_threePrime_table = threePrime_table.sort_values(by=['3pG1'],ascending=False).transpose()
sorted_fivePrime_table  = fivePrime_table.sort_values(by=['5pC1'],ascending=False).transpose()


# In[8]:


sorted_threePrime_table_stats = sorted_threePrime_table.describe()
sorted_fivePrime_table_stats  = sorted_fivePrime_table.describe()


# In[ ]:


fig, axes = plt.subplots(1, 1, figsize=(10, 5))

plot_data = (sorted_threePrime_table.loc['3pG1'] / sorted_threePrime_table_stats.loc['mean']).to_frame().sort_values(by=0,ascending=False)
above_threshold = len(plot_data[plot_data[0] > 2.0])
axes.scatter(plot_data.index.values, plot_data, label='{} > 2.0'.format(above_threshold), alpha=1.00)
axes.plot([0,37],[2.0,2.0])
axes.legend()
plt.xticks(rotation=90, ha='right')
fig.savefig('3pG_foldchange.pdf')


# In[ ]:


fig, axes = plt.subplots(1, 1, figsize=(10, 5))

plot_data = (sorted_fivePrime_table.loc['5pC1'] / sorted_fivePrime_table_stats.loc['mean']).to_frame().sort_values(by=0,ascending=False)
above_threshold = len(plot_data[plot_data[0] > 2.0])
axes.scatter(plot_data.index.values, plot_data, label='{} > 2.0'.format(above_threshold), alpha=1.00)
axes.plot([0,37],[2.0,2.0])
axes.legend()
plt.xticks(rotation=90, ha='right')
fig.savefig('5pC_foldchange.pdf')


# In[14]:
fig, axes = plt.subplots(1, 1, figsize=(10, 5))

threePrime_plot_data = (sorted_threePrime_table.loc['3pG1'] / sorted_threePrime_table_stats.loc['mean']).to_frame().sort_values(by=0,ascending=False)
above_threshold = len(threePrime_plot_data[threePrime_plot_data[0] > 2.0])
axes.scatter(threePrime_plot_data.index.values, threePrime_plot_data, label='{} > 2.0'.format(above_threshold), alpha=1.00)


fivePrime_plot_data = (sorted_fivePrime_table.loc['5pC1'] / sorted_fivePrime_table_stats.loc['mean']).to_frame().sort_values(by=0,ascending=False)
above_threshold = len(fivePrime_plot_data[fivePrime_plot_data[0] > 2.0])
axes.scatter(fivePrime_plot_data.index.values, fivePrime_plot_data, label='{} > 2.0'.format(above_threshold), alpha=1.00)


axes.plot([0,37],[2.0,2.0])
axes.legend()
plt.xticks(rotation=90, ha='right')
axes.set_ylabel('Fold change')
fig.savefig('3pG_5pC_foldchange.pdf')

