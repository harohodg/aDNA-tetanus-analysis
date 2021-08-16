#!/usr/bin/env python
# coding: utf-8

# In[1]:


ROOT = './tetanus_aDNA_analysis'
GFF_COLUMNS = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

import pandas as pd
import pysam
from glob import glob
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm

#from IPython.core.display import display, HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))


# In[2]:


plasmid_gff    = pd.read_csv('{}/E88_plasmid.gff3'.format(ROOT), header=None, skiprows=2, delimiter='\t', names=GFF_COLUMNS)
chromosome_gff = pd.read_csv('{}/E88_chromosome.gff3'.format(ROOT), header=None, skiprows=2, delimiter='\t', names=GFF_COLUMNS)


# In[3]:


meta_data           = pd.read_csv('./metadata.txt', delimiter='\t', index_col='Run')
ancient_DNA_samples = meta_data[ meta_data['aDNA'] == 'Y']


# In[4]:


plasmid_bam_files    = glob('{}/plasmid/singles/*.sorted.bam'.format(ROOT))
chromosome_bam_files = glob('{}/chromosome/singles/*.sorted.bam'.format(ROOT))


# In[5]:


filtered_plasmid_bam_files   = [ fname for fname in plasmid_bam_files if fname.split('/')[-1].split('_')[0] in ancient_DNA_samples.index.values ]
filtered_chromosome_bam_files = [ fname for fname in chromosome_bam_files if fname.split('/')[-1].split('_')[0] in ancient_DNA_samples.index.values ]


# In[6]:


merged_plasmid_files = glob('{}/plasmid/merged/*.sorted.bam'.format(ROOT))
merged_chromosome_files = glob('{}/chromosome/merged/*.sorted.bam'.format(ROOT))


# In[7]:


len(plasmid_bam_files), len(chromosome_bam_files), len(filtered_plasmid_bam_files), len(filtered_chromosome_bam_files), len(merged_plasmid_files), len(merged_chromosome_files)


# In[8]:


def calculate_coverage_stats(fname, target, region_start, region_end, coverage_type='total'):
    assert coverage_type in ['total', 'percent']
    
    coverage_data = np.zeros( region_end - region_start + 1)
    
    samfile = pysam.AlignmentFile(fname, "rb")#Currently assumed to be a bam file
    for pileupcolumn in samfile.pileup(target, region_start, region_end):
        if pileupcolumn.pos >= region_start and pileupcolumn.pos <= region_end:
            coverage_data[pileupcolumn.pos - region_start] = pileupcolumn.n if coverage_type =='total' else 1
    samfile.close()
    
    return coverage_data


# In[9]:


def get_coverage_stats(coverage_targets, coverage_type='total'):
    coverage_stats = {}
    for coverage_target in coverage_targets:
        coverage_stats[ coverage_target['label'] ] = {}
        for fname in tqdm(coverage_target['files']):
            fname_label   = fname.split('/')[-1].split('_')[0]
            coverage_stats[ coverage_target['label'] ][fname_label] = np.mean( calculate_coverage_stats(fname, coverage_target['target'], coverage_target['region_start'], coverage_target['region_end'], coverage_type) )
    return coverage_stats


# In[10]:


plasmid_seqid    = plasmid_gff['seqid'].values[0]
chromosome_seqid = chromosome_gff['seqid'].values[0]
coverage_targets = [{'label' : 'TeNT', 'target':plasmid_seqid, 'region_start':68640, 'region_end':72587, 'files':filtered_plasmid_bam_files},
                    {'label' : 'CoIT', 'target':plasmid_seqid, 'region_start':39438, 'region_end':42413, 'files':filtered_plasmid_bam_files},
                    {'label' : 'Plasmid', 'target':plasmid_seqid, 'region_start':1, 'region_end':74082, 'files':filtered_plasmid_bam_files},
                    {'label' : 'Chromosome', 'target':chromosome_seqid, 'region_start':1, 'region_end':2799251, 'files':filtered_chromosome_bam_files}
                   ]

merged_coverage_targets = [{'label' : 'TeNT', 'target':plasmid_seqid, 'region_start':68640, 'region_end':72587, 'files':merged_plasmid_files},
                    {'label' : 'CoIT', 'target':plasmid_seqid, 'region_start':39438, 'region_end':42413, 'files':merged_plasmid_files},
                    {'label' : 'Plasmid', 'target':plasmid_seqid, 'region_start':1, 'region_end':74082, 'files':merged_plasmid_files},
                    {'label' : 'Chromosome', 'target':chromosome_seqid, 'region_start':1, 'region_end':2799251, 'files':merged_chromosome_files}
                   ]


# In[11]:


total_coverage_stats        = pd.DataFrame( get_coverage_stats(coverage_targets, 'total') )
total_merged_coverage_stats = pd.DataFrame( get_coverage_stats(merged_coverage_targets, 'total') )

percent_coverage_stats        = pd.DataFrame( get_coverage_stats(coverage_targets, 'percent') )
percent_merged_coverage_stats = pd.DataFrame( get_coverage_stats(merged_coverage_targets, 'percent') )


# In[14]:


def plot_coverage(data_table, sorting_target, bar_width=1, low_to_high=False):
    data_table.sort_values(by=sorting_target, ascending=False, inplace=True)

    fig, axes = plt.subplots(len(data_table.columns),1, figsize=(20,10), sharex=True, sharey=True)
    for row_index, row_label in enumerate(data_table.columns.values):
        data_table[row_label].plot.bar(width=bar_width, ax=axes[row_index])
        axes[row_index].set_ylabel(row_label)
        
    return fig, axes


# In[18]:


fig, axes = plot_coverage(total_coverage_stats, 'TeNT')
fig.suptitle('Average Number of Reads per Base')
fig.savefig('aDNA_Tetanus_coverage.pdf')


# In[19]:


fig, axes = plot_coverage(total_merged_coverage_stats, 'TeNT')
fig.suptitle('Average Number of Reads per Base')
fig.savefig('aDNA_Tetanus_bioSamples_merged.pdf')


# In[23]:


fig, axes = plot_coverage(percent_coverage_stats*100, 'TeNT')
fig.suptitle('Percent Coverage')
fig.savefig('aDNA_Tetanus_percent_coverage.pdf')


# In[24]:


fig, axes = plot_coverage(percent_merged_coverage_stats*100, 'TeNT')
fig.suptitle('Percent Coverage')
fig.savefig('aDNA_Tetanus_bioSamples_merged_percent_coverage.pdf')


# In[ ]:





# In[ ]:


hatching_options = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']

fig, axes = plt.subplots(2,1, figsize=(20,5), sharex=True, sharey=True)

column_labels = coverage_stats.index.values
x = np.arange(len(column_labels))  # the label locations
width = 1 # the width of the bars


cm = plt.get_cmap('gist_rainbow')
num_colors = len( ancient_DNA_samples['BioProject'].unique() )
axes[0].set_prop_cycle(color=[cm(1.*i/num_colors) for i in range(num_colors)])
for bio_project in ancient_DNA_samples['BioProject'].unique():
    x_values = []
    sraIDS = ancient_DNA_samples[ ancient_DNA_samples['BioProject'] == bio_project].index
    for column_index, sraID in enumerate(column_labels):
        if sraID in sraIDS:
            x_values.append(x[column_index])
    y_values = [1]*len(x_values)
    axes[-2].bar(x_values, y_values, width, edgecolor='black', label=bio_project)


cm = plt.get_cmap('viridis')
num_colors = len( ancient_DNA_samples['BioSample'].unique() )
axes[1].set_prop_cycle(color=[cm(1.*i/num_colors) for i in range(num_colors)])
for bio_sample in ancient_DNA_samples['BioSample'].unique():
    x_values = []
    sraIDS = ancient_DNA_samples[ ancient_DNA_samples['BioSample'] == bio_sample].index
    for column_index, sraID in enumerate(column_labels):
        if sraID in sraIDS:
            x_values.append(x[column_index])
    y_values = [1]*len(x_values)
    axes[-1].bar(x_values, y_values, width, edgecolor='black', label=bio_sample)
    
for ax in axes:
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    
fig.savefig('aDNA_Tetanus_coverage_bioSamples_bioProject_labels.pdf')


# In[ ]:





# In[ ]:


hatching_options = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']

fig, axes = plt.subplots(1,1, figsize=(20,5), sharex=True, sharey=True)

column_labels = merged_coverage_stats.index.values
x = np.arange(len(column_labels))  # the label locations
width = 1 # the width of the bars


cm = plt.get_cmap('gist_rainbow')
num_colors = len( ancient_DNA_samples['BioProject'].unique() )
axes.set_prop_cycle(color=[cm(1.*i/num_colors) for i in range(num_colors)])
for bio_project in ancient_DNA_samples['BioProject'].unique():
    x_values = []
    bioSamples = ancient_DNA_samples[ ancient_DNA_samples['BioProject'] == bio_project]['BioSample'].values
    for column_index, label in enumerate(column_labels):
        if label in bioSamples:
            x_values.append(x[column_index])
    y_values = [1]*len(x_values)
    axes.bar(x_values, y_values, width, edgecolor='black', label=bio_project)


    

axes.set_yticklabels([])
axes.set_xticklabels([])
    

fig.savefig('aDNA_Tetanus_coverage_merged_bioProject_labels.pdf')


# In[27]:


total_coverage_stats.to_csv('total_coverage.csv', index_label='sraID')
total_merged_coverage_stats.to_csv('total_merged_coverage.csv', index_label='BioSample')

(percent_coverage_stats*100).to_csv('percent_coverage.csv', index_label='sraID')
(percent_merged_coverage_stats*100).to_csv('percent_merged_coverage.csv', index_label='BioSample')


# In[ ]:


merged_coverage_stats.to_csv('merged_biosamples_coverage.csv',index_label='BioSample')


# In[ ]:




