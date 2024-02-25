import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from pyscenic.export import export2loom, add_scenic_metadata
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns
import scanpy as sc
import anndata
import re
from imp import reload
from distributed import Client, LocalCluster
#import umap

#from geosketch import gs
#import pynndescent
import matplotlib.pyplot as pl
#import PyIOH5 as myh5
import matplotlib.pyplot as plt
ncores = 8
nthreads = 8

def _derive_threshold(auc_mtx: pd.DataFrame, regulon_name: str) -> float:
    assert auc_mtx is not None and not auc_mtx.empty
    assert regulon_name in auc_mtx.columns
    # Fit a two component Gaussian Mixture model on the AUC distribution using an Expectation-Maximization algorithm.
    data = auc_mtx[regulon_name].values.reshape(-1, 1)
    gmm = mixture.GaussianMixture(n_components=2, covariance_type='full').fit(data)
    avgs = gmm.means_
    stds = np.sqrt(gmm.covariances_.reshape(-1, 1))

    # The threshold is based on the distribution with the highest mean and is defined as (mu - 2 x std)
    idx = np.argmax(avgs)
    threshold = max(avgs[idx] - 2 * stds[idx], 0)
    # This threshold cannot be lower than (mu + 2 x std) based on the distribution with the lowest mean.
    idx = np.argmin(avgs)
    lower_bound = avgs[idx] + 2 * stds[idx]

    return max(lower_bound, threshold)

def binarize(auc_mtx: pd.DataFrame, threshold_overides:Optional[Mapping[str,float]]=None) -> (pd.DataFrame, pd.Series):
    """
    "Binarize" the supplied AUC matrix, i.e. decide if for each cells in the matrix a regulon is active or not based
    on the bimodal distribution of the AUC values for that regulon.
    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param threshold_overides: A dictionary that maps name of regulons to manually set thresholds.
    :return: A "binarized" dataframe and a series containing the AUC threshold used for each regulon.
    """
    def derive_thresholds(auc_mtx):
        return pd.Series(index=auc_mtx.columns, data=[_derive_threshold(auc_mtx, name) for name in auc_mtx.columns])
    thresholds = derive_thresholds(auc_mtx)
    if threshold_overides is not None:
        thresholds[list(threshold_overides.keys())] = list(threshold_overides.values())
    return (auc_mtx > thresholds).astype(int), thresholds

def plot_binarization(auc_mtx: pd.DataFrame, regulon_name: str, bins: int=200, threshold=None, ax=None) -> None:
    """
    Plot the "binarization" process for the given regulon.
    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param regulon_name: The name of the regulon.
    :param bins: The number of bins to use in the AUC histogram.
    :param threshold: The threshold to use for binarization. If None then this will be derived automatically.
    """
    if ax is None:
        ax=plt.gca()
    auc_mtx[regulon_name].hist(bins=bins,ax=ax)
    if threshold is None:
        threshold = _derive_threshold(auc_mtx, regulon_name)

    ylim = ax.get_ylim()
    ax.plot([threshold]*2, ylim, 'r:')
    ax.set_ylim(ylim)
    ax.set_xlabel('AUC')
    ax.set_ylabel('#')
    ax.set_title(regulon_name)



RESOURCES_FOLDER="pyscenic/resources/"
DATABASE_FOLDER = "pyscenic/databases/"
#SCHEDULER="113.105.131.192:8176"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm10_*.mc9nr.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_mgi_tfs.txt')
#SC_EXP_FNAME = os.path.join(DATA_FOLDER, "WT.hvg8k.norm.data.csv")
REGULONS_FNAME = "lee_astr/regulons.p"
MOTIFS_FNAME = "lee_astr/motifs.csv"
REGULONS_DF_FNAME = "lee_astr/regulons.csv"
AUCMTX_FNAME = "lee_astr/auc_mtx.csv"

LOOM_FILE =  "lee_astr/pyscenic.auc.loom"
adata_fname = "lee_astr/auc.scanpy.h5ad"

###load single cell astrocyte normalized expression matrix
ex_matrix = pd.read_csv("L2.Astrocyte.nC2500.rmMt_Rp.SCTregICMR.doscale.pc16.res035.k40.SCT.normalized.data.csv",index_col=0).T
ex_matrix.shape
#load all tf
tf_names = load_tf_names(MM_TFS_FNAME)
#retain genes both at tf database and our hvg3000
tf_names = [i for i in tf_names if i in ex_matrix.columns.values]
len(tf_names)

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True)

# Save adjacencies matrix
adjacencies.to_csv("lee_astr/adjacencies.csv")
adjacencies.shape

modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

modules=np.array(modules)
np.save("lee_astr/modules.npy",modules)

modules=np.load("lee_astr/modules.npy",allow_pickle=True)
modules=modules.tolist()

# Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)

#with open(REGULONS_FNAME,"rb") as f:
#    regulons=pickle.load(f)
    
name=[];tfs=[];targets = [];score = [];state=[];motif=[]
for i in regulons:
    name.append(i.name)
    tfs.append(i.transcription_factor)
    targets.append(','.join(i.genes))
    score.append(i.score)
    ct=list(i.context)
    if 'png' in ct[0]:
        motif.append(ct[0].split('.')[0])
    elif 'png' in ct[1]:
        motif.append(ct[1].split('.')[0])
    else:
        motif.append('')

regulons_df = pd.DataFrame(data={'name':name,'tfs':tfs, 'score':score,'targets':targets,'motif':motif})
regulons_df.to_csv(REGULONS_DF_FNAME, index=False)

#Phase III: cellular regulon enrichment matrix (aka AUCell)
auc_mtx = aucell(ex_matrix, regulons, num_workers=4)

auc_mtx.to_csv(AUCMTX_FNAME)

motif=load_motifs(MOTIFS_FNAME)

regulons = [r.rename(r.name.replace('(+)',' ('+str(len(r))+'g)')) for r in regulons] 

export2loom(ex_matrix, regulons, 
                LOOM_FILE,
                #title = "Zeisel et al.",
                #nomenclature = "MGI"
           )
#RELOADING THE ENRICHED MOTIFS AND REGULONS FROM FILE SHOULD BE DONE AS FOLLOWS
#auc_mtx=pd.read_csv(AUCMTX_FNAME,index_col=0)

###