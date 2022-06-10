import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS
from rdkit.ML.Cluster import Butina

import hdbscan
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.cluster import DBSCAN



def get_largest_fragment(mol):
    '''
    Source: https://gist.github.com/PatWalters/3bb9f1091d989665caf37460fec969f3
    A simple function to return largest fragment in a molecule. Useful for stripping counterions.
    
    Input
    -----
    mol : RDKit mol
        RDKit molecule object
    
    Output
    ------
    frags : RDKit mol
        largest fragment from input RDKit mol object
    '''
    frags = list(Chem.GetMolFrags(mol, asMols=True))
    frags.sort(key=lambda x: x.GetNumAtoms(), reverse=True)
    return frags[0]

def display_cluster_members(df, sel, align_mcs=True, strip_counterions=True, hdbscan=False):
    '''
    Source: https://gist.github.com/PatWalters/3bb9f1091d989665caf37460fec969f3
    A function to generate an image of molecules in a selected cluster.
    
    Input
    -----
    df : dataframe
        cluster information for selected cluster
    sel : selection
        selection when working with mol2grid graphical interface
    align_mcs : bool
        set whether to aligned based on maximum common substructure; does not work too well with ring systems
    strip_counterions : bool
        set whether to remove counterions; if True, calls `get_largest_fragment`
    hdbscan : bool
        indicate whether the results will be fed into HDBSCAN clustering
        if True, prints strength of cluster membership along with P(active) value when visualizing compounds in clusters
        
    Output
    ------
    img : Image
        returns visualization of compounds for the cluster selected via the mols2grid graphical interface
    '''
    mol_list = []
    img = 'Nothing selected'
    if len(sel):
        sel_df = df.query('Cluster in @sel')
        sel_df = sel_df.sort_values(by='P(active)', ascending=False)
        # order by ascending P(active)
        mol_list = [Chem.MolFromSmiles(smi) for smi in sel_df['stereo_SMILES']]
        # strip counterions
        if strip_counterions:
            mol_list = [get_largest_fragment(x) for x in mol_list]
        # align structures on the MCS
        if align_mcs and len(mol_list) > 1:
            mcs = rdFMCS.FindMCS(mol_list)
            mcs_query = Chem.MolFromSmarts(mcs.smartsString)
            AllChem.Compute2DCoords(mcs_query)
            for m in mol_list:
                AllChem.GenerateDepictionMatching2DStructure(m, mcs_query)
        legends = [f'P(active): {x:.4f}' for x in sel_df['P(active)']]
        if hdbscan:
            legends = [f'P(active): {x:.4f}\n strength of membership: {y:.3f}' for x,y in zip(sel_df['P(active)'], sel_df['cluster_prob'])]
        img = Draw.MolsToGridImage(mol_list, molsPerRow=5, maxMols=40, legends=legends, useSVG=True)
    return img


def cluster_results(mat, cluster_labels, BB_info):
    '''
    Input
    -----
    mat : array
        matrix of pairwise distances of compounds
    cluster_labels : array
        cluster labels assigned to each data point
    BB_info : dataframe
        contains the SMILES and P(active) values for the BBs corresponding to the inputted distance matrix
    
    Output
    ------
    cluster_data : dataframe
        SMILES and cluster id for all compounds
    cluster_df : dataframe 
        SMILES of cluster centroid, cluster label and the total number of compounds for each cluster
    '''
    # Index dictionary to access cluster labels of each point and calculate silhouette scores
    ind = np.where(cluster_labels != -1)[0]
    sil_samples = silhouette_samples(mat, cluster_labels, metric='precomputed')
    
    # Print relevant information for the specific clustering initialization chosen
    avg_score = np.mean(sil_samples[ind]) * len(sil_samples[ind])/len(BB_info)
    print('Number of clusters:', len(np.unique(cluster_labels[ind])))
    print(f'Fraction of clustered compounds: {len(ind)}/{len(BB_info)}')
    print(f'Adjusted avg silhouette score: {avg_score:.3f}')
    
    cluster_data = BB_info.copy(deep=True)
    cluster_data['Cluster'] = cluster_labels + 1

    # Create a dataframe that groups labels together -- this is the data that we visualize later
    cluster_rows = []
    for index, (k,v) in enumerate(cluster_data.groupby('Cluster')):
        cluster_rows.append([v.SMILES.values[0], k, len(v)])
    cluster_df = pd.DataFrame(cluster_rows, columns=['SMILES', 'Cluster', 'Num'])
    
    return cluster_data, cluster_df

def stacked_plots(data, bins_x, logy):
    '''
    A function to visualize the distribution of P(active) values across the different clusters formed
    
    Input
    -----
    data : dataframe
        SMILES of each compound and the cluster it is assigned to
    bb_pos : str
        specify which building block position is being analyzed (changes P(active) bins)
    bins_x : bool
        set whether activity bins are on the x-axis
    logy : bool
        set whether the y-axis is on a log scale
        
    Output
    ------
    Stacked bar plot of the distribution of P(active) values in each cluster
    '''
    # Set bins to correspond to the range of P(active) values sampled for that building block position
    span = np.linspace(np.min(data['P(active)']), np.max(data['P(active)']), 6)
    icr = span[1] - span[0]
    bins = np.arange(-icr, np.max(data['P(active)'])+icr-0.001, icr)
    labels = ['0.00'] + [f'({bins[i]:.2f}, {bins[i+1]:.2f}]' for i in range(1, len(bins)-1)]
    data['bins'] = pd.cut(data['P(active)'], bins)

    stacks = data.groupby(['Cluster', 'bins'], as_index=True).size().unstack()
    
    if bins_x:
        fig, axs = plt.subplots(figsize=(12,13))
        # Do not plot the cluster corresponding to noise
        stacks[1:].T.plot(kind='bar', stacked=True, color=plt.cm.rainbow(np.linspace(0,1,len(stacks[1:]))), ax=axs, logy=logy)
        axs.set_xticklabels(labels=labels)
        #axs.set_xticklabels(labels=['0.0', '(0.0, 0.2]', '(0.2, 0.4]', '(0.4, 0.6]', '(0.6, 0.8]', '(0.8, 1.0]'])
        axs.set_ylabel('Number of BBs')
        axs.set_title('Proportion of cluster in each activity bin')
        axs.legend(ncol=7, loc='upper right')
        
    else:
        fig, axs = plt.subplots(figsize=(20,12))
        stacks[1:].plot(kind='bar', stacked=True, logy=logy, ax=axs)
        axs.legend(labels=labels)
        #axs.legend(labels=['0.0', '(0.0, 0.2]', '(0.2, 0.4]', '(0.4, 0.6]', '(0.6, 0.8]', '(0.8, 1.0]'])
        axs.set_ylabel('Number of BBs')
        axs.set_title('Distribution of activity of BBs in each cluster')
    return None
 
def HDBSCAN(mat, BB_info, **kwargs):
    '''
    Runs HDBSCAN clustering for a given list of compounds and their distances to each other
    
    Input
    -----
    mat : array
        matrix of pairwise distances of compounds
    BB_info : dataframe
        contains the SMILES and P(active) values for the BBs corresponding to the inputted distance matrix
    **kwargs
        optional keyword arguments to tune the hdbscan algorithm; see hdbscan.HDBSCAN documentation for more info
        
    Output
    ------
    data : dataframe
        SMILES and cluster id for all compounds
    cluster : dataframe 
        SMILES of cluster centroid, cluster label and the total number of compounds for each cluster
    '''
    # initialize clustering algorithm and save cluster assignments and probabilities
    clusterer = hdbscan.HDBSCAN(metric='precomputed', **kwargs)
    clusterer.fit(mat)
    #clusterer.condensed_tree_.plot()
    data, clusters = cluster_results(mat, clusterer.labels_, BB_info)
    data['cluster_prob'] = clusterer.probabilities_
    
    # Plot two different plots for clustering analysis
    stacked_plots(data, bins_x=True, logy=True)
    stacked_plots(data, bins_x=False, logy=False)
    
    return data, clusters

def dbscan_gridsearch(mat, eps_range, min_samp_range):
    '''
    Performs a grid search to help determine the optimal set of DBSCAN parameters
    
    Input
    -----
    mat : array
        2D matrix containing the pairwise distances of all compounds
    eps_range : array
        values of the eps parameter to loop over
    min_samp_range : array
        values of the min_samples parameter to loop over
        
    Output
    ------
    dbscan_dict : dictionary
        contains the cluster labels for each compound for a given clustering parameterization
        keys for the dictionary are strings of the pair of parameter values used (e.g. dict['eps=0.20, m=5'])  
    '''
    # Grid search for optimal DBSCAN parameters for 3D sim scoring
    params = []
    sil_scores = []
    n_clusters = []
    avg_cluster_size = []
    dbscan_dict = {}

    for index, eps in enumerate(eps_range):
        if index % 10 == 0:
            print(f"{eps:.2f}")
        for min_samp in min_samp_range:
            try:
                np.fill_diagonal(mat, 0)
                # Run the clustering algorithm with specified parameters
                db = DBSCAN(eps=eps, min_samples=min_samp, metric='precomputed').fit(mat)
                # Calculate silhouette scores for each point in each formed cluster
                sil_samples = silhouette_samples(mat, db.labels_, metric='precomputed')
                # Extract indices of only points that are clustered
                ind = np.where(db.labels_ != -1)[0]
                # Store DBSCAN parameters
                params.append([eps, min_samp])
                # Store the average silhouette score across clustered points only 
                # and weight against fraction of the total number of compounds clustered  
                sil_scores.append(np.mean(sil_samples[ind]))
                # Store the number of clusters formed and average cluster size
                cluster_labels, counts = np.unique(db.labels_, return_counts=True)
                n_clusters.append(len(cluster_labels) - 1)
                avg_cluster_size.append(np.mean(counts[1:]))
                # Store cluster labels for each point in each different parameterization of DBSCAN
                dbscan_dict[f'eps={eps:.2f}, m={min_samp}'] = db.labels_
            except ValueError:
                pass
    
    # Generate scatter plot for the results of the grid search
    eps_= [x[0] for x in params]
    min_samp_ = [x[1] for x in params]
    
    # Plot the results of the parameter search
    fig, axs = plt.subplots(figsize=(20,12), dpi=150)
    im = axs.scatter(eps_, min_samp_, s=avg_cluster_size, c=sil_scores, cmap='cool')
    axs.set_yticks(np.arange(0, 20, 2))
    axs.set_ylim([np.min(min_samp_range)-1, np.max(min_samp_range)+1])
    axs.set_xlabel('eps')
    axs.set_ylabel('min_samples')
    axs.set_title('DBSCAN grid search of 3D sim scoring')
    fig.colorbar(im, ax=axs)

    for i, txt in enumerate(sil_scores):
        axs.annotate(f'{txt:.3f}', (eps_[i]+0.003, min_samp_[i]))

    plt.show()
    
    return dbscan_dict

def butina_search(mat, BB_info, dist_range):
    '''
    Performs a search to find the most optimal distance threshold
    
    Input
    -----
    mat : array
        2D matrix containing the pairwise distances of all compounds
    BB_info : dataframe
        contains the SMILES and P(active) values for the BBs corresponding to the inputted distance matrix
   dist_range : array
        values of the distance threshold to loop over

    Output
    ------
    butina_dict : dictionary
        contains the cluster labels for each compound for a given clustering parameterization
        keys for the dictionary are strings of the distance value (e.g. dict['dist=0.20'])
    '''
    dists_ = mat[np.tril_indices(len(mat), k=-1)]
    butina_dict = {}
    for index, thresh in enumerate(dist_range):
        if index % 10 == 0:
            print('{:.2f}'.format(thresh))
        # Run clustering
        but = Butina.ClusterData(data=dists_, nPts=len(mat), distThresh=thresh, isDistData=True)
        # Extract cluster labels; similar to DBSCAN, assign all unclustered compounds the label of "-1"
        but_labels = np.zeros(len(mat), dtype=int)
        cluster_id = 0
        for cluster in but:
            if len(cluster) == 1:
                but_labels[cluster[0]] = -1
            else:
                for ind in cluster:
                    but_labels[ind] = cluster_id
                cluster_id += 1
        butina_dict[f'dist={thresh:.2f}'] = but_labels
    
    raw_sil_scores = []
    adj_sil_scores = []
    
    for key in butina_dict: 
        cluster_id, counts = np.unique(butina_dict[key], return_counts=True)
        avg_cluster_size = np.mean(counts)
        scores = silhouette_samples(mat, butina_dict[key], metric='precomputed')
        ind = np.where(butina_dict[key] != -1)[0]
        # Get the average silhouette score for all clustered points
        raw_sil_scores.append(np.mean(scores[ind]))
        # Weight silhouette score by the number of points clustered
        adj_sil_scores.append(np.mean(scores[ind] * np.sum(counts[1:])/len(BB_info)))
        
    plt.plot(dist_range, raw_sil_scores, label='raw scores')
    plt.plot(dist_range, adj_sil_scores, label='adj scores')
    plt.xlabel('distance threshold')
    plt.ylabel('silhouette score')
    plt.title('Butina grid search of distance threshold')
    plt.legend(loc='best')
    plt.show()
    
    return butina_dict

