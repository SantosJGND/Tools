import numpy as np 
import collections
import itertools as it
import scipy 

def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)

from sklearn.neighbors import KernelDensity
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.cluster import MeanShift, estimate_bandwidth

from structure_tools.AMOVA_func import amova_cofactor, AMOVA_FM42
from plotly import tools
import plotly.graph_objs as go
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

import pandas as pd


def avoid_idx(data,char_avoid= [9,1],ind_thresh= 0.01):
    """
    get index of samples with proportion of features in char_avoid > ind_thresh.
    """
    nstat= np.zeros(data.shape)
    for ix in char_avoid:
        nm= data == ix
        nstat+= nm

    nstat= np.array(nstat,dtype= int)

    nstat= np.sum(nstat,axis= 1)
    nstat= nstat / float(data.shape[1])

    nstat= nstat >= ind_thresh

    nstat= np.array(nstat,dtype= int)
    #
    return nstat


def Windows_KDE_amova(SequenceStore,admx_lib,refs_lib,ncomps= 4,clsize= 15,
                    supervised= True,amova= True,Bandwidth_split= 20,
                    char_avoid= [9,1],ind_thresh= 0.01):
    '''
    Perform PCA + Mean Shift across windows. Extract Meanshift p-value vectors. Perform amova (optional).
    '''
    from IPython.display import clear_output
    Geneo= admx_lib

    Geneo_order= list(Geneo.keys())
    ref_order= list(refs_lib.keys())

    Whose= list(range(sum([len(x) for x in Geneo.values()])))
    Sup_labels= list(np.repeat(Geneo_order,[len(Geneo[x]) for x in Geneo_order]))

    ### Define parameters and libraries of analyses.

    Results = {x:recursively_default_dict() for x in SequenceStore.keys()}

    Construct = recursively_default_dict()
    PC_var= recursively_default_dict()

    for CHR in SequenceStore.keys():
        print('going on CHR: '+ str(CHR))
        for c in SequenceStore[CHR].keys():

            ### PCA and MeanShift of information from each window copied from *FM36_Galaxy.py.
            Sequences= [SequenceStore[CHR][c][x] for x in Whose]
            Sequences= np.array(Sequences) 

            Sequences= np.nan_to_num(Sequences)
            nstat= avoid_idx(Sequences,char_avoid= char_avoid,ind_thresh= ind_thresh)
            avoid= [x for x in range(len(nstat)) if nstat[x] == 1]
            keep= [x for x in range(len(nstat)) if nstat[x] == 0]
            
            if len(keep) < 10:
                continue

            pca = PCA(n_components=ncomps, whiten=False,svd_solver='randomized').fit(Sequences[keep])
            data = pca.transform(Sequences)
            PC_var[CHR][c]= [x for x in pca.explained_variance_]

            params = {'bandwidth': np.linspace(np.min(data), np.max(data),Bandwidth_split)}
            grid = GridSearchCV(KernelDensity(algorithm = "ball_tree",breadth_first = False), params,verbose=0)

            ######################################
            ####### TEST global Likelihood #######
            ######################################
            Focus_labels = [z for z in it.chain(*refs_lib.values()) if z not in avoid]

            #### Mean Shift approach
            ## from sklearn.cluster import MeanShift, estimate_bandwidth

            bandwidth = estimate_bandwidth(data, quantile=0.2, n_samples=len(Focus_labels))
            if bandwidth <= 1e-3:
                bandwidth = 0.1

            ms = MeanShift(bandwidth=bandwidth, cluster_all=False, min_bin_freq=clsize)
            ms.fit(data[Focus_labels,:])
            labels = ms.labels_


            Tree = {x:[Focus_labels[y] for y in range(len(labels)) if labels[y] == x] for x in [g for g in list(set(labels)) if g != -1]}
            Keep= [x for x in Tree.keys() if len(Tree[x]) > clsize]

            Tree= {x:Tree[x] for x in Keep}
            Ngps= len(Tree)
            SpaceX = {x:data[Tree[x],:] for x in Tree.keys()}

            ### Extract MScluster likelihood by sample

            for hill in SpaceX.keys():

                grid.fit(data[Tree[hill],:])

                # use the best estimator to compute the kernel density estimate
                kde = grid.best_estimator_

                # normalize kde derived log-likelihoods, derive sample p-values
                P_dist = kde.score_samples(data[Tree[hill],:])
                Dist = kde.score_samples(data)
                P_dist= np.nan_to_num(P_dist)
                Dist= np.nan_to_num(Dist)
                if np.std(P_dist) == 0:
                    Dist= np.array([int(Dist[x] in P_dist) for x in range(len(Dist))])
                else:
                    Dist = scipy.stats.norm(np.mean(P_dist),np.std(P_dist)).cdf(Dist)
                Dist= np.nan_to_num(Dist)
                Dist[avoid]= 0
                Construct[CHR][c][hill] = Dist


                ######################################### 
            ############# AMOVA ################
                #########################################
            if supervised:
                labels= Sup_labels
                Who= [z for z in it.chain(*[refs_lib[x] for x in ref_order])]
                labels= [labels[x] for x in Who]
                Ngps= len(refs_lib)

            else:
                Who = [x for x in range(len(labels)) if labels[x] != -1 and labels[x] in Keep]
                labels = [labels[x] for x in Who]
                Who= [Focus_labels[x] for x in Who]

            if amova:
                clear_output()
                Bool_set= Sequences[Who,:].astype(bool)
                print('chr {}, where: {}, supervised: {}, n clusters: {}'.format(CHR,c,str(supervised),Ngps))
                Amova1,Ciggy= AMOVA_FM42(Bool_set,labels,n_boot=0,metric= 'jaccard')
                Amova2,Ciggy= AMOVA_FM42(data[Who,:],labels,n_boot=0,metric= 'euclidean')
                Amova3,Ciggy= AMOVA_FM42(Bool_set,labels,n_boot=0,metric= 'hamming')
                print('old: ; jaccard: {}; PCA euc: {}; nHam: {}'.format(Amova1,Amova2,Amova3))
                Results[CHR][c] = [Ngps,Amova1,Amova2,Amova3]
    
    return Results, Construct, PC_var



def Haplotype_MSlabs(Haplotypes,label_vector,ncomps= 3):

    print('Number of markers: {}'.format(Haplotypes.shape[1]))

    print('Number of individuals: {}'.format(Haplotypes.shape[0]))

    pca2 = PCA(n_components=ncomps, whiten=False,svd_solver='randomized')

    feats= pca2.fit_transform(Haplotypes)

    coords= {z:[x for x in range(len(label_vector)) if label_vector[x] == z] for z in list(set(label_vector))}

    ## perform MeanShift clustering.
    bandwidth = estimate_bandwidth(feats, quantile=0.25)

    ms = MeanShift(bandwidth=bandwidth, bin_seeding=False, cluster_all=False, min_bin_freq=35)
    ms.fit(feats)
    labels1 = ms.labels_
    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}

    ##

    PCA_color_ref= ['darkseagreen','crimson', 'cyan', 'darkblue', 'darkcyan',
                'darkgoldenrod', 'darkgray', 'darkgrey', 'darkgreen',
                'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange',
                'darkorchid', 'darkred', 'darksalmon', 'darkseagreen',
                'darkslateblue', 'darkslategray', 'darkslategrey',
                'darkturquoise', 'darkviolet', 'deeppink']

    from plotly import tools

    fig_pca_subplots = tools.make_subplots(rows=2, cols=2)

    for subp in range(4):

        n_plot= subp

        if subp >= 2:
            coords= label_select
            subp -= 2
            Col_vec= PCA_color_ref

        else:
            Col_vec= ['red','yellow','blue','black','green']


        for i in coords.keys():
            trace= go.Scatter(
            x = feats[coords[i],0],
            y = feats[coords[i],subp + 1],
            mode= "markers",
            name= str(i),
            marker= {
            'color': Col_vec[i],
            'line': {'width': 0},
            'size': 6,
            'symbol': 'circle',
            "opacity": .8})

            fig_pca_subplots.append_trace(trace, int(n_plot/float(2)) + 1, subp + 1)

        fig_pca_subplots['layout']['yaxis' + str(n_plot + 1)].update(title='PC{}'.format(subp + 2))
        fig_pca_subplots['layout']['xaxis' + str(n_plot + 1)].update(title='PC1')

    layout = go.Layout()

    fig= go.Figure(data=fig_pca_subplots, layout=layout)
    
    return fig, feats, label_select


def MAC_process(Construct,Out,Cl_store,refs_lib,Fam,Names= [],target_var= [],Dr_var= 'all',focus_subset= False,Focus= [],Dr_dim= 4,threshold= 0.1,Method= 'MeanShift'):

    Coordinates = [[[[CHR,bl,Out[CHR][bl],x] for x in Construct[CHR][bl].keys()] for bl in sorted(Construct[CHR].keys())] for CHR in sorted(Construct.keys())]
    Coordinates = [z for z in it.chain(*[y for y in it.chain([x for x in it.chain(*Coordinates)])])]


    Coordinates= np.array(Coordinates)
    

    Clover= [[[Construct[CHR][bl][x] for x in Construct[CHR][bl].keys()] for bl in sorted(Construct[CHR].keys())] for CHR in sorted(Construct.keys())]
    Clover= [z for z in it.chain(*[y for y in it.chain(*Clover)])]
    Clover= np.array(Clover)
    Clover.shape
    
    Membership=[]
    
    for CHR in sorted(Construct.keys()):
        for bl in sorted(Construct[CHR].keys()):
            
            Bls= sorted(list(Construct[CHR][bl].keys()))
            pVals= np.array([Construct[CHR][bl][y] for y in Bls])
            
            max_vals= np.amax(pVals,axis= 0)
            max_indx= np.argmax(pVals,axis= 0)
            
            inlier= [x for x in range(pVals.shape[1]) if max_vals[x] >= threshold]
            
            BL_select= list(set([max_indx[x] for x in inlier]))
            
            #print('clusters {} selected. {} %'.format(BL_select,len(BL_select)/float(len(Bls))))
            
            if not BL_select:
                Empty.append([CHR,bl])
                continue
            
            BL_select= { 
                x: pVals[x] for x in BL_select
                }
            
            Assignment= {
                    Bls[b]: [x for x in inlier if max_indx[x] == b] for b in BL_select.keys()
                }
            
            for cl in Bls:
                if cl not in BL_select.keys():
                    vector= ''
                else:
                    vector= '.'.join([str(x) for x in Assignment[cl]])
                
                Membership.append(vector)

    
    #Membership= np.array(Membership)
    
    Coordinates= pd.DataFrame(Coordinates,columns= ['chrom','start','end','bl'])
    Coordinates['members']= Membership
    

    from sklearn import preprocessing

    Clover = np.nan_to_num(Clover)
    preProc_Clover = Clover

    print('Clover shape: ', Clover.shape)

    Clover = preprocessing.scale(Clover,axis = 1)
    #

    print("Clover shape: ", Clover.shape)


    reefer= [g for g in it.chain(*[refs_lib[y] for y in sorted(refs_lib.keys())])]

    if not focus_subset:
        Subset= list(range(Clover.shape[1]))
    else:
        Subset= [Names.index(x) for x in Focus]

    ## apply pca to reference accessions, transform the rest.

    Dr_processes= ['target','focus_inc','all']

    if Dr_var not in Dr_processes:
        print('Dr_process selected: {}, Dr_var processes available: {}'.format(Dr_var,Dr_processes))
        Dr_var= 'target'

    print('focusing Dr on {}'.format(Dr_var))

    if Dr_var== 'target':
        variation_focus= [Names.index(Fam[x]) for x in it.chain(*[refs_lib[z] for z in target_var])]

    if Dr_var== 'focus_inc':
        variation_focus= [Names.index(x) for x in Focus]
        variation_focus.extend([Names.index(Fam[x]) for x in it.chain(*[refs_lib[z] for z in target_var])])

    if Dr_var== 'all':
        variation_focus= list(range(Clover.shape[1]))


    ### PCA
    pca = PCA(n_components=Dr_dim, whiten=False).fit(Clover[:,variation_focus].T)
    X_se = pca.transform(Clover[:,Subset].T)
    COMPS = pca.components_.T*np.sqrt(pca.explained_variance_)


    ###############################################################################
    ########################### PAINTING SHIT!! ###################################
    ###############################################################################

    ## 
    ## CLUSTER EIGENVALUES
    ##

    bandwidth = estimate_bandwidth(COMPS, quantile=0.1)
    if bandwidth==0:
        bandwidth = 0.1

    func_cl= Cl_store[Method]['Clusterfunc']
    func_kwargs= Cl_store[Method]['cluster_kwargs']


    Clusterfunck= func_cl(**func_kwargs,bandwidth= bandwidth) 
    Clusterfunck.fit(COMPS)

    labels1 = Clusterfunck.labels_
    label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}


    ###############################################################################
    #### Average normalized likelihhod among clustered eigenvectors by haplotype #####
    ###############################################################################


    Cameo = []

    for cramp in sorted(label_select.keys()):
        Clamp = np.mean(preProc_Clover[label_select[cramp],:],axis = 0)
        Fry = [Clamp[x] for x in Subset]
        Cameo.append(Fry)

    Cameo = np.array(Cameo).T


    ###########################################################################
    ### cosine of the clustered eigenvectors with haplotype coordinates ######## DEPRECATED
    ###########################################################################

    #cos_threshold = .6
    #
    #
    #from numpy import dot
    #from numpy.linalg import norm
    #
    #SpaceY = recursively_default_dict()
    #
    #for g in label_select.keys():
    #    Green = COMPS[label_select[g],:]
    #    SpaceY[g] = [mean(dot(X_se[x],Green.T)/norm(Green,axis=1)/norm(X_se[x])) for x in range(X_se.shape[0])]
    #
    #Globe = np.array([SpaceY[x] for x in sorted(SpaceY.keys())]).T
    #
    #

    ######## Reducing the number of cluster profiles to print:
    new_labs= labels1

    return preProc_Clover, Cameo, Coordinates, COMPS, X_se, label_select, Subset, labels1



def KDE_pca(feats= [],Cameo= [],label_vector= [],Subset= [],Col_vec= [],height= 2000,width= 1000):
    
    Ncols= 2
    titles=['Global']
    titles.extend(['cluster ' + str(x + 1) for x in range(Cameo.shape[1])])
    titles= list(np.repeat(titles,2))
    print(titles)

    fig_pca_subplots = tools.make_subplots(rows= int(len(titles) / float(Ncols)) + (len(titles) % Ncols > 0), cols=Ncols,
                             subplot_titles=tuple(titles))
    
    #####
    for subp in range(len(titles)):
        #print(subp)
        
        pos1= int(float(subp) / Ncols) + 1

        pos2= subp % Ncols + 1
        
        n_plot= subp
        
        
        if subp >= 2:
            gradient= Cameo[:,pos1-2]

            trace= go.Scatter(
            x = feats[:,0],
            y = feats[:,pos2],
            mode= "markers",
            marker= {
                'color': gradient,
                'colorscale':'Viridis',
                'line': {'width': 0},
                'size': 6,
                'symbol': 'circle',
                "opacity": .6
            })

            fig_pca_subplots.append_trace(trace, row=pos1, col=pos2)

        else:
            coords= {z:[x for x in Subset if label_vector[x] == z] for z in list(set(label_vector))}

            for i in coords.keys():
                if coords[i]:
                    trace= go.Scatter(
                    x = feats[coords[i],0],
                    y = feats[coords[i],pos2],
                    mode= "markers",
                    name= str(i),
                    marker= {
                    'color': Col_vec[i],
                    'line': {'width': 0},
                    'size': 6,
                    'symbol': 'circle',
                    "opacity": .8})

                    fig_pca_subplots.append_trace(trace, row=pos1, col=pos2)

        fig_pca_subplots['layout']['yaxis' + str(n_plot + 1)].update(title='PC{}'.format(pos2+1))
        fig_pca_subplots['layout']['xaxis' + str(n_plot + 1)].update(title='PC1')
    
    fig_pca_subplots['layout'].update(height= height,width= width)
    
    #fig= go.Figure(data=fig_pca_subplots, layout=layout)
    iplot(fig_pca_subplots)


def MS_ideogram(gp,MS_threshold,Clover,Coordinates,label_select,Out,ideo_order= [],Chr= 1,height_chrom= .5,height= 10,width= 5):
    
    Where= Coordinates[label_select[gp-1]]

    Clover_select= Clover[label_select[gp-1]]
    
    if ideo_order:
        Clover_select= Clover_select[:,ideo_order]

    Clover_select= Clover_select > MS_threshold

    Clover_select.astype(np.int)
    
    ##
    trace0 = go.Scatter(
        x=[0, max(Out[Chr].keys())],
        y=[0, height_chrom * Clover_select.shape[1]],
        mode='text',
    )

    data = [trace0]
    layout = {
        'xaxis': {
            'showgrid': False,
        },
        'yaxis': {
        },
        'shapes': []
    }

    for row in range(Where.shape[0]):
            site= Where[row]
            # filled Rectangle
            CHR= site[0]
            start= site[1]
            end= site[2]
            cl= site[3]

            for v in np.where(Clover_select[row,:] == 1)[0]:
                v= Clover_select.shape[1]-v-1
                rekt= {
                    'type': 'rect',
                    'y0': v * height_chrom,
                    'x0': start,
                    'y1': (v + 1)*height_chrom,
                    'x1': end,
                    'line': {
                        'color': 'rgba(128, 0, 128, 1)',
                        'width': 2,
                    },
                    'fillcolor': 'rgba(128, 0, 128, 0.7)',
                }

                layout['shapes'].append(rekt)



    fig = {
        'data': data,
        'layout': layout,
    }
    iplot(fig)


def recover_MSgp(Construct,COMPS,Out,labels1,refs_lib,Fam,nneighbours= 5,Names= [],target_var= [],Dr_var= 'all',focus_subset= False,Focus= [],Dr_dim= 3,Method= 'MeanShift'):
    Coordinates = [[[[CHR,bl,Out[CHR][bl],x] for x in Construct[CHR][bl].keys()] for bl in sorted(Construct[CHR].keys())] for CHR in sorted(Construct.keys())]
    Coordinates = [z for z in it.chain(*[y for y in it.chain([x for x in it.chain(*Coordinates)])])]


    Coordinates= np.array(Coordinates)

    Clover= [[[Construct[CHR][bl][x] for x in Construct[CHR][bl]] for bl in sorted(Construct[CHR].keys())] for CHR in sorted(Construct.keys())]
    Clover= [z for z in it.chain(*[y for y in it.chain(*Clover)])]
    Clover= np.array(Clover)
    Clover.shape

    from sklearn import preprocessing

    Clover = np.nan_to_num(Clover)
    preProc_Clover = Clover

    print('Clover shape: ', Clover.shape)

    Clover = preprocessing.scale(Clover,axis = 1)
    #

    print("Clover shape: ", Clover.shape)


    reefer= [g for g in it.chain(*[refs_lib[y] for y in sorted(refs_lib.keys())])]

    if not focus_subset:
        Subset= list(range(Clover.shape[1]))
    else:
        Subset= [Names.index(x) for x in Focus]

    ## apply pca to reference accessions, transform the rest.

    Dr_processes= ['target','focus_inc','all']

    if Dr_var not in Dr_processes:
        print('Dr_process selected: {}, Dr_var processes available: {}'.format(Dr_var,Dr_processes))
        Dr_var= 'target'

    print('focusing Dr on {}'.format(Dr_var))

    if Dr_var== 'target':
        variation_focus= [Names.index(Fam[x]) for x in it.chain(*[refs_lib[z] for z in target_var])]

    if Dr_var== 'focus_inc':
        variation_focus= [Names.index(x) for x in Focus]
        variation_focus.extend([Names.index(Fam[x]) for x in it.chain(*[refs_lib[z] for z in target_var])])

    if Dr_var== 'all':
        variation_focus= list(range(Clover.shape[1]))


    ### PCA
    pca = PCA(n_components=Dr_dim, whiten=False).fit(Clover[:,variation_focus].T)
    new_X_se = pca.transform(Clover[:,Subset].T)
    new_COMPS = pca.components_.T*np.sqrt(pca.explained_variance_)
    print(new_COMPS.shape)


    from sklearn.neighbors import KNeighborsClassifier
    neigh = KNeighborsClassifier(n_neighbors=nneighbours)
    neigh.fit(COMPS, labels1) 

    new_labs= neigh.predict(new_COMPS)
    nlab_select= {z:[x for x in range(new_COMPS.shape[0]) if new_labs[x] == z] for z in list(set(labels1))}

    MS_ideogram(gp,
                MS_threshold,
                preProc_Clover,
                Coordinates,
                nlab_select,
                Out,
                ideo_order= ideo_order,
                Chr= Chr,
                height_chrom= height_chrom,
                height= height,
                width= width)
