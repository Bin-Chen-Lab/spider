import numpy as np
import torch.nn as nn
import pandas as pd
import matplotlib.pyplot as plt
import torch.nn.functional as F
import torch
import torch.optim as optim
import sys
from scipy import stats
import os
import random
import copy
import time
from sklearn.linear_model import LinearRegression
import pickle
from numpy.linalg import norm
import collections


def predict_unseen(tissue, disease, save_path, SPIDER_model_file_path, use_pretrain, cell_type_file,
                   all_trainable_proteins_gene_names_6_training_sets, training_epoch, all_test_proteins_gene_names,
                   file_A_B_C_matching_table):
    #--------------------------------------------------------------------------------------
    #select ensemble members:
    if use_pretrain == 'T':
        match_training_protein_gene_name = pd.read_csv(SPIDER_model_file_path + 'protein_gene_names_union_289_DNNs_from_combined_6_training_sets_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_internal_val_threshold_0.6_20230115.csv')
        training_protein_DNN_internal_val = pd.read_csv(SPIDER_model_file_path + 'cor_per_pro_internal_val_6_combined_training_sets_cell_features_protein_specific_DNN_onehot_celltype_tissue_disease_SCANVI_128dim_64_32_16_0.0001_seen_proteins_20230115.csv', index_col = ['Unnamed: 0'])
    else:
        match_training_protein_gene_name = all_trainable_proteins_gene_names_6_training_sets
        training_protein_DNN_internal_val = pd.read_csv(SPIDER_model_file_path + 'retrain_internal_val_performance.csv')
    combined_6_training_sets_gene_cor_mat = pd.read_csv(SPIDER_model_file_path + 'unseen_file/shared_gene_features_20230115.csv')
    threshold = 0.6
    training_protein_DNN_internal_val = training_protein_DNN_internal_val.iloc[:training_protein_DNN_internal_val.shape[0]-1, :]
    select_protein_DNN_good = training_protein_DNN_internal_val.loc[training_protein_DNN_internal_val['pearson'] > threshold,:].index
    match_training_protein_gene_name = match_training_protein_gene_name[{'consistent_protein_name', 'gene_name'}]
    match_training_protein_gene_name = match_training_protein_gene_name.drop_duplicates() #289 proteins
    match_training_protein_gene_name.index = match_training_protein_gene_name['consistent_protein_name']
    #Continue to select training proteins with existed z features in the gene co-expression matrix.
    select_final_ensemble_members = []
    for i in range(0, len(select_protein_DNN_good)):
        if(match_training_protein_gene_name.loc[select_protein_DNN_good[i], 'gene_name'] in combined_6_training_sets_gene_cor_mat['gene'].tolist()):
            select_final_ensemble_members.append(select_protein_DNN_good[i])
    match_training_protein_gene_name_2 = match_training_protein_gene_name.loc[select_final_ensemble_members, :]
    match_training_protein_gene_name_2 = match_training_protein_gene_name_2.drop(match_training_protein_gene_name_2.loc[match_training_protein_gene_name_2['gene_name'].duplicated(),:].index) #delete duplicated genes among ensemble members
    #--------------------------------------------------------------------------------------
    use_n_ensemble = 8
    threshold_DNN_internal_val_acc = '0.6'
    RNA_SCANVI_latent_representation_file = save_path + 'query_embeddings.csv'
    test_gene_coexp_matrix_file = save_path + 'query_gene_coexpression.csv'
    save_linear_model_pred_vs_truth_y_unseen_file_path = save_path
    #--------------------------------------------------------------------------------------
    #Prepare zero-shot learning:
    #-------------------------------------------------------------
    #Get each one of the ensemble members' ~10000-dim normalized original features (i.e., Z1, Z2, ... Z22 for P1, P2, ... P22, respectively), save them:
    x = pd.read_csv(test_gene_coexp_matrix_file, index_col = 'Unnamed: 0')
    all_trainable_proteins_gene_names_6_training_sets.index = all_trainable_proteins_gene_names_6_training_sets['consistent_protein_name']
    all_trainable_proteins_gene_names_6_training_sets = all_trainable_proteins_gene_names_6_training_sets.loc[match_training_protein_gene_name_2['consistent_protein_name'], :]
    all_trainable_proteins_gene_names_6_training_sets_2 = all_trainable_proteins_gene_names_6_training_sets
    all_trainable_proteins_gene_names_6_training_sets_2.index = all_trainable_proteins_gene_names_6_training_sets_2['gene_name']
    #Select ensemble members with corresponding genes existing in external validation set's coexpression file:
    shared_gene_ensemble_members_test_set_coexp = list(set(x.columns) & set(all_trainable_proteins_gene_names_6_training_sets_2['gene_name']))
    #Select ensemble members with corresponding genes existing in external validation set's coexpression file
    shared_training_test_set_coexp_features = list(set(combined_6_training_sets_gene_cor_mat['gene']) & set(x.columns))
    z = x[shared_training_test_set_coexp_features]
    z_ensemble = z.loc[shared_gene_ensemble_members_test_set_coexp, :] #113 x 10702
    #-------------------------------------------------------------
    #Get each one of the tested unseen proteins' ~10000-dim normalized original features (i.e., Z1, Z2, ... Z22 for P1, P2, ... P22, respectively), save them:
    shared_gene_unseen_proteins_test_set_coexp = list(set(x.columns) & set(all_test_proteins_gene_names))
    if len(shared_gene_unseen_proteins_test_set_coexp) == 0:
        sys.exit('No predictable unseen protein')
    z_unseen_proteins = z.loc[shared_gene_unseen_proteins_test_set_coexp, :] #3927 x 10702
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
    #Impute unseen proteins:
    all_protein_list = shared_gene_unseen_proteins_test_set_coexp
    def softmax(x):
        """Compute softmax values for each sets of scores in x."""
        return np.exp(x) / np.sum(np.exp(x), axis=1)
    summary_max_coef = []
    all_linear_model_pred_y_unseen = None
    for protein_name in all_protein_list: #among tested unseen proteins
        gene_exp_file = pd.read_csv(RNA_SCANVI_latent_representation_file, index_col = 'Unnamed: 0')
        celltype = pd.read_csv(cell_type_file, index_col = 'Unnamed: 0')
        protein_list = [protein_name]   
        #----------------------------------------    
        X = gene_exp_file
        test_set_unseen_y_pred_combine = None
        coef_col_names = list(set(z_ensemble.index).difference(set({protein_name})))
        #onehot encoding of cell type, tissue and disease:
        test_set_cell_type_class = np.array([100000] * X.shape[0], dtype=object).T
        #onehot encoding of tissue:
        all_tissue_type = list(dict(collections.Counter(file_A_B_C_matching_table['tissue_type'])))
        if not tissue in all_tissue_type:
            test_set_onehot_tissue = np.array([100000] * X.shape[0], dtype=object).T
        else:
            for t in all_tissue_type:
                if tissue == t:
                    test_set_onehot_tissue = np.array([list(file_A_B_C_matching_table.loc[file_A_B_C_matching_table['tissue_type'] == t]['tissue_class'])[0]] * X.shape[0], dtype=object).T
        #onehot encoding of disease:
        all_disease_type = list(dict(collections.Counter(file_A_B_C_matching_table['disease_type'])))
        if not disease in all_disease_type:
            test_set_onehot_disease = np.array([100000] * X.shape[0], dtype=object).T
        else:
            for d in all_disease_type:
                if disease == d:
                    test_set_onehot_disease = np.array([list(file_A_B_C_matching_table.loc[file_A_B_C_matching_table['disease_type'] == d]['disease_class'])[0]] * X.shape[0], dtype=object).T
        #onehot encoding of cell type:
        all_cell_type = list(dict(collections.Counter(file_A_B_C_matching_table['cell_type'])))
        for c in range(0, celltype.shape[0]):
            if celltype.iloc[c, :]['final_celltype'] in all_cell_type:
                test_set_cell_type_class[c] = list(file_A_B_C_matching_table.loc[file_A_B_C_matching_table['cell_type'] == celltype.iloc[c, :]['final_celltype']]['cell_type_class'])[0]
        #
        #-------------------------------------------------------------
        #Get cosine similarity
        Dot_product = None
        dot_product_col_names = coef_col_names
        for protein_name2 in dot_product_col_names:
            Zt = pd.DataFrame(z_unseen_proteins.loc[protein_name, :])
            Zt = Zt.dropna()
            Zi = pd.DataFrame(z_ensemble.loc[protein_name2, :])
            Zi = Zi.dropna()
            use_common_genes = list(set(Zt.index) & set(Zi.index))
            Zt = Zt.loc[use_common_genes, :]
            Zi = Zi.loc[Zt.index, :]
            if not Dot_product is None:
                Dot_product.append((np.dot(Zi.T, Zt)/(norm(Zi)*norm(Zt))).reshape(-1)[0])
            else:
                Dot_product = [(np.dot(Zi.T, Zt)/(norm(Zi)*norm(Zt))).reshape(-1)[0]]
        Dot_product = pd.DataFrame(Dot_product).T
        Dot_product.index = {protein_name}
        Dot_product.columns = dot_product_col_names
        dot_product = Dot_product
        dot_product = dot_product[coef_col_names]
        dot_product.index = {protein_name}
        dot_product = dot_product.sort_values([protein_name], axis = 1, ascending = False)
        #-------------------------------------------------------------
        #Get Sc(gj(X)):
        dot_product = dot_product.iloc[:, :use_n_ensemble]   
        dot_product_original = copy.deepcopy(dot_product)
        x0 = np.array(dot_product)
        x0 = softmax(x0)
        dot_product.iloc[0, :] = x0
        summary_max_coef.append(dot_product_original.max(axis = 1).iloc[0])
        my_intercepts = np.zeros(1)
        my_coefficients = np.array(dot_product)
        reg = LinearRegression()
        reg.intercept_ = my_intercepts
        reg.coef_ = my_coefficients
        ##############################################################################
        #get g(xi):
        for j in dot_product.columns:
            test_set_onehot = pd.concat([pd.DataFrame(test_set_onehot_tissue), pd.DataFrame(test_set_onehot_disease), pd.DataFrame(test_set_cell_type_class)], axis=1)
            if use_pretrain == 'T':
                with open(SPIDER_model_file_path + 'onehot_celltype_tissue_disease_all_training_samples_20230115_' + all_trainable_proteins_gene_names_6_training_sets_2.loc[j, 'consistent_protein_name'], 'rb') as pickle_file:
                    encoder = pickle.load(pickle_file)
            else:
                with open(SPIDER_model_file_path + 'onehot_retrain_' + all_trainable_proteins_gene_names_6_training_sets_2.loc[j, 'consistent_protein_name'], 'rb') as pickle_file:
                    encoder = pickle.load(pickle_file)
            test_set_onehot = pd.DataFrame(encoder.transform(test_set_onehot).toarray(), index = gene_exp_file.index)
            test_set_onehot.columns = encoder.get_feature_names_out()
            X = pd.concat([gene_exp_file, test_set_onehot], axis=1)
            cells=X.index
            proteins = ['protein']
            X = torch.tensor(np.vstack(X.values).astype(np.float))
            X = X.type(torch.FloatTensor)
            protein_list=['protein']
            def set_seed(seed):
                np.random.seed(seed)
                random.seed(seed)
                torch.manual_seed(seed) # cpu
                torch.cuda.manual_seed_all(seed)  # gpu
                torch.backends.cudnn.deterministic = True
            # consistent results on the cpu and gpu
            set_seed(10)
            class Net(nn.Module):
                def __init__(self):
                    super(Net, self).__init__()
                    self.fc1 = nn.Linear(X.shape[1], 64)
                    self.fc2 = nn.Linear(64, 32)
                    self.fc3 = nn.ModuleDict({})
                    self.fc3[protein_list[0]]=nn.Linear(32, 16)
                    self.fc4 = nn.ModuleDict({})
                    self.fc4[protein_list[0]]=nn.Linear(16, 1)  
                def forward(self, x):
                    x = F.relu(self.fc1(x))
                    x = F.relu(self.fc2(x))
                    outputs={}
                    outputs[protein_list[0]]=self.fc4[protein_list[0]](F.relu(self.fc3[protein_list[0]](x)))
                    return outputs 
            net = Net()
            if use_pretrain == 'T':
                net.load_state_dict(torch.load(SPIDER_model_file_path + all_trainable_proteins_gene_names_6_training_sets_2.loc[j, 'consistent_protein_name'] + '_all_training_samples_cell_features_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_DNN_SCANVI_128dim_20230115_ep' + str(training_epoch.loc[all_trainable_proteins_gene_names_6_training_sets_2.loc[j, 'consistent_protein_name'], 'training_epoch'])))
            else:
                net.load_state_dict(torch.load(SPIDER_model_file_path + all_trainable_proteins_gene_names_6_training_sets_2.loc[j, 'consistent_protein_name'] + '_model_ep' + str(training_epoch.loc[all_trainable_proteins_gene_names_6_training_sets_2.loc[j, 'consistent_protein_name'], 'training_epoch'])))                
            y_pred = net(X)
            tmp2 = torch.transpose(torch.stack(list(y_pred[protein_list[0]])),0,1).view(X.shape[0],-1)
            y_pred = pd.DataFrame(tmp2.detach().numpy().T)
            y_pred.index=proteins
            y_pred.columns=cells
            y_pred = y_pred.T
            y_pred.columns = ['y_pred']
        #-------------------------------------------------------------
            if not test_set_unseen_y_pred_combine is None:         
                tmp = y_pred
                tmp = pd.DataFrame(tmp['y_pred'])
                test_set_unseen_y_pred_combine = pd.concat([test_set_unseen_y_pred_combine, tmp], axis=1)
            else:          
                tmp = y_pred
                test_set_unseen_y_pred_combine = pd.DataFrame(tmp['y_pred'])
            print(j)
        test_set_unseen_y_pred_combine.columns = dot_product.columns
        tmp = y_pred
        linear_model_pred_y_unseen = reg.predict(test_set_unseen_y_pred_combine[dot_product.columns])
        linear_model_pred_y_unseen = pd.DataFrame(linear_model_pred_y_unseen, index = tmp.index, columns = ['Predicted'])
        #linear_model_pred_y_unseen.to_csv(save_linear_model_pred_vs_truth_y_unseen_file_path + protein_name + '_predicted_unseen.csv')
        if(all_linear_model_pred_y_unseen is None):
            all_linear_model_pred_y_unseen = linear_model_pred_y_unseen
        else:
            all_linear_model_pred_y_unseen = pd.concat([all_linear_model_pred_y_unseen, linear_model_pred_y_unseen], axis=1)
    all_linear_model_pred_y_unseen.columns = all_protein_list
    all_linear_model_pred_y_unseen.to_csv(save_path + 'all_unseen_proteins_predicted.csv')
    summary_max_coef = pd.DataFrame(summary_max_coef, index = all_protein_list)
    summary_max_coef.columns = {'max_inferred_coef'}
    summary_max_coef.to_csv(save_linear_model_pred_vs_truth_y_unseen_file_path + 'confidence_score_all_unseen_proteins.csv')




    





        





        
































