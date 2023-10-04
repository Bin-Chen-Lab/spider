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
import time
import re
import pickle
import collections


def predict_seen(tissue, disease, save_path, all_protein_list, SPIDER_model_file_path, training_epoch, use_cell_type, query_cell_type, use_pretrain, file_A_B_C_matching_table):
    all_y_pred = None
    for protein_name in all_protein_list:
        #----------------------------------------
        gene_exp_file = pd.read_csv(save_path + 'query_embeddings.csv', index_col = 'Unnamed: 0')
        if use_cell_type == 'SingleR':
            celltype = pd.read_csv(save_path + 'query_celltype_SingleR.csv', index_col = 'Unnamed: 0')
        else:
            celltype = query_cell_type
        protein_list = [protein_name]   
        #----------------------------------------    
        X = gene_exp_file
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
        test_set_onehot = pd.concat([pd.DataFrame(test_set_onehot_tissue), pd.DataFrame(test_set_onehot_disease), pd.DataFrame(test_set_cell_type_class)], axis=1)
        if use_pretrain == 'T':
            with open(SPIDER_model_file_path + 'onehot_celltype_tissue_disease_all_training_samples_20230115_' + protein_name, 'rb') as pickle_file:
                encoder = pickle.load(pickle_file)
        else:
            with open(SPIDER_model_file_path + 'onehot_retrain_' + protein_name, 'rb') as pickle_file:
                encoder = pickle.load(pickle_file)
        test_set_onehot = pd.DataFrame(encoder.transform(test_set_onehot).toarray(), index = X.index)
        test_set_onehot.columns = encoder.get_feature_names_out()
        X = pd.concat([X, test_set_onehot], axis=1)
        cells=X.index
        proteins = ['protein']
        X = torch.tensor(np.vstack(X.values).astype(np.float))
        X = X.type(torch.FloatTensor)#torch.Size([8000, 178])
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
            net.load_state_dict(torch.load(SPIDER_model_file_path + protein_name + '_all_training_samples_cell_features_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_DNN_SCANVI_128dim_20230115_ep' + str(training_epoch.loc[protein_name, 'training_epoch'])))
        else:
            net.load_state_dict(torch.load(SPIDER_model_file_path + protein_name + '_model_ep' + str(training_epoch.loc[protein_name, 'training_epoch'])))
        y_pred = net(X)
        tmp2 = torch.transpose(torch.stack(list(y_pred[protein_list[0]])),0,1).view(X.shape[0],-1)
        y_pred = pd.DataFrame(tmp2.detach().numpy().T)
        y_pred.index=proteins
        y_pred.columns=cells
        y_pred = y_pred.T
        y_pred.columns = ['y_pred']
        sys.stdout.flush()
        if(all_y_pred is None):
            all_y_pred = y_pred
        else:
            all_y_pred = pd.concat([all_y_pred, y_pred], axis=1)
    all_y_pred.columns = all_protein_list
    pd.DataFrame({'protein_name':all_protein_list}).to_csv(save_path + 'all_seen_protein_names.csv')
    all_y_pred.to_csv(save_path + 'all_seen_proteins_predicted.csv')





















