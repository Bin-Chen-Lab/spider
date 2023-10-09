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
from sklearn.preprocessing import OneHotEncoder
import pickle

def train_model(file_A_B_C_matching_table, all_protein_list, SPIDER_model_file_path): #First train on 90% samples, record the epoch, then train on all samples
    #-------------------------------------------------------------
    #all_protein_list: Names of all proteins in your reference set.
    #file_A_B_C_matching_table: A table containing column names: 'file_B', 'file_C', 'tissue_class', 'disease_class', 'cell_type_class'.
    #-------------------------------------------------------------
    #train on 90% samples:
    all_epochs = []
    all_internal_val_performance = []
    for protein_name in all_protein_list:
     select_file_A_B_C_matching_table = file_A_B_C_matching_table
     X = None
     y = None
     y_1 = None #y_1 represents protein abundance in each trained cell.
     y_2 = None #y_2 represents tissue of each trained cell.
     y_3 = None #y_3 represents disease of each trained cell.
     y_4 = None #y_4 represents celltype of each trained cell.
     for i in range(0, len(select_file_A_B_C_matching_table)):
        gene_exp_file = pd.read_csv(select_file_A_B_C_matching_table.iloc[i, :]['file_B'], index_col = 'Unnamed: 0')
        ADT = pd.read_csv(select_file_A_B_C_matching_table.iloc[i, :]['file_C'], index_col = 'Unnamed: 0')
        tissue = select_file_A_B_C_matching_table.iloc[i, :]['tissue_class']
        disease = select_file_A_B_C_matching_table.iloc[i, :]['disease_class']
        celltype = select_file_A_B_C_matching_table.iloc[i, :]['cell_type_class'] 
        if protein_name in ADT.index:
             ADT = pd.DataFrame(ADT.loc[protein_name, :]).T
             protein_list = [protein_name]
             X1 = gene_exp_file
             y1 = pd.DataFrame(np.array(ADT).reshape((1, -1)))
             y1.columns = X1.index
             y1 = y1.T #y1 represents protein abundance in each trained cell.
             y2 = pd.DataFrame({'tissue': [tissue] * X1.shape[0]}) #y2 represents tissue of each trained cell.
             y2.index = X1.index
             y3 = pd.DataFrame({'disease': [disease] * X1.shape[0]}) #y3 represents disease of each trained cell.
             y3.index = X1.index
             y4 = pd.DataFrame({'celltype': [celltype] * X1.shape[0]}) #y4 represents celltype of each trained cell.
             y4.index = X1.index
             if not X is None:
                 X = X.append(X1)
                 y_1 = y_1.append(y1)
                 y_2 = y_2.append(y2)
                 y_3 = y_3.append(y3)
                 y_4 = y_4.append(y4)
             else:
                 X = X1
                 y_1 = y1
                 y_2 = y2
                 y_3 = y3
                 y_4 = y4
     #One-hot encoding of cell type labels:
     enc = OneHotEncoder(handle_unknown='ignore')
     enc.fit(pd.concat([y_2, y_3, y_4], axis=1))
     onehot = pd.DataFrame(enc.transform(pd.concat([y_2, y_3, y_4], axis=1)).toarray(), index = y_2.index)
     onehot.columns = enc.get_feature_names_out()
     os.chdir(SPIDER_model_file_path)
     #with open('onehot_retrain_' + protein_name, "wb") as f:
     #    pickle.dump(enc, f)
     X = pd.concat([X, onehot], axis=1)
     #
     y = y_1
     X=X.sample(frac=1,random_state=4905) #shuffle
     X_test=X.sample(n=round((X.shape[0]/10)),random_state=4905) #random 10% holdout validation
     X = X.drop(X_test.index)
     y_test=y.loc[X_test.index, :]
     y_test_index = y_test.index
     y = y.loc[X.index, :]
     #pd.DataFrame(y.index, columns = ['train_cell_id']).to_csv(protein_name + '_cell_features_training_combined_6_sets_random_0.9_internal_training_cell_id_20230115.csv') #Save training data's cell id, so as to compare with baseline models.
     X = torch.tensor(np.vstack(X.values).astype(np.float))
     X = X.type(torch.FloatTensor)#torch.Size([8000, 178])
     y = torch.tensor(y.values)
     y = y.type(torch.FloatTensor)#torch.Size([8000, 1])
     X_test = torch.tensor(np.vstack(X_test.values).astype(np.float))#torch.Size([1000, 178])
     X_test = X_test.type(torch.FloatTensor)
     y_test = torch.tensor(y_test.values)
     y_test = y_test.type(torch.FloatTensor)#torch.Size([1000, 1])
     n_batches = 32 #?
     protein_list=['protein']
     def set_seed(seed):
            np.random.seed(seed)
            random.seed(seed)
            torch.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)  # gpu
            torch.backends.cudnn.deterministic = True
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
     criterion = nn.MSELoss()
     optimizer = optim.Adam(net.parameters(), lr=0.0001,amsgrad=True, weight_decay=0.001)
     max_epochs=200
     train_loss=pd.DataFrame(np.zeros(shape=(len(protein_list),max_epochs)),index=protein_list)
     test_loss=pd.DataFrame(np.zeros(shape=(len(protein_list),max_epochs)),index=protein_list)
     patience=30
     best_score=None
     Dy=len(protein_list)
     estop_counter=pd.Series(np.zeros(Dy),index=protein_list)
     early_stop=pd.Series([False]*Dy,index=protein_list)
     os.chdir(SPIDER_model_file_path)
     torch.cuda.synchronize() # wait for warm-up to finish
     times = []
     for epoch in range(max_epochs):
        if all(early_stop):
            break
        start_epoch = time.time()
        running_loss=pd.Series(np.zeros(Dy),index=protein_list)
        for i in range(int(y.shape[0]/n_batches)):
            # Local batches and labels
            local_X, local_y = X[i*n_batches:min((i+1)*n_batches,X.shape[0]-1),], y[i*n_batches:min((i+1)*n_batches,y.shape[0]-1),]
            # zero the parameter gradients
            optimizer.zero_grad()
            # forward + backward + optimize
            outputs_dict = net(local_X)
            loss=None
            loss_count=0.0
            p = protein_list[0]
            notNaN=(local_y[:,protein_list.index(p):(protein_list.index(p)+1)]==local_y[:,protein_list.index(p):(protein_list.index(p)+1)])
            loss_p=criterion(outputs_dict[p][notNaN],local_y[:,protein_list.index(p):(protein_list.index(p)+1)][notNaN])
            if not torch.isnan(loss_p):
                loss_count+=1.0
                running_loss[p]+=loss_p.item()
                if loss is None:
                    loss=loss_p
                else:
                    loss=loss+loss_p
            loss.backward()
            optimizer.step()
            if(i==(int(y.shape[0]/n_batches)-1)):
                train_loss.iloc[:,epoch]=(running_loss / 150)
            if i % 150 == 149:    # print every mini-batches
                print('[%d, %5d] loss: %.3f' % (epoch + 1, i + 1, sum(running_loss / 150)))
                running_loss=pd.Series(np.zeros(Dy),index=protein_list)
                sys.stdout.flush()
        test_outputs = net(X_test)
        test_outputs = [test_outputs[p] for p in protein_list]
        test_outputs=torch.transpose(torch.stack(test_outputs),0,1).view(X_test.shape[0],-1)
        test_loss_i=pd.Series([criterion(test_outputs[:,pi][y_test[:,pi]==y_test[:,pi]], y_test[:,pi][y_test[:,pi]==y_test[:,pi]]).item() for pi in range(Dy)],index=protein_list)
        test_loss.iloc[:,epoch]=test_loss_i
        # Implement early stopping
        if best_score is None:
            best_score=test_loss_i
        else:
            p = protein_list[0]
            if test_loss_i[p]>(best_score[p]-0.001) and (not early_stop[p]):
                estop_counter[p]+=1
                if estop_counter[p]>=patience:
                    early_stop[p]=True
            else:
                best_score[p]=test_loss_i[p]
                estop_counter[p]=0
            print(estop_counter)
        torch.cuda.synchronize()
        end_epoch = time.time()
        elapsed = end_epoch - start_epoch
        times.append(elapsed)
     print('Finished Training')
     df = pd.DataFrame({"y_pred":test_outputs.detach().numpy().flatten(),'y_truth':y_test.detach().numpy().flatten()})
     df.index = y_test_index
     #df.to_csv(protein_name + '_retrain_internal_val.csv')
     all_internal_val_performance.append(df.corr().iloc[0,1]) #record internal val performance.
     train_loss.index=['train_'+p for p in protein_list]
     test_loss.index=['test_'+p for p in protein_list]
     all_epochs.append(epoch)
    all_epochs = pd.DataFrame({'protein_name': all_protein_list, 'training_epoch': all_epochs})
    all_epochs.to_csv('retrain_record_epochs.csv')
    all_internal_val_performance = pd.DataFrame({'pearson': all_internal_val_performance}, index = all_protein_list)
    all_internal_val_performance.to_csv('retrain_internal_val_performance.csv')
    #-----------------------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------------------
    #train on all samples:
    training_epoch = pd.read_csv(SPIDER_model_file_path + 'retrain_record_epochs.csv', index_col = ['protein_name'])
    for protein_name in all_protein_list:
     select_file_A_B_C_matching_table = file_A_B_C_matching_table
     X = None
     y = None
     y_1 = None #y_1 represents protein abundance in each trained cell.
     y_2 = None #y_2 represents tissue of each trained cell.
     y_3 = None #y_3 represents disease of each trained cell.
     y_4 = None #y_4 represents celltype of each trained cell.
     for i in range(0, len(select_file_A_B_C_matching_table)):
        gene_exp_file = pd.read_csv(select_file_A_B_C_matching_table.iloc[i, :]['file_B'], index_col = 'Unnamed: 0')
        ADT = pd.read_csv(select_file_A_B_C_matching_table.iloc[i, :]['file_C'], index_col = 'Unnamed: 0')
        tissue = select_file_A_B_C_matching_table.iloc[i, :]['tissue_class']
        disease = select_file_A_B_C_matching_table.iloc[i, :]['disease_class']
        celltype = select_file_A_B_C_matching_table.iloc[i, :]['cell_type_class'] 
        if protein_name in ADT.index:
             ADT = pd.DataFrame(ADT.loc[protein_name, :]).T
             protein_list = [protein_name]
             X1 = gene_exp_file
             y1 = pd.DataFrame(np.array(ADT).reshape((1, -1)))
             y1.columns = X1.index
             y1 = y1.T #y1 represents protein abundance in each trained cell.
             y2 = pd.DataFrame({'tissue': [tissue] * X1.shape[0]}) #y2 represents tissue of each trained cell.
             y2.index = X1.index
             y3 = pd.DataFrame({'disease': [disease] * X1.shape[0]}) #y3 represents disease of each trained cell.
             y3.index = X1.index
             y4 = pd.DataFrame({'celltype': [celltype] * X1.shape[0]}) #y4 represents celltype of each trained cell.
             y4.index = X1.index
             if not X is None:
                 X = X.append(X1)
                 y_1 = y_1.append(y1)
                 y_2 = y_2.append(y2)
                 y_3 = y_3.append(y3)
                 y_4 = y_4.append(y4)
             else:
                 X = X1
                 y_1 = y1
                 y_2 = y2
                 y_3 = y3
                 y_4 = y4
     #One-hot encoding of cell type labels:
     enc = OneHotEncoder(handle_unknown='ignore')
     enc.fit(pd.concat([y_2, y_3, y_4], axis=1))
     onehot = pd.DataFrame(enc.transform(pd.concat([y_2, y_3, y_4], axis=1)).toarray(), index = y_2.index)
     onehot.columns = enc.get_feature_names_out()
     os.chdir(SPIDER_model_file_path)
     with open("onehot_retrain_" + protein_name, "wb") as f:
         pickle.dump(enc, f)
     X = pd.concat([X, onehot], axis=1)
     #
     y = y_1
     X=X.sample(frac=1,random_state=4905) #shuffle
     y = y.loc[X.index, :]
     #pd.DataFrame(y.index, columns = ['train_cell_id']).to_csv(protein_name + '_cell_features_training_combined_6_sets_all_training_samples_cell_id_20230115.csv') #Save training data's cell id, so as to compare with baseline models.
     X = torch.tensor(np.vstack(X.values).astype(np.float))
     X = X.type(torch.FloatTensor)#torch.Size([8000, 178])
     y = torch.tensor(y.values)
     y = y.type(torch.FloatTensor)#torch.Size([8000, 1])
     n_batches = 32 #?
     protein_list=['protein']
     def set_seed(seed):
            np.random.seed(seed)
            random.seed(seed)
            torch.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)  # gpu
            torch.backends.cudnn.deterministic = True
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
     criterion = nn.MSELoss()
     optimizer = optim.Adam(net.parameters(), lr=0.0001,amsgrad=True, weight_decay=0.001)
     #max_epochs=200
     max_epochs = training_epoch.loc[protein_name, 'training_epoch'] + 1
     train_loss=pd.DataFrame(np.zeros(shape=(len(protein_list),max_epochs)),index=protein_list)
     test_loss=pd.DataFrame(np.zeros(shape=(len(protein_list),max_epochs)),index=protein_list)
     patience=30
     best_score=None
     Dy=len(protein_list)
     estop_counter=pd.Series(np.zeros(Dy),index=protein_list)
     early_stop=pd.Series([False]*Dy,index=protein_list)
     os.chdir(SPIDER_model_file_path)
     torch.cuda.synchronize() # wait for warm-up to finish
     times = []
     for epoch in range(max_epochs):
        if all(early_stop):
            break
        start_epoch = time.time()
        running_loss=pd.Series(np.zeros(Dy),index=protein_list)
        for i in range(int(y.shape[0]/n_batches)):
            # Local batches and labels
            local_X, local_y = X[i*n_batches:min((i+1)*n_batches,X.shape[0]-1),], y[i*n_batches:min((i+1)*n_batches,y.shape[0]-1),]
            # zero the parameter gradients
            optimizer.zero_grad()
            # forward + backward + optimize
            outputs_dict = net(local_X)
            loss=None
            loss_count=0.0
            p = protein_list[0]
            notNaN=(local_y[:,protein_list.index(p):(protein_list.index(p)+1)]==local_y[:,protein_list.index(p):(protein_list.index(p)+1)])
            loss_p=criterion(outputs_dict[p][notNaN],local_y[:,protein_list.index(p):(protein_list.index(p)+1)][notNaN])
            if not torch.isnan(loss_p):
                loss_count+=1.0
                running_loss[p]+=loss_p.item()
                if loss is None:
                    loss=loss_p
                else:
                    loss=loss+loss_p
            loss.backward()
            optimizer.step()
            if(i==(int(y.shape[0]/n_batches)-1)):
                train_loss.iloc[:,epoch]=(running_loss / 150)
            if i % 150 == 149:    # print every mini-batches
                print('[%d, %5d] loss: %.3f' % (epoch + 1, i + 1, sum(running_loss / 150)))
                running_loss=pd.Series(np.zeros(Dy),index=protein_list)
                sys.stdout.flush()
        torch.cuda.synchronize()
        end_epoch = time.time()
        elapsed = end_epoch - start_epoch
        times.append(elapsed)
     print('Finished Training')
     torch.save(net.state_dict(), protein_name + '_model_ep' + str(epoch))
     train_loss.index=['train_'+p for p in protein_list]
     log=train_loss
     log.to_csv(protein_name + '_log.csv')
     #pd.DataFrame(times).to_csv(protein_name + '_all_training_samples_cell_features_training_combined_6_sets_64_32_16_0.0001_onehot_celltype_tissue_disease_DNN_SCANVI_128dim_record_model_time_20230115.csv')

    
     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    







        






