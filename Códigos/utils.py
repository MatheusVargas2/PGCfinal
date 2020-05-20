import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import time
import pickle
from statistics import stdev
from tqdm.notebook import tqdm_notebook
from itertools import combinations
from sklearn.model_selection import KFold
from sklearn.tree import DecisionTreeClassifier  
from sklearn.svm import SVC 
from sklearn.metrics import accuracy_score


def var_to_pickle(var, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(var, f)
    return 0    

def pickle_to_var(name):
    try:
        with open(name + '.pkl', 'rb') as f:
            var = pickle.load(f)
        return var
    except:
        return 0


class Cache_ngrafo:
    def __init__(self, data, function, n=2):
        self.content = {}
        self.data = data
        self.function = function
        self.n = n
            
    def get(self, var):             
        aux = self.content.get(var)
        if not aux:
            aux = self.function(self.data[self.data.Sample_Label == var].drop(["Sample_Label", "User_Label"], axis=1).iloc[0, :], self.n)
            self.content[var] = aux
        return aux
    
    
class Cache_distacia:
    def __init__(self, function):
        self.content = {}
        self.function = function
            
    def get(self, var, index, t=0):
        
        index.sort()
        concat_index = ""
        for i in index:
            concat_index = concat_index + " / " + str(i)
        concat_index = concat_index[3:]
            
        aux = self.content.get(concat_index)
        if not aux:
            if t==0:
                aux = self.function(var[0], var[1])
            else:
                aux = self.function(var[0], var[1], t)
            self.content[concat_index] = aux
        return aux
    
    
def k_usuarios(id_list, n_folds=5):
    
    import numpy as np

    interno = []
    externo = []
    
    id_list = np.unique(id_list)
    np.random.shuffle(id_list)
    n = len(id_list)
    
    for i in range(n_folds):
        index_externo = [int(n/n_folds * i), int(n/n_folds * i + n/n_folds)]
        interno.append(np.concatenate([id_list[0 : index_externo[0]], id_list[index_externo[1] : n]]))
        externo.append(id_list[index_externo[0] : index_externo[1]])
          
    return interno, externo

def k_amostras(interno, externo, dados, qtd_amostras_usuario=4, n_folds=4, method="LOO"):
    
    treino = []
    autentico = []
    imp_interno = []
    treino_aux = []
    imp_externo = []
    count = 0
    
    if method == "LOO":
    
        # itera entre k_fold de usuários internos x externos
        for index_interno in interno:

            # itera entre cada usuário u
            for u in index_interno:

                samples_u = dados["Sample_Label"][dados["User_Label"] == u]    
                # itera entre as amostras do usuário u
                for j in samples_u:

                    autentico.append(j)

                    aux = list(samples_u.copy())
                    aux.remove(j)
                    treino.append(aux)



                # itera entre cada amostra dos usuários u2 (impostores internos)
                index_interno_2 = list(index_interno.copy())
                index_interno_2.remove(u)
                samples_u2 = dados["Sample_Label"][dados["User_Label"].isin(index_interno_2)]
                samples_u2 = np.array(samples_u2).reshape((int(len(samples_u2)/qtd_amostras_usuario), qtd_amostras_usuario))

                for i in range(0, qtd_amostras_usuario):
                    aux = []
                    for j in samples_u2:
                        aux.append(j[i])
                    imp_interno.append(aux)

                    # treino impostores_internos
                    aux = []
                    for j in samples_u2:
                        aux_2 = []
                        index_treino_imp_interno = list(range(0, qtd_amostras_usuario))
                        index_treino_imp_interno.remove(i)
                        for i2 in index_treino_imp_interno:
                            aux_2.append(j[i2])
                        aux.append(aux_2)
                    treino_aux.append(aux)


                    # obtem index das amostras dos impostores externos
                    index_ext = list(dados["Sample_Label"][dados["User_Label"].isin(externo[count])])
                    imp_externo.append(index_ext)

            count += 1



    elif method == "DifferentTexts":

        # itera entre k_fold de usuários internos x externos
        for index_interno in interno:

            # itera entre cada usuário u
            for u in index_interno:

                samples_u = dados["Sample_Label"][dados["User_Label"] == u]    
                # itera entre as amostras do usuário u
                for j in samples_u:

                    # Seleciona amostra autentica de teste
                    autentico.append(j)
                    
                    aux = list(samples_u.copy())
                    aux.remove(j)
                    
                    # Remove do treino amostras que contém o mesmo texto que a amostra de validação (TrueEssay ou FakeEssay)
                    if ("True" in j) or ("Copy_1" in j):
                        for txt in aux:
                            if ("True" in txt) or ("Copy_1" in txt):
                                aux.remove(txt)
                    else:
                        for txt in aux:
                            if ("Fake" in txt) or ("Copy_2" in txt):
                                aux.remove(txt)
                    
                    treino.append(aux)

                    
                # itera entre cada amostra dos usuários u2 (impostores internos)
                index_interno_2 = list(index_interno.copy())
                index_interno_2.remove(u)
                samples_u2 = dados["Sample_Label"][dados["User_Label"].isin(index_interno_2)]
                samples_u2 = np.array(samples_u2).reshape((int(len(samples_u2)/qtd_amostras_usuario), qtd_amostras_usuario))

                for i in range(0, qtd_amostras_usuario):
                    aux = []
                    for j in samples_u2:
                        aux.append(j[i])
                    imp_interno.append(aux)

                    # treino impostores_internos
                    aux = []
                    for j in samples_u2:
                        aux_2 = []
                        index_treino_imp_interno = list(range(0, qtd_amostras_usuario))
                        index_treino_imp_interno.remove(i)
                        for i2 in index_treino_imp_interno:
                            aux_2.append(j[i2])
                        aux.append(aux_2)
                    treino_aux.append(aux)


                    # obtem index das amostras dos impostores externos
                    index_ext = list(dados["Sample_Label"][dados["User_Label"].isin(externo[count])])
                    imp_externo.append(index_ext)

            count += 1
            
            
    
    elif method == "DifferentPatterns":

        # itera entre k_fold de usuários internos x externos
        for index_interno in interno:

            # itera entre cada usuário u
            for u in index_interno:

                samples_u = dados["Sample_Label"][dados["User_Label"] == u]    
                # itera entre as amostras do usuário u
                for j in samples_u:

                    # Seleciona amostra autentica de teste
                    autentico.append(j)
                    
                    aux = list(samples_u.copy())
                    aux.remove(j)
                    
                    # Remove do treino amostras que contém o mesmo padrão que a amostra de validação (Original ou Copy)
                    if ("True" in j) or ("Fake" in j):
                        for txt in aux:
                            if ("True" in txt) or ("Fake" in txt):
                                aux.remove(txt)
                    else:
                        for txt in aux:
                            if ("Copy_1" in txt) or ("Copy_2" in txt):
                                aux.remove(txt)
                    
                    treino.append(aux)

                    
                # itera entre cada amostra dos usuários u2 (impostores internos)
                index_interno_2 = list(index_interno.copy())
                index_interno_2.remove(u)
                samples_u2 = dados["Sample_Label"][dados["User_Label"].isin(index_interno_2)]
                samples_u2 = np.array(samples_u2).reshape((int(len(samples_u2)/qtd_amostras_usuario), qtd_amostras_usuario))

                for i in range(0, qtd_amostras_usuario):
                    aux = []
                    for j in samples_u2:
                        aux.append(j[i])
                    imp_interno.append(aux)

                    # treino impostores_internos
                    aux = []
                    for j in samples_u2:
                        aux_2 = []
                        index_treino_imp_interno = list(range(0, qtd_amostras_usuario))
                        index_treino_imp_interno.remove(i)
                        for i2 in index_treino_imp_interno:
                            aux_2.append(j[i2])
                        aux.append(aux_2)
                    treino_aux.append(aux)


                    # obtem index das amostras dos impostores externos
                    index_ext = list(dados["Sample_Label"][dados["User_Label"].isin(externo[count])])
                    imp_externo.append(index_ext)

            count += 1
    
    else:
        print("Método", method, "não disponível.")
        return 0
        
    return [treino, autentico, imp_interno, imp_externo, treino_aux]