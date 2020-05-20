from utils import *

def extrai_duracao_ngrafo(vetor,n):
    features = {}

    i = 0
    while not vetor[i] == "":
        if vetor[i].split()[1] == "KeyDown":
            count = 1
            ngrafo = vetor[i].split()[2]    
            t0 = int(vetor[i].split()[0])
            
            if n == 1:
                if features.get(ngrafo):
                    features[ngrafo].append(t0)
                else:
                    features[ngrafo] = [t0]
                    
                    
            else:
                j = i + 1
                while j > 0 and not vetor[j] == "":

                    if vetor[j].split()[1] == "KeyDown":
                        count += 1
                        ngrafo = ngrafo + " " + vetor[j].split()[2]

                    if count >= n:                        
                        k = j + 1
                        while k > 0 and not vetor[k] == "":
                            if vetor[k].split()[2] == vetor[j].split()[2]:
                                t1 = int(vetor[k].split()[0])
                                delta_tempo = t1 - t0
                                j = -1
                                k = -1
                                if features.get(ngrafo):
                                    features[ngrafo].append(delta_tempo)
                                else:
                                    features[ngrafo] = [delta_tempo]
                            k+=1
                    j+=1
        i+=1    
    
    vetor_ngrafos = []
    for i in features:
        vetor_ngrafos.append([sum(features[i])/len(features[i]), i])
    
    vetor_ngrafos.sort()
    return vetor_ngrafos


def extrai_features_digrafo(vetor):
    H1 = {}
    H2 = {}
    DD = {}
    UU = {}
    UD = {}

    i = 0
    while not vetor[i] == "":

        if vetor[i].split()[1] == "KeyDown":

            key_1 = vetor[i].split()[2]
            t1_down = int(vetor[i].split()[0])

            key_1_released = False
            key_2_pressed = False
            key_2_released = False
            flag_erro = False
            key_2 = "nan"

            j = i + 1
            while j > 0 and not vetor[j] == "":

                if vetor[j].split()[1] == "KeyUp" and vetor[j].split()[2] == key_1:
                    t1_up = int(vetor[j].split()[0])
                    key_1_released = True
                elif vetor[j].split()[1] == "KeyDown":
                    key_2 = vetor[j].split()[2] 
                    t2_down = int(vetor[j].split()[0])
                    key_2_pressed = True
                elif vetor[j].split()[1] == "KeyUp" and vetor[j].split()[2] == key_2 and key_2_pressed:
                    t2_up = int(vetor[j].split()[0])
                    key_2_released = True

                if key_1_released and key_2_pressed and key_2_released:
                    j = -1

                # limite de busca
                if j > i+10:
                    j = -1
                    flag_erro = True


                j+=1

            if not flag_erro and not key_2 == "nan":
                digrafo = key_1 + " " + key_2

                if H1.get(digrafo):
                    H1[digrafo].append(t1_up - t1_down)
                else:
                    H1[digrafo] = [t1_up - t1_down]

                if H2.get(digrafo):
                    H2[digrafo].append(t2_up - t2_down)
                else:
                    H2[digrafo] = [t2_up - t2_down]

                if DD.get(digrafo):
                    DD[digrafo].append(t2_down - t1_down)
                else:
                    DD[digrafo] = [t2_down - t1_down]

                if UU.get(digrafo):
                    UU[digrafo].append(t2_up - t1_up)
                else:
                    UU[digrafo] = [t2_up - t1_up]

                if UD.get(digrafo):
                    UD[digrafo].append(t2_down - t1_up)
                else:
                    UD[digrafo] = [t2_down - t1_up]

        i+=1    

    features_digrafos = []
    for i in H1:
        features_digrafos.append([i, sum(H1[i])/len(H1[i]), sum(H2[i])/len(H2[i]), sum(DD[i])/len(DD[i]), sum(UU[i])/len(UU[i]), sum(UD[i])/len(UD[i])])

    
    return features_digrafos



def ngrafos_comuns(V1, V2):
    
    ngrafos = []
    V1_aux = []
    V2_aux = []
    
    for i in range(len(V1)):
        for j in range(len(V2)):
            if V1[i][1] == V2[j][1]:
                V1_aux.append((V1[i][0], V1[i][1]))
                V2_aux.append((V2[j][0], V2[j][1]))
                break
    V2_aux.sort()

    return V1_aux, V2_aux


def distancia_R(V1, V2):
    
    V1, V2 = ngrafos_comuns(V1, V2)
    
    dist = 0
    for i in range(len(V1)):
        for j in range(len(V2)):
            if V1[i][1] == V2[j][1]:
                dist = dist + abs(i-j)
                break
              
    if len(V1) % 2 == 0:
        d_m = (len(V1) ** 2) / 2
    else:
        d_m = (len(V1) ** 2 - 1) / 2
                
    return dist / d_m



def distancia_A(V1, V2, t):     
    
    tamanho = 0
    similares = 0
    for i in range(len(V1)):
        for j in range(len(V2)):
            if V1[i][1] == V2[j][1]:
                tamanho += 1
                eq = max(V1[i][0],V2[j][0]) / min(V1[i][0],V2[j][0])
                if eq >= 1 and eq < t:
                    similares += 1
                    break
    
    return 1 - (similares / tamanho)

    


def wpm(vetor, qtd_words=50):
    
    # define tempo inicial
    i = 0
    while not vetor[i] == "":            
        if vetor[i].split()[1] == "KeyDown":
            t0 = int(vetor[i].split()[0])
            break          
        i+=1    
    

    i = 0
    words = 0
    medias_wpm = []
    while not vetor[i] == "":            

        if vetor[i].split()[1] == "KeyDown" and vetor[i].split()[2] == "32":
            words += 1
            
            if words == qtd_words:
                delta_t = (int(vetor[i].split()[0]) - t0) / 60000
                medias_wpm.append(words/delta_t)
                
                t0 = int(vetor[i].split()[0])
                words = 0
                
        i+=1    
    
    if len(medias_wpm) != 0:
        return sum(medias_wpm) / len(medias_wpm)
    else:
        return 0


def neg_UD_UU(vetor):
    
    i = 0
    count = 0
    neg_ud = 0
    neg_uu = 0
    medias_wpm = []
    while not vetor[i] == "":            

        if vetor[i].split()[1] == "KeyDown" and vetor[i+1] != "" and vetor[i+2] != "":
            count += 1
            
            if vetor[i+1].split()[1] == "KeyDown" and vetor[i+2].split()[2] == vetor[i].split()[2]:
                neg_ud += 1
            elif vetor[i+1].split()[1] == "KeyDown" and vetor[i+2].split()[2] == vetor[i+1].split()[2]:
                neg_uu += 1
                
        i+=1    
    
    if count != 0:
        return neg_ud/count, neg_uu/count
    else:
        return 0, 0
    




def error_rate(vetor):
    
    i = 0
    count = 0
    error = 0
    while not vetor[i] == "":            

        if vetor[i].split()[1] == "KeyDown":
            count += 1
            
            if vetor[i].split()[2] == "8" or vetor[i].split()[2] == "46":
                error += 1
                count -= 1
                
        i+=1    
    
    if count != 0:
        return error/count
    else:
        return 0


def caps_usage(vetor):
    
    i = 0
    count = 0
    caps = 0
    while not vetor[i] == "":            

        if vetor[i].split()[1] == "KeyDown":
            count += 1
            
            if vetor[i].split()[2] == "20":
                caps += 1
                
        i+=1    
    
    if count != 0:
        return caps/count
    else:
        return 0


def shift_after_before(vetor):
    
    i = 0
    count = 0
    after = 0
    before = 0
    vazio = 0
    
    while not vetor[i] == "":            

        if vetor[i].split()[1] == "KeyDown" and vetor[i].split()[2] == "16" and vetor[i+1] != "" and vetor[i+2] != "":
            count += 1
            
            if vetor[i+2].split()[1] == "KeyUp" and vetor[i+2].split()[2] != "16":
                after += 1
            elif vetor[i+2].split()[1] == "KeyUp" and vetor[i+2].split()[2] == "16":
                before += 1
            else:
                vazio += 1
                
        i+=1  
    
    if count != 0:
        return after/count, before/count, vazio/count
    else:
        return 0, 0, 0



def nonConventionalFeatures(vetor):
    wpm_ = wpm(vetor)
    neg_UD_, neg_UU_ = neg_UD_UU(vetor)
    error_rate_ = error_rate(vetor)
    caps_usage_ = caps_usage(vetor)
    shift_after_, shift_before_, _ = shift_after_before(vetor)
    return [wpm_, neg_UD_, neg_UU_, error_rate_, caps_usage_, shift_after_, shift_before_]



def features_adjacent_keys(vetor, adjacencia_teclas):
    features = []
    for i in vetor:
        keys = i[0].split()
        if int(keys[0]) in list(adjacencia_teclas.Key_Code) and int(keys[1]) in list(adjacencia_teclas.Key_Code):

            if keys[1] in list(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[0])].Adjacent)[0]:
                adj = "Adjacent"
            elif keys[1] in list(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[0])].SecondAdjacent)[0]:
                adj = "SecondAdjacent"
            elif keys[1] in list(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[0])].ThirdAdjacent)[0]:
                adj = "ThirdAdjacent"
            elif keys[1] in list(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[0])].FourthAdjacent)[0]:
                adj = "FourthAdjacent"
            elif keys[1] in list(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[0])].NonAdjacent)[0]:
                adj = "NonAdjacent"

            if int(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[0])].Side == "Right") & int(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[1])].Side == "Right"):
                side = "Right"
            elif int(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[0])].Side == "Left") & int(adjacencia_teclas[adjacencia_teclas.Key_Code == int(keys[1])].Side == "Left"):
                side = "Left"
            else:
                side = "Different"

            features.append([i[0], adj, side, i[1], i[2], i[3], i[4], i[5]])

    H1 = {}
    H2 = {}
    DD = {}
    UU = {}
    UD = {}
    for i in features:

        grupo = i[1]+i[2]

        if H1.get(grupo):
            H1[grupo].append(i[3])
        else:
            H1[grupo] = [i[3]]

        if H2.get(grupo):
            H2[grupo].append(i[4])
        else:
            H2[grupo] = [i[4]]

        if DD.get(grupo):
            DD[grupo].append(i[5])
        else:
            DD[grupo] = [i[5]]

        if UU.get(grupo):
            UU[grupo].append(i[6])
        else:
            UU[grupo] = [i[6]]

        if UD.get(grupo):
            UD[grupo].append(i[7])
        else:
            UD[grupo] = [i[7]]

    media_features = {}
    for i in H1:
        media_features[i + "_H1"] = sum(H1[i])/len(H1[i])
        media_features[i + "_H2"] = sum(H2[i])/len(H2[i])
        media_features[i + "_DD"] = sum(DD[i])/len(DD[i])
        media_features[i + "_UU"] = sum(UU[i])/len(UU[i])
        media_features[i + "_UD"] = sum(UD[i])/len(UD[i])
        
    return media_features




















def modelo(mod, dados, kf, params):
    
    dados_x = dados.drop(["Sample_Label", "User_Label"], axis=1)
    dados_y = dados[["Sample_Label", "User_Label"]]
    
    FNM = []
    FM_Int = []
    FM_Ext = []
    
    
    # Mods 1 e 2
    salva_cache = {}
    k_list1 = []
    k_list2 = []
    k_list3 = []
    flag_list1 = []
    flag_list2 = []
    flag_list3 = []
    
    if mod in [1, 2, 2.1]:
        
        # Inicializa Caches
        cache_ngrafos = {}
        cache_distancias_R = {}
        cache_distancias_A = {}

        #Verifica se existe cache dos ngrafos utilizados para calcular as distâncias R e A (se não existe, cria)
        for n in np.unique(params["dist_R"] + params["dist_A"]):
            cache_ngrafos[n] = pickle_to_var("../Caches/" + params["nome_cache"] + str(n) + "_n")
            if not cache_ngrafos[n]:
                cache_ngrafos[n] = Cache_ngrafo(dados, extrai_duracao_ngrafo, n)
                salva_cache["../Caches/" + params["nome_cache"] + str(n) + "_n"] = cache_ngrafos[n]
                print("Novo cache criado:", params["nome_cache"] + str(n) + "_n")
            else:
                cache_ngrafos[n].data = dados
                
        #Verifica se existe cache das distâncias R (se não existe, calcula e cria)
        for n in params["dist_R"]:
            cache_distancias_R[n] = pickle_to_var("../Caches/" + params["nome_cache"] + str(n) + "_R_n")
            if not cache_distancias_R[n]:
                cache_distancias_R[n] = Cache_distacia(distancia_R)
                salva_cache["../Caches/" + params["nome_cache"] + str(n) + "_R_n"] = cache_distancias_R[n]
                print("Novo cache criado:", params["nome_cache"] + str(n) + "_R_n")
         
        #Verifica se existe cache das distâncias A com respectivo threshold t (se não existe, calcula e cria)
        t_str = str(round(params["t"]*100))[1:3]
        for n in params["dist_A"]:
            cache_distancias_A[n] = pickle_to_var("../Caches/" + params["nome_cache"] + str(n) + "_A_p" + t_str + "_n")
            if not cache_distancias_A[n]:
                cache_distancias_A[n] = Cache_distacia(distancia_A)
                salva_cache["../Caches/" + params["nome_cache"] + str(n) + "_A_p" + t_str + "_n"] = cache_distancias_A[n]
                print("Novo cache criado:", params["nome_cache"] + str(n) + "_A_p" + t_str + "_n")
                
        # Lista de ngrafos utilzados
        n_list = list(cache_ngrafos.keys())
    

    
    
    # Mod 3
    if mod == 3:  
        if params["action"] == "new":
            # Extrai features de todos os usuários
            features = []
            j=0
            for i in tqdm_notebook(dados_x.itertuples()):
                features.append(nonConventionalFeatures(i[1:]) + [dados_y["Sample_Label"][j]] + [dados_y["User_Label"][j]])
                j+=1
            features = pd.DataFrame(features, columns = ["wpm", "neg_UD", "neg_UU","error_rate", "caps_usage", "shift_after", "shift_before", "Sample_Label", "User_Label"])
            var_to_pickle(features, "../Caches/features_mod3")
        elif params["action"] == "load":
            features = pickle_to_var("../Caches/features_mod3")
            
            
            
    # Mod 4
    if mod == 4:
        adjacencia_teclas = pd.read_csv("../Outros/adjacencia_teclas.csv")
        if params["action"] == "new":
            # Extrai features de todos os usuários
            features = []
            j=0
            for i in tqdm_notebook(dados_x.itertuples()):
                features_digrafos = extrai_features_digrafo(i[1:])
                adjacent_features = features_adjacent_keys(features_digrafos, adjacencia_teclas)
                
                non_conv_features = nonConventionalFeatures(i[1:])
                adjacent_features["wpm"] = non_conv_features[0]
                adjacent_features["neg_UD"] = non_conv_features[1]
                adjacent_features["neg_UU"] = non_conv_features[2]
                adjacent_features["error_rate"] = non_conv_features[3]
                adjacent_features["caps_usage"] = non_conv_features[4]
                adjacent_features["shift_after"] = non_conv_features[5]
                adjacent_features["shift_before"] = non_conv_features[6]
                adjacent_features["Sample_Label"] = dados_y["Sample_Label"][j]
                adjacent_features["User_Label"] = dados_y["User_Label"][j]
      
                features.append(adjacent_features)
                j+=1
                
            features = pd.DataFrame(features)
            features = features.fillna(0)
            var_to_pickle(features, "../Caches/features_mod4")
        elif params["action"] == "load":
            features = pickle_to_var("../Caches/features_mod4")
            
        
        
######################################################################################################################    

    for i in tqdm_notebook(range(len(kf[0]))):
        
        index_train = kf[0][i]
        index_genuino = kf[1][i] 
        index_imp_interno = kf[2][i] 
        index_imp_externo = kf[3][i] 
        index_train_imp = kf[4][i]
        
        if mod == 1:
            
            #### Cria modelo m
            m = []
            for j in combinations(index_train, 2):
            
                # Recupera (ou salva) vetor de ngrafos no dicionário
                v1 = {}
                v2 = {}
                for n in n_list:
                    v1[n] = cache_ngrafos[n].get(j[0])
                    v2[n] = cache_ngrafos[n].get(j[1])

                # Calcula distância entre v1 e v2
                dist = 0
                for n in params["dist_R"]:
                    dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [j[0], j[1]])
                  
                m.append(dist)
                
            m = sum(m)/len(m)         

            
            
            #### Valida usuário legítimo v1 (com cada amostra de treino v2)
            md = []
            v1 = {}
            v2 = {}
            for n in n_list:
                v1[n] = cache_ngrafos[n].get(index_genuino)
            for a in index_train:
                for n in n_list:
                    v2[n] = cache_ngrafos[n].get(a)

                # Calcula distância entre v1 e v2
                dist = 0
                for n in params["dist_R"]:
                    dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [index_genuino, a])
                md.append(dist)
                
            k_list1.append((sum(md)/len(md))/m)
                
            if sum(md)/len(md) < m * params["k"]:
                FNM.append(0)
            else:
                FNM.append(1)


            
            #### Valida usuários impostores internos v1 (com cada amostra de treino v2) 
            for j in index_imp_interno:
                
                md = []
                v1 = {}
                v2 = {}
                for n in n_list:
                    v1[n] = cache_ngrafos[n].get(j)
             
                for a in index_train:
                    for n in n_list:
                        v2[n] = cache_ngrafos[n].get(a)
            
                    # Calcula distância entre v1 e v2
                    dist = 0
                    for n in params["dist_R"]:
                        dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [j, a])
                    md.append(dist)
                
                k_list2.append((sum(md)/len(md))/m)
                
                if sum(md)/len(md) < m * params["k"]:
                    FM_Int.append(1)
                else:
                    FM_Int.append(0)
            


            

            #### Valida usuários impostores externos v1 (com cada amostra de treino v2) 
            for j in index_imp_externo:
                
                md = []
                v1 = {}
                v2 = {}
                for n in n_list:
                    v1[n] = cache_ngrafos[n].get(j)
             
                for a in index_train:
                    for n in n_list:
                        v2[n] = cache_ngrafos[n].get(a)
            
                    # Calcula distância entre v1 e v2
                    dist = 0
                    for n in params["dist_R"]:
                        dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [j, a])
                    md.append(dist)
                
                k_list3.append((sum(md)/len(md))/m)
                
                if sum(md)/len(md) < m * params["k"]:
                    FM_Ext.append(1)
                else:
                    FM_Ext.append(0)
            
        
        
        
######################################################################################################################


        elif mod == 2:
            
            #### Cria modelo m
            m = []
            for j in combinations(index_train, 2):
            
                # Recupera (ou salva) vetor de ngrafos no dicionário
                v1 = {}
                v2 = {}
                for n in n_list:
                    v1[n] = cache_ngrafos[n].get(j[0])
                    v2[n] = cache_ngrafos[n].get(j[1])

                # Calcula distância entre v1 e v2
                dist = 0
                for n in params["dist_R"]:
                    dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [j[0], j[1]])
                for n in params["dist_A"]:
                    dist = dist + cache_distancias_A[n].get([v1[n], v2[n]], [j[0], j[1]], params["t"])

                m.append(dist)
                
            m = sum(m)/len(m)

            
            #### Calcula md(A,x) entre amostras de todos os usuários A cadastrados (legítimo e imp internos) e 
            #### tentativas de autenticação x (legítimo, imp interno e imp externo)
            md = {}
            for x in ([index_genuino] + index_imp_interno + index_imp_externo):
                
                md[x] = {}
                md_aux = []
                v1 = {}
                v2 = {}
                for n in n_list:
                    v1[n] = cache_ngrafos[n].get(x)
                for A in ([index_train] + index_train_imp):
                    for a in A:
                        for n in n_list:
                            v2[n] = cache_ngrafos[n].get(a)
                        
                        # Calcula distância entre v1 e v2
                        dist = 0
                        for n in params["dist_R"]:
                            dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [x, a])
                        for n in params["dist_A"]:
                            dist = dist + cache_distancias_A[n].get([v1[n], v2[n]], [x, a], params["t"])
                            
                        md_aux.append(dist)
                        
                    md[x][str(A)] = sum(md_aux)/len(md_aux)
                
       
                
            ##### Valida usuário legítimo   
            # Define o menor e segundo menor md
            md_values = list(md[index_genuino].values())
            md_values.sort()
            first_min = md_values[0]
            second_min = md_values[1]

            for A, md_A in md[index_genuino].items():
                if md_A == first_min:
                    index_first_min = A
                    break
            
            # 1) Verifica se md(A, x) é o menor
            flag_menor = 0
            if md[index_genuino][str(index_train)] == first_min:
                flag_menor = 1
            else:
                second_min = first_min
            
            # 2) Verifica se md(A, x) é próximo à m(A)
            flag_proximo_m = 0
            if md[index_genuino][str(index_train)] < m + params["k"]*(second_min - m):
                flag_proximo_m = 1
            
            # Validação das condições
            k_list1.append((md[index_genuino][str(index_train)] - m) / (second_min - m))
            flag_list1.append([flag_menor, flag_proximo_m])

            if flag_menor == 1 and flag_proximo_m == 1:
                FNM.append(0)
            else:
                FNM.append(1)
                
               
                
                
            ##### Valida impostor interno (iteração entre cada um)
            for j in index_imp_interno:
                # Define o menor e segundo menor md
                md_values = list(md[j].values())
                md_values.sort()
                first_min = md_values[0]
                second_min = md_values[1]

                for A, md_A in md[j].items():
                    if md_A == first_min:
                        index_first_min = A
                        break

                # 1) Verifica se md(A, x) é o menor
                flag_menor = 0
                if md[j][str(index_train)] == first_min:
                    flag_menor = 1
                else:
                    second_min = first_min

                # 2) Verifica se md(A, x) é próximo à m(A)
                flag_proximo_m = 0
                if md[j][str(index_train)] < m + params["k"]*(second_min - m):
                    flag_proximo_m = 1
                    
                    
                # Validação das condições
                k_list2.append((md[j][str(index_train)] - m) / (second_min - m))
                flag_list2.append([flag_menor, flag_proximo_m])

                if flag_menor == 1 and flag_proximo_m == 1:
                    FM_Int.append(1)
                else:
                    FM_Int.append(0)
                    
                    
               
            ##### Valida impostor externo (iteração entre cada um)
            for j in index_imp_externo:
                # Define o menor e segundo menor md
                md_values = list(md[j].values())
                md_values.sort()
                first_min = md_values[0]
                second_min = md_values[1]

                for A, md_A in md[j].items():
                    if md_A == first_min:
                        index_first_min = A
                        break

                # 1) Verifica se md(A, x) é o menor
                flag_menor = 0
                if md[j][str(index_train)] == first_min:
                    flag_menor = 1
                else:
                    second_min = first_min

                # 2) Verifica se md(A, x) é próximo à m(A)
                flag_proximo_m = 0
                if md[j][str(index_train)] < m + params["k"]*(second_min - m):
                    flag_proximo_m = 1
                    
                
                # Validação das condições
                k_list3.append((md[j][str(index_train)] - m) / (second_min - m))
                flag_list3.append([flag_menor, flag_proximo_m])

                if flag_menor == 1 and flag_proximo_m == 1:
                    FM_Ext.append(1)
                else:
                    FM_Ext.append(0)
            

        
######################################################################################################################
        
        
        
        
        elif mod == 2.1:
            
            #### Cria modelo m
            m = []
            for j in combinations(index_train, 2):
            
                # Recupera (ou salva) vetor de ngrafos no dicionário
                v1 = {}
                v2 = {}
                for n in n_list:
                    v1[n] = cache_ngrafos[n].get(j[0])
                    v2[n] = cache_ngrafos[n].get(j[1])

                # Calcula distância entre v1 e v2
                dist = 0
                for n in params["dist_R"]:
                    dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [j[0], j[1]])
                for n in params["dist_A"]:
                    dist = dist + cache_distancias_A[n].get([v1[n], v2[n]], [j[0], j[1]], params["t"])

                  
                m.append(dist)
                
            m = sum(m)/len(m)         

            
            
            #### Valida usuário legítimo v1 (com cada amostra de treino v2)
            md = []
            v1 = {}
            v2 = {}
            for n in n_list:
                v1[n] = cache_ngrafos[n].get(index_genuino)
            for a in index_train:
                for n in n_list:
                    v2[n] = cache_ngrafos[n].get(a)

                # Calcula distância entre v1 e v2
                dist = 0
                for n in params["dist_R"]:
                    dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [index_genuino, a])
                for n in params["dist_A"]:
                    dist = dist + cache_distancias_A[n].get([v1[n], v2[n]], [index_genuino, a], params["t"])
                md.append(dist)
                
            k_list1.append((sum(md)/len(md))/m)
                
            if sum(md)/len(md) < m * params["k"]:
                FNM.append(0)
            else:
                FNM.append(1)


            
            #### Valida usuários impostores internos v1 (com cada amostra de treino v2) 
            for j in index_imp_interno:
                
                md = []
                v1 = {}
                v2 = {}
                for n in n_list:
                    v1[n] = cache_ngrafos[n].get(j)
             
                for a in index_train:
                    for n in n_list:
                        v2[n] = cache_ngrafos[n].get(a)
            
                    # Calcula distância entre v1 e v2
                    dist = 0
                    for n in params["dist_R"]:
                        dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [j, a])
                    for n in params["dist_A"]:
                        dist = dist + cache_distancias_A[n].get([v1[n], v2[n]], [j, a], params["t"])
                    md.append(dist)
                
                k_list2.append((sum(md)/len(md))/m)
                
                if sum(md)/len(md) < m * params["k"]:
                    FM_Int.append(1)
                else:
                    FM_Int.append(0)
            


            

            #### Valida usuários impostores externos v1 (com cada amostra de treino v2) 
            for j in index_imp_externo:
                
                md = []
                v1 = {}
                v2 = {}
                for n in n_list:
                    v1[n] = cache_ngrafos[n].get(j)
             
                for a in index_train:
                    for n in n_list:
                        v2[n] = cache_ngrafos[n].get(a)
            
                    # Calcula distância entre v1 e v2
                    dist = 0
                    for n in params["dist_R"]:
                        dist = dist + cache_distancias_R[n].get([v1[n], v2[n]], [j, a])
                    for n in params["dist_A"]:
                        dist = dist + cache_distancias_A[n].get([v1[n], v2[n]], [j, a], params["t"])
                    md.append(dist)
                
                k_list3.append((sum(md)/len(md))/m)
                
                if sum(md)/len(md) < m * params["k"]:
                    FM_Ext.append(1)
                else:
                    FM_Ext.append(0)

            
            
            
            
######################################################################################################################


            
        elif mod == 3:
               
            ##### Define conjuntos de treino/validação
            # Converte label_sample -> label_user
            index_train_full = []
            index_train_full_ = []
            index_train_ = []
            index_train_imp_ = []
            index_genuino_ = index_genuino.split('|')[1]
            index_imp_interno_ = []
            index_imp_externo_ = []
            
            for j in index_train:
                index_train_.append(j.split('|')[1])
                index_train_full_.append(j.split('|')[1])
                index_train_full.append(j)
                
            for j in index_train_imp:
                for k in j:
                    index_train_imp_.append(k.split('|')[1])
                    index_train_full_.append(k.split('|')[1])
                    index_train_full.append(k)

            for j in index_imp_interno:
                index_imp_interno_.append(j.split('|')[1])
                index_imp_interno_.append(j.split('|')[1])
            
            for j in index_imp_externo:
                index_imp_externo_.append(j.split('|')[1])
                index_imp_externo_.append(j.split('|')[1])
                           
            # Features X
            x_train = features[features["Sample_Label"].isin(index_train_full)].drop(["Sample_Label", "User_Label"], axis=1).reset_index(drop=True)
            x_legit = features[features["Sample_Label"] == index_genuino].drop(["Sample_Label", "User_Label"], axis=1).reset_index(drop=True)  #.values.reshape(1, -1)
            x_imp_int = features[features["Sample_Label"].isin(index_imp_interno)].drop(["Sample_Label", "User_Label"], axis=1).reset_index(drop=True)
            x_imp_ext = features[features["Sample_Label"].isin(index_imp_externo)].drop(["Sample_Label", "User_Label"], axis=1).reset_index(drop=True)
            
            # Label Multiclasse Y
            y_train_mc = features[features["Sample_Label"].isin(index_train_full)]["User_Label"].reset_index(drop=True)
            y_legit_mc = features[features["Sample_Label"] == index_genuino]["User_Label"].reset_index(drop=True)
            y_imp_int_mc = features[features["Sample_Label"].isin(index_imp_interno)]["User_Label"].reset_index(drop=True)
            y_imp_ext_mc = features[features["Sample_Label"].isin(index_imp_externo)]["User_Label"].reset_index(drop=True)
            
            # Label Binário Y
            y_train_bin = y_train_mc.copy()
            y_train_bin[y_train_bin.isin(index_train_)] = 1
            y_train_bin[y_train_bin.isin(index_train_imp_)] = 0
            y_train_bin = list(y_train_bin)
            
            y_legit_bin = 1
            
            y_imp_int_bin = y_imp_int_mc.copy()
            y_imp_int_bin[:] = 0
            
            y_imp_ext_bin = y_imp_ext_mc.copy()
            y_imp_ext_bin[:] = 0

            
            ##### Treina Modelos
            if params["algorithm"] == "tree":
                
                # Árvore de decisão multiclasse
                mod_mc = DecisionTreeClassifier(random_state=0)
                mod_mc.fit(x_train, y_train_mc)

                # Árvore de decisão binária
                mod_bin = DecisionTreeClassifier(random_state=0)
                mod_bin.fit(x_train, y_train_bin)
                
            elif params["algorithm"] == "svm":
                
                # SVM multiclasse
                mod_mc = SVC(random_state=0, gamma = "auto", probability=True)
                mod_mc.fit(x_train, y_train_mc)

                # SVM binário
                mod_bin = SVC(random_state=0, gamma = "auto", probability=True)
                mod_bin.fit(x_train, y_train_bin)
                
                
                
            ##### Validação
            # Usuários legítimos
            pred_legit_mc = mod_mc.predict(x_legit)
            pred_legit_bin = mod_bin.predict(x_legit) 
            
            if params["decision"] == "binary":
                if pred_legit_bin[0] == y_legit_bin:
                    FNM = np.concatenate((FNM, [0]))
                else:
                    FNM = np.concatenate((FNM, [1]))
                    
            elif params["decision"] == "multiclass":
                if pred_legit_mc[0] == y_legit_mc[0]:
                    FNM = np.concatenate((FNM, [0]))
                else:
                    FNM = np.concatenate((FNM, [1]))
                
            elif params["decision"] == "merge":
                if pred_legit_mc[0] == y_legit_mc[0] or pred_legit_bin[0] == y_legit_bin:
                    FNM = np.concatenate((FNM, [0]))
                else:
                    FNM = np.concatenate((FNM, [1]))
                


            # Impostores internos
            pred_imp_int_mc = mod_mc.predict(x_imp_int)
            pred_imp_int_bin = mod_bin.predict(x_imp_int)
            
            erros = []
            for j in range(len(pred_imp_int_mc)):

                if params["decision"] == "binary":
                    if pred_imp_int_bin[j] == y_legit_bin:
                        erros.append(1)
                    else:
                        erros.append(0)
                        
                elif params["decision"] == "multiclass":
                    if pred_imp_int_mc[j] == y_legit_mc[0]:
                        erros.append(1)
                    else:
                        erros.append(0)
                
                elif params["decision"] == "merge":
                    if pred_imp_int_mc[j] == y_legit_mc[0] or pred_imp_int_bin[j] == y_legit_bin:
                        erros.append(1)
                    else:
                        erros.append(0)
                    
            FM_Int = np.concatenate((FM_Int, erros))

            
            
            # Impostores externos
            pred_imp_ext_mc = mod_mc.predict(x_imp_ext)
            pred_imp_ext_bin = mod_bin.predict(x_imp_ext)
            
            erros = []
            for j in range(len(pred_imp_ext_mc)):
                
                if params["decision"] == "binary":
                    
                    if pred_imp_ext_bin[j] == y_legit_bin:
                        erros.append(1)
                    else:
                        erros.append(0)
                        
                elif params["decision"] == "multiclass":
                    if pred_imp_ext_mc[j] == y_legit_mc[0]:
                        erros.append(1)
                    else:
                        erros.append(0)
                    
                elif params["decision"] == "merge":
                    if pred_imp_ext_mc[j] == y_legit_mc[0] or pred_imp_ext_bin[j] == y_legit_bin:
                        erros.append(1)
                    else:
                        erros.append(0)
                    
            FM_Ext = np.concatenate((FM_Ext, erros))
                
            
            
######################################################################################################################
        
        elif mod == 4:
                           
            ##### Define conjuntos de treino/validação
            # Converte label_sample -> label_user
            index_train_full = []
            index_train_full_ = []
            index_train_ = []
            index_train_imp_ = []
            index_genuino_ = index_genuino.split('|')[1]
            index_imp_interno_ = []
            index_imp_externo_ = []
            
            for j in index_train:
                index_train_.append(j.split('|')[1])
                index_train_full_.append(j.split('|')[1])
                index_train_full.append(j)
                
            for j in index_train_imp:
                for k in j:
                    index_train_imp_.append(k.split('|')[1])
                    index_train_full_.append(k.split('|')[1])
                    index_train_full.append(k)

            for j in index_imp_interno:
                index_imp_interno_.append(j.split('|')[1])
                index_imp_interno_.append(j.split('|')[1])
            
            for j in index_imp_externo:
                index_imp_externo_.append(j.split('|')[1])
                index_imp_externo_.append(j.split('|')[1])

                           
            # Features X
            x_train = features[features["Sample_Label"].isin(index_train_full)].drop(["Sample_Label", "User_Label"], axis=1).reset_index(drop=True)            
            x_legit = features[features["Sample_Label"] == index_genuino].drop(["Sample_Label", "User_Label"], axis=1).reset_index(drop=True)  #.values.reshape(1, -1)
            x_imp_int = features[features["Sample_Label"].isin(index_imp_interno)].drop(["Sample_Label", "User_Label"], axis=1).reset_index(drop=True)
            x_imp_ext = features[features["Sample_Label"].isin(index_imp_externo)].drop(["Sample_Label", "User_Label"], axis=1).reset_index(drop=True)
            
            # Label Multiclasse Y
            y_train_mc = features[features["Sample_Label"].isin(index_train_full)]["User_Label"].reset_index(drop=True)
            y_legit_mc = features[features["Sample_Label"] == index_genuino]["User_Label"].reset_index(drop=True)
            y_imp_int_mc = features[features["Sample_Label"].isin(index_imp_interno)]["User_Label"].reset_index(drop=True)
            y_imp_ext_mc = features[features["Sample_Label"].isin(index_imp_externo)]["User_Label"].reset_index(drop=True)
            
            # Label Binário Y
            y_train_bin = y_train_mc.copy()
            y_train_bin[y_train_bin.isin(index_train_)] = 1
            y_train_bin[y_train_bin.isin(index_train_imp_)] = 0
            y_train_bin = list(y_train_bin)
            
            y_legit_bin = 1
            
            y_imp_int_bin = y_imp_int_mc.copy()
            y_imp_int_bin[:] = 0
            
            y_imp_ext_bin = y_imp_ext_mc.copy()
            y_imp_ext_bin[:] = 0
            

            
            ##### Treina Modelos
            if params["algorithm"] == "tree":
                
                # Árvore de decisão multiclasse
                mod_mc = DecisionTreeClassifier(   criterion=params["model_params"]["criterion"],
                                                   splitter=params["model_params"]["splitter"],
                                                   max_depth=params["model_params"]["max_depth"],
                                                   min_samples_split=params["model_params"]["min_samples_split"],
                                                   min_samples_leaf=params["model_params"]["min_samples_leaf"],
                                                   min_weight_fraction_leaf=params["model_params"]["min_weight_fraction_leaf"],
                                                   max_features=params["model_params"]["max_features"],
                                                   random_state=params["model_params"]["random_state"],
                                                   max_leaf_nodes=params["model_params"]["max_leaf_nodes"],
                                                   min_impurity_decrease=params["model_params"]["min_impurity_decrease"],
                                                   min_impurity_split=params["model_params"]["min_impurity_split"],
                                                   class_weight=params["model_params"]["class_weight"])
                mod_mc.fit(x_train, y_train_mc)

                # Árvore de decisão binária
                mod_bin = DecisionTreeClassifier(  criterion=params["model_params"]["criterion"],
                                                   splitter=params["model_params"]["splitter"],
                                                   max_depth=params["model_params"]["max_depth"],
                                                   min_samples_split=params["model_params"]["min_samples_split"],
                                                   min_samples_leaf=params["model_params"]["min_samples_leaf"],
                                                   min_weight_fraction_leaf=params["model_params"]["min_weight_fraction_leaf"],
                                                   max_features=params["model_params"]["max_features"],
                                                   random_state=params["model_params"]["random_state"],
                                                   max_leaf_nodes=params["model_params"]["max_leaf_nodes"],
                                                   min_impurity_decrease=params["model_params"]["min_impurity_decrease"],
                                                   min_impurity_split=params["model_params"]["min_impurity_split"],
                                                   class_weight=params["model_params"]["class_weight"])
                mod_bin.fit(x_train, y_train_bin)
                
            elif params["algorithm"] == "svm":
                
                # SVM multiclasse
                mod_mc = SVC(random_state=0, gamma="auto", probability=True)
                mod_mc.fit(x_train, y_train_mc)

                # SVM binário
                mod_bin = SVC(random_state=0, gamma="auto", probability=True)
                mod_bin.fit(x_train, y_train_bin)
                
                
                
            ##### Validação
            # Usuários legítimos
            pred_legit_mc = mod_mc.predict(x_legit)
            pred_legit_bin = mod_bin.predict(x_legit)

            if params["decision"] == "binary":
                if pred_legit_bin[0] == y_legit_bin:
                    FNM = np.concatenate((FNM, [0]))
                else:
                    FNM = np.concatenate((FNM, [1]))
                    
            elif params["decision"] == "multiclass":
                if pred_legit_mc[0] == y_legit_mc[0]:
                    FNM = np.concatenate((FNM, [0]))
                else:
                    FNM = np.concatenate((FNM, [1]))
                
            elif params["decision"] == "merge":
                if pred_legit_mc[0] == y_legit_mc[0] or pred_legit_bin[0] == y_legit_bin:
                    FNM = np.concatenate((FNM, [0]))
                else:
                    FNM = np.concatenate((FNM, [1]))
                


            # Impostores internos
            pred_imp_int_mc = mod_mc.predict(x_imp_int)
            pred_imp_int_bin = mod_bin.predict(x_imp_int)
            
            erros = []
            for j in range(len(pred_imp_int_mc)):
                if params["decision"] == "binary":
                    if pred_imp_int_bin[j] == y_legit_bin:
                        erros.append(1)
                    else:
                        erros.append(0)
                        
                elif params["decision"] == "multiclass":
                    if pred_imp_int_mc[j] == y_legit_mc[0]:
                        erros.append(1)
                    else:
                        erros.append(0)
                
                elif params["decision"] == "merge":
                    if pred_imp_int_mc[j] == y_legit_mc[0] or pred_imp_int_bin[j] == y_legit_bin:
                        erros.append(1)
                    else:
                        erros.append(0)
                    
            FM_Int = np.concatenate((FM_Int, erros))

            
            
            # Impostores externos
            pred_imp_ext_mc = mod_mc.predict(x_imp_ext)
            pred_imp_ext_bin = mod_bin.predict(x_imp_ext)
            
            erros = []
            for j in range(len(pred_imp_ext_mc)):
                
                if params["decision"] == "binary":
                    
                    if pred_imp_ext_bin[j] == y_legit_bin:
                        erros.append(1)
                    else:
                        erros.append(0)
                        
                elif params["decision"] == "multiclass":
                    if pred_imp_ext_mc[j] == y_legit_mc[0]:
                        erros.append(1)
                    else:
                        erros.append(0)
                    
                elif params["decision"] == "merge":
                    if pred_imp_ext_mc[j] == y_legit_mc[0] or pred_imp_ext_bin[j] == y_legit_bin:
                        erros.append(1)
                    else:
                        erros.append(0)
                    
            FM_Ext = np.concatenate((FM_Ext, erros))

            
        elif mod == 5:
            pass
        

    
    for key, item in salva_cache.items():
        var_to_pickle(item, key)
    
    FMR_Int = sum(FM_Int)/len(FM_Int)
    FMR_Ext = sum(FM_Ext)/len(FM_Ext)
    FMR = (sum(FM_Int) + sum(FM_Ext)) / (len(FM_Int) + len(FM_Ext))
    FNMR = sum(FNM)/len(FNM)
    BAcc = 1 - (FMR + FNMR)/2
    
    return [FNMR, FMR, FMR_Int, FMR_Ext, BAcc, k_list1, k_list2, k_list3, flag_list1, flag_list2, flag_list3]
