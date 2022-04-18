from gekko import GEKKO
from resources.msc_rafael_pavan import *
import copy

def balanco_potencia_ativa(sep,pg_sgen,pg, pc, barras_origem, barras_destino, barra_atual, gkm_linhas, bkm_linhas, tensoes, angulos, to, td, tap, bkmt,gkmt,pshunt):
    
    soma = 0
    
    fluxos = []
    
    
    linhas = np.arange(0,len(barras_origem),1)
    
    baux = []
    baux.append(33333333)
    baux.append(33333331)
    
    for bd in barras_destino[barras_origem==barra_atual]:
        
        baux.append(bd)
        posi = linhas[(barras_destino==bd) & (barras_origem==barra_atual)][0]
        
        if baux[-2]==bd:
            
            posi = linhas[(barras_destino==bd) & (barras_origem==barra_atual)][1]
            
        soma = soma + gkm_linhas[posi]*(tensoes[barra_atual]**2) - tensoes[barra_atual]*tensoes[bd]*(gkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(gkm_linhas[posi]*(tensoes[barra_atual]**2) - tensoes[barra_atual]*tensoes[bd]*(gkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    for bd in barras_origem[barras_destino==barra_atual]:
        
        posi = linhas[(barras_destino==barra_atual) & (barras_origem==bd)][0]
        
        baux.append(bd)
        
        if baux[-2]==bd:

            posi = linhas[(barras_destino==barra_atual) & (barras_origem==bd)][1]
        
        soma = soma + gkm_linhas[posi]*(tensoes[barra_atual]**2) - tensoes[barra_atual]*tensoes[bd]*(gkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(gkm_linhas[posi]*(tensoes[barra_atual]**2) - tensoes[barra_atual]*tensoes[bd]*(gkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
    
    linhas = np.arange(0,len(to),1)
    
    for bd in td[to==barra_atual]:
        
        posi = linhas[(td==bd) & (to==barra_atual)][0]

        soma = soma + (gkmt[posi]*tensoes[barra_atual]*tensoes[barra_atual]/(tap[posi]/100)**2 - (tensoes[barra_atual]*tensoes[bd]/(tap[posi]/100))*(gkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
        fluxos.append(gkmt[posi]*tensoes[barra_atual]*tensoes[barra_atual]/(tap[posi]/100)**2 - (tensoes[barra_atual]*tensoes[bd]/(tap[posi]/100))*(gkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
    
    for bd in to[td==barra_atual]:
        
        posi = linhas[(td==barra_atual) & (to==bd)][0]

        soma = soma + (gkmt[posi]*tensoes[barra_atual]*tensoes[barra_atual]  - (tensoes[barra_atual]*tensoes[bd]/(tap[posi]/100))*(gkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
    
        fluxos.append(gkmt[posi]*tensoes[barra_atual]*tensoes[barra_atual] - (tensoes[barra_atual]*tensoes[bd]/(tap[posi]/100))*(gkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))

    return  pg[barra_atual] - pc[barra_atual] - soma + pg_sgen[barra_atual] - pshunt[barra_atual]*tensoes[barra_atual]**2



def balanco_potencia_reativa(sep,qg, qc, barras_origem, barras_destino, barra_atual, gkm_linhas, bkm_linhas, tensoes, angulos, to, td, tap, bkmt,gkmt, bshl, bsht,qshunt):
    
    soma = 0
    
    fluxos = []
    
    linhas = np.arange(0,len(barras_origem),1)
    
    baux = []
    baux.append(10101010101)
    baux.append(1000101010)
    
    for bd in barras_destino[barras_origem==barra_atual]:
        
        baux.append(bd)
        posi = linhas[(barras_destino==bd) & (barras_origem==barra_atual)][0]
        
        if baux[-2]==bd:
            
            posi = linhas[(barras_destino==bd) & (barras_origem==barra_atual)][1]
            
        
        soma = soma + -(-bkm_linhas[posi]+bshl[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*tensoes[bd]*(-bkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(-(-bkm_linhas[posi]+bshl[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*tensoes[bd]*(-bkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    for bd in barras_origem[barras_destino==barra_atual]:
        
        posi = linhas[(barras_destino==barra_atual) & (barras_origem==bd)][0]
        
        baux.append(bd)
        
        if baux[-2]==bd:

            posi = linhas[(barras_destino==barra_atual) & (barras_origem==bd)][1]
        
        soma = soma + -(-bkm_linhas[posi]+bshl[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*tensoes[bd]*(-bkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(-(-bkm_linhas[posi]+bshl[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*tensoes[bd]*(-bkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    linhas = np.arange(0,len(to),1)
    
    for bd in td[to==barra_atual]:
        
        posi = linhas[(td==bd) & (to==barra_atual)][0]         

        soma = soma + -(-bkmt[posi]/((tap[posi]/100)**2)+bsht[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*(1/(tap[posi]/100))*tensoes[bd]*(-bkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(-(-bkmt[posi]/((tap[posi]/100)**2)+bsht[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*(1/(tap[posi]/100))*tensoes[bd]*(-bkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    for bd in to[td==barra_atual]:
        
        posi = linhas[(td==barra_atual) & (to==bd)][0]

        soma = soma + -(-bkmt[posi]+bsht[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*(1/(tap[posi]/100))*tensoes[bd]*(-bkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(-(-bkmt[posi]+bsht[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*(1/(tap[posi]/100))*tensoes[bd]*(-bkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    return  qg[barra_atual] - qc[barra_atual] - soma - (qshunt[barra_atual]/100)*tensoes[barra_atual]**2


# In[4]:


def perdas(sep,gkml, gkmt, angulos, tensoes, tap, origem, destino, hv, lv):
    
    i = 0
    
    eq = []
   
    for bus in zip(origem,destino):
        
        
        perdas = gkml[i]*(tensoes[bus[0]]**2 + tensoes[bus[1]]**2 - 2*tensoes[bus[1]]*tensoes[bus[0]]*sep.cos(angulos[bus[0]]-angulos[bus[1]]))
        i=i+1
        
        eq.append(perdas)
    
    j = 0
    

    
    for bus in zip(hv,lv):
        
        perdas = gkmt[j]*((tensoes[bus[0]]/(tap[j]/100))**2 + tensoes[bus[1]]**2 - 2*tensoes[bus[1]]*tensoes[bus[0]]*(1/(tap[j]/100))*sep.cos(angulos[bus[0]]-angulos[bus[1]]))
        
        eq.append(perdas)
            
        j=j+1
        
    return perdas, eq
    


# In[5]:


def voltage_dev(tensoes):
    
    
    sum_dev = 0
    
    for tensao in tensoes:
        dev = (1-tensao)**2

        sum_dev = dev + sum_dev
        
    return sum_dev



def Sequential_Quadratic_Programming_Branch_and_Bound(sep_teste, n_vizinhos_tap, n_vizinhos_shunt, verbose=True, travado=False):
    
    
    sep_modelo = copy.copy(sep_teste) 
    sep_modelo.res_line = sep_modelo.res_line.sort_index()
    sep_modelo.line = sep_modelo.line.sort_index()
    origem = sep_modelo.line[['from_bus']].values
    destino = sep_modelo.line[['to_bus']].values
    
        ########################################################################### Vetores de condutância e susceptância série

    m_z = np.zeros((5,len(sep_modelo.line)))

    gkm = np.zeros(len(sep_modelo.line))

    bkm = np.zeros(len(sep_modelo.line))

    bo = np.zeros(len(sep_modelo.line))

    bd = np.zeros(len(sep_modelo.line))

    sep_modelo.line = sep_modelo.line.sort_index()

    sep_modelo.bus = sep_modelo.bus.sort_index()

    vbus = sep_modelo.bus.vn_kv.to_numpy(dtype=np.float64)

    zbase = np.power(np.multiply(vbus,1000), 2)/(100*1e6)

    m_z[0,:] = sep_modelo.line.from_bus.to_numpy()

    m_z[1,:] = sep_modelo.line.to_bus.to_numpy()

    bsh = 1e-9*(2*np.pi*60*sep_modelo.line.c_nf_per_km.to_numpy())

    m_z[4,:] = bsh


    for i in range(len(sep_modelo.line.index.ravel())):    

        m_z[2,i] = sep_modelo.line.r_ohm_per_km[i]/zbase[int(m_z[0,i])]

        m_z[3,i] = sep_modelo.line.x_ohm_per_km[i]/zbase[int(m_z[0,i])]

        m_z[4,i] =  m_z[4,i] * zbase[int(m_z[0,i])]


    gkm = np.array(np.divide(m_z[2,:], np.power(m_z[2,:],2)+np.power(m_z[3],2)))

    bo = m_z[0,:]

    bd = m_z[1,:]

    ########################################################################### Vetor de susceptância

    bkm = np.array(np.divide(m_z[3,:], np.power(m_z[2,:],2)+np.power(m_z[3],2)))


    ########################################################################### Vetor de susceptância shunt

    bsh = m_z[4,:]
    

    
    sep_modelo.trafo = sep_modelo.trafo.sort_index()

    barras = sep_modelo.trafo['hv_bus'].to_numpy()
    
    zkm = (sep_modelo.trafo['vk_percent'].to_numpy()/100)*(1000/sep_modelo.trafo['sn_mva'].to_numpy())
    
    rkm = (sep_modelo.trafo['vkr_percent'].to_numpy()/100)*(1000/sep_modelo.trafo['sn_mva'].to_numpy())
    
    #a = (sep_modelo.trafo['vn_lv_kv'].to_numpy()*sep_modelo.trafo['vn_lv_kv'].to_numpy()*1000/sep_modelo.trafo['sn_mva'].to_numpy())/(sep_modelo.trafo['vn_lv_kv'].to_numpy()*sep_modelo.trafo['vn_lv_kv'].to_numpy()/1000)
    
    a = 1
    
    zkm=zkm/10
    
    rkm=rkm/10
    
    xkm = np.sqrt(zkm**2-rkm**2)
    

#     xkm[91] = 0.0231
        
    gkmt = (rkm*a/((a*rkm)**2+(a*xkm)**2))
    
    bkmt = (xkm*a/((a*rkm)**2+(a*xkm)**2))
    
    bsht = np.sqrt(np.power(sep_modelo.trafo['i0_percent'].to_numpy()/100,2))
    
    bsht = bsht*99
    

        ########################################################################### Vetor de tap

    tap_pos = sep_modelo.trafo[~pd.isnull(sep_modelo.trafo['tap_pos'])]['tap_pos'].to_numpy(dtype=np.float64)

    tap_neutral = sep_modelo.trafo[~pd.isnull(sep_modelo.trafo['tap_neutral'])]['tap_neutral'].to_numpy(dtype=np.float64)

    tap_step_percent = sep_modelo.trafo[~pd.isnull(sep_modelo.trafo['tap_step_percent'])]['tap_step_percent'].to_numpy(dtype=np.float64)

    valor_percentual = (tap_pos-tap_neutral)*(tap_step_percent/100) + 1

#     valor_percentual = np.resize(valor_percentual,(len(sep_modelo.trafo)))


    to = np.zeros(len(sep_modelo.trafo))
    td = np.zeros(len(sep_modelo.trafo))
    
    for i in range(len(sep_modelo.trafo)):
        
        if sep_modelo.trafo['tap_side'].iloc[i] == None or sep_modelo.trafo['tap_side'].iloc[i] == 'hv':
        
            to[i] = int(sep_modelo.trafo['hv_bus'].iloc[i])


            td[i] = int(sep_modelo.trafo['lv_bus'].iloc[i])

        if sep_modelo.trafo['tap_side'].iloc[i] == 'lv':
        
            to[i] = int(sep_modelo.trafo['lv_bus'].iloc[i])


            td[i] = int(sep_modelo.trafo['hv_bus'].iloc[i])
            
    to = to.astype(int)
    td = td.astype(int)

    i = 0

    for i in range(len(valor_percentual)):

        if i < len(tap_pos):

            valor_percentual[i] = valor_percentual[i]

        else:

            valor_percentual[i] = 1


    tap = valor_percentual
    
    sep_modelo.trafo['tap_pos'][~pd.isnull(sep_modelo.trafo['tap_pos'])] = valor_percentual
    sep_modelo.trafo['tap_pos'][pd.isnull(sep_modelo.trafo['tap_pos'])] = 1
    
    tap = sep_modelo.trafo['tap_pos'].values
    
    
    ########################################################################### Vetor de tensões das barras
    sep_modelo.line = sep_modelo.line.sort_index()

    sep_modelo.res_bus = sep_modelo.res_bus.sort_index()
    
    sep_modelo.sgen = sep_modelo.sgen.sort_index()

    v = sep_modelo.res_bus['vm_pu'].to_numpy()
    Sbase=100
    ########################################################################### Vetor de ângulos das barras

    theta = np.radians(sep_modelo.res_bus['va_degree'].to_numpy())

    ########################################################################### Vetor de potência ativa gerada
    sep_modelo.gen = sep_modelo.gen.sort_index()
    
    pg = np.zeros(len(sep_modelo.bus))
    
    pgs_max = np.zeros(len(sep_modelo.bus))
    
    pgs_min = np.zeros(len(sep_modelo.bus))
    
    pg_sgen = np.zeros(len(sep_modelo.bus))
    qg_sgen = np.zeros(len(sep_modelo.bus))
    
    i = 0

    sep_modelo.gen = sep_modelo.gen.sort_index()

    sep_modelo.res_gen = sep_modelo.res_gen.sort_index()

    for bus in sep_modelo.gen['bus'].to_numpy():

        pg[bus] = sep_modelo.gen['p_mw'].to_numpy()[i]/Sbase
        pgs_max[bus] = sep_modelo.gen['max_p_mw'].to_numpy()[i]/Sbase
        pgs_min[bus] = sep_modelo.gen['min_p_mw'].to_numpy()[i]/Sbase
        i=i+1
        

    i=0
    
    qg = np.zeros(len(pg))
   
    for bus in sep_modelo.sgen['bus'].to_numpy():
        
        pg_sgen[bus] = sep_modelo.sgen['p_mw'].to_numpy()[i]/Sbase
        
        qg_sgen[bus] = sep_modelo.sgen['q_mvar'].to_numpy()[i]/Sbase

        i = i+1
    
    if len(sep_modelo.bus)==118:

        pg[68] = sep_modelo.res_ext_grid['p_mw'].to_numpy()/100
        qg[68] = sep_modelo.res_ext_grid['q_mvar'].to_numpy()/100
        pgs_max[68] = sep_modelo.ext_grid['max_p_mw'].to_numpy()/100
        pgs_min[68] = sep_modelo.ext_grid['min_p_mw'].to_numpy()/100
        
    if len(sep_modelo.bus)==300:

        pg[256] = sep_modelo.res_ext_grid['p_mw'].to_numpy()/100
        qg[256] = sep_modelo.res_ext_grid['q_mvar'].to_numpy()/100
        

    
    if len(sep_modelo.bus)==14 or len(sep_modelo.bus)==30 :

        pg[0] = sep_modelo.res_ext_grid['p_mw'].to_numpy()/100
        pgs_max[0] = sep_modelo.ext_grid['max_p_mw'].to_numpy()/100
        pgs_min[0] = sep_modelo.ext_grid['min_p_mw'].to_numpy()/100
        
    pg_ls = sep_modelo.ext_grid['max_p_mw'].to_numpy()/100

    pg_li = sep_modelo.ext_grid['min_p_mw'].to_numpy()/100
    
    sep_modelo.load = sep_modelo.load.sort_index()
    pc = np.zeros(len(sep_modelo.bus))

    i = 0

    sep_modelo.load = sep_modelo.load.sort_index()

    for bus in sep_modelo.load['bus'].to_numpy():

        pc[bus] = sep_modelo.load['p_mw'].to_numpy()[i]/Sbase

        i=i+1

    qc = np.zeros(len(sep_modelo.bus))

    i = 0

    for bus in sep_modelo.load['bus'].to_numpy():

        qc[bus] = sep_modelo.load['q_mvar'].to_numpy()[i]/Sbase

        i=i+1


########################################################################### Vetores de condutância e susceptância série

    m_z = np.zeros((5,len(sep_modelo.line)))

    gkm = np.zeros(len(sep_modelo.line))

    bkm = np.zeros(len(sep_modelo.line))

    bo = np.zeros(len(sep_modelo.line))

    bd = np.zeros(len(sep_modelo.line))

    sep_modelo.line = sep_modelo.line.sort_index()

    sep_modelo.bus = sep_modelo.bus.sort_index()

    vbus = sep_modelo.bus.vn_kv.to_numpy(dtype=np.float64)

    zbase = np.power(np.multiply(vbus,1000), 2)/(100*1e6)

    m_z[0,:] = sep_modelo.line.from_bus.to_numpy()

    m_z[1,:] = sep_modelo.line.to_bus.to_numpy()

    bsh = 1e-9*(2*np.pi*60*sep_modelo.line.c_nf_per_km.to_numpy())

    m_z[4,:] = bsh


    for i in range(len(sep_modelo.line.index.ravel())):    

        m_z[2,i] = sep_modelo.line.r_ohm_per_km[i]/zbase[int(m_z[0,i])]

        m_z[3,i] = sep_modelo.line.x_ohm_per_km[i]/zbase[int(m_z[0,i])]

        m_z[4,i] =  m_z[4,i] * zbase[int(m_z[0,i])]


    gkm = np.array(np.divide(m_z[2,:], np.power(m_z[2,:],2)+np.power(m_z[3],2)))

    bo = m_z[0,:]

    bd = m_z[1,:]

    ########################################################################### Vetor de susceptância

    bkm = np.array(np.divide(m_z[3,:], np.power(m_z[2,:],2)+np.power(m_z[3],2)))


    ########################################################################### Vetor de susceptância shunt

    bsh = m_z[4,:]
    
    qg = np.zeros(len(pg))

    if len(sep_modelo.bus)==14 or len(sep_modelo.bus)==30:

        qg[0] = sep_modelo.res_ext_grid['q_mvar'].values/100

    if len(sep_modelo.bus)==118:

        qg[68] = sep_modelo.res_ext_grid['q_mvar'].values/100


    if len(sep_modelo.bus)==300:

        qg[256] = sep_modelo.res_ext_grid['q_mvar'].values/100


    sepaux = sep_modelo.gen['bus']

    for barra in sepaux:

        qg[barra]=sep_modelo.res_gen[sep_modelo.gen['bus']==barra]['q_mvar'].values/100


    barras = sep_modelo.shunt['bus'].to_numpy()
    qshunt = np.zeros(np.shape(qg))
    pshunt = np.zeros(np.shape(pg))

    for barra in barras:

        qshunt[barra]=sep_modelo.shunt[sep_modelo.shunt['bus']==barra]['q_mvar'].values/100
        pshunt[barra]=sep_modelo.shunt[sep_modelo.shunt['bus']==barra]['p_mw'].values/100

    hv=sep_modelo.trafo['hv_bus'].values
    lv=sep_modelo.trafo['lv_bus'].values

    if len(sep_modelo.bus) == 118:
        


        tensoes = []
        angulos = []

        sep_modelo.res_bus = sep_modelo.res_bus.sort_index()

        v = sep_modelo.res_bus['vm_pu'].to_numpy()

        theta = np.radians(sep_modelo.res_bus['va_degree'].to_numpy())

        sep = GEKKO(remote=True)

        for bus in range(len(sep_modelo.bus)):
        
            if len(sep_modelo.bus)==118:
                
                tensoes.append(sep.Var(v[bus],0.94,1.06))
            if bus>0:
                
                angulos.append(sep.Var(theta[bus],-np.pi,np.pi))
            else: angulos.append(0)
                
        minimum = pgs_min.copy()
        maximum = pgs_max.copy()

        if len(sep_modelo.bus)==118:

            posicoes = [1,4,6,15,34,70]

            for val in posicoes:
                
                val=val-1
                
                if minimum[val]<= pg[val]<=20/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=20/100

                elif 30/100 <= pg[val] <= 60/100:


                    pgs_min[val] = 30/100
                    pgs_max[val] = 60/100
                    
                elif 85/100<= pg[val]<=maximum[val]:


                    pgs_min[val] = 85/100
                    pgs_max[val] = maximum[val].copy()

            posicoes = [10]

            for val in posicoes:
                val=val-1
                if minimum[val]<= pg[val]<=15/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=15/100

                elif 45/100 <= pg[val] <= 165/100:


                    pgs_min[val] = 45/100
                    pgs_max[val] = 165/100
                    

                elif 200/100 <= pg[val] <= 395/100:


                    pgs_min[val] = 200/100
                    pgs_max[val] = 395/100
                    
                    
                elif 410/100<= pg[val]<=maximum[val]:


                    pgs_min[val] = 410/100
                    pgs_max[val] = maximum[val].copy()

            posicoes = [25]
            

            for val in posicoes:
                val=val-1
                if minimum[val]<= pg[val]<=40/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=40/100

                elif 65/100 <= pg[val] <= 190/100:


                    pgs_min[val] = 65/100
                    pgs_max[val] = 190/100
                    
                    
                elif 200/100<= pg[val]<=maximum[val]:


                    pgs_min[val] = 200/100
                    pgs_max[val] = maximum[val].copy()
                    
                    
            posicoes = [26]


            for val in posicoes:
                val=val-1
                if minimum[val]<= pg[val]<=75/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=75/100

                elif 95/100 <= pg[val] <= 260/100:


                    pgs_min[val] = 95/100
                    pgs_max[val] = 260/100
                    
                    
                elif 280/100<= pg[val]<=maximum[val]:


                    pgs_min[val] = 280/100
                    pgs_max[val] = maximum[val].copy()
                    
            posicoes = [49]
            
            for val in posicoes:
                val=val-1
                if minimum[val]<= pg[val]<=45/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=45/100

                elif 60/100 <= pg[val] <= 185/100:


                    pgs_min[val] = 60/100
                    pgs_max[val] = 185/100
                    
                    
                elif 200/100<= pg[val]<=maximum[val]:


                    pgs_min[val] = 200/100
                    pgs_max[val] = maximum[val].copy()
                    
                    
            posicoes = [59]
            

            for val in posicoes:
                val=val-1
                if minimum[val]<= pg[val]<=95/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=95/100

                elif 105/100 <= pg[val] <= 140/100:


                    pgs_min[val] = 105/100
                    pgs_max[val] = 140/100
                    
                    
                elif 155/100<= pg[val]<=maximum[val]:


                    pgs_min[val] = 155/100
                    pgs_max[val] = maximum[val].copy()
                    
            posicoes = [61]

            for val in posicoes:
                
                val=val-1

                if minimum[val]<= pg[val]<=145/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=145/100

                elif 155/100 <= pg[val] <= 210/100:


                    pgs_min[val] = 155/100
                    pgs_max[val] = 210/100
                    
                    
                elif 230/100<= pg[val]<=maximum[val]:


                    pgs_min[val] = 230/100
                    pgs_max[val] = maximum[val].copy()

                    
            posicoes = [65]

            for val in posicoes:
                
                val=val-1
                
                
                if minimum[val]<= pg[val]<=180/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=180/100

                elif 200/100<= pg[val] <= 350/100:


                    pgs_min[val] = 200/100
                    pgs_max[val] = 350/100
                    
                    
                elif 360/100<= pg[val]<=maximum[val]:


                    pgs_min[val] = 360/100
                    pgs_max[val] = maximum[val].copy()

            posicoes = [89]
            
            


            for val in posicoes:
                
                val=val-1

                if minimum[val]<= pg[val]<=120/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=120/100

                elif 145/100 <= pg[val] <= 410/100:


                    pgs_min[val] = 145/100
                    pgs_max[val] = 410/100
                    

                elif 460/100 <= pg[val] <= 500/100:


                    pgs_min[val] = 460/100
                    pgs_max[val] = 500/100
                    
                elif 525/100 <= pg[val] <= maximum[val]:

                    pgs_min[val] = 525/100
                    pgs_max[val] = maximum[val].copy()
                    
            posicoes = [40,42,85,99,104,116]
            
            
            for val in posicoes:
                val=val-1
                if minimum[val]<= pg[val]<=20/100:
                        
                    pgs_min[val]=minimum[val].copy()
                    pgs_max[val]=20/100

                elif 30/100 <= pg[val] <= 45/100:


                    pgs_min[val] = 30/100
                    pgs_max[val] = 45/100
                    

                elif 55/100 <= pg[val] <= maximum[val]:

                    pgs_min[val] = 55/100
                    pgs_max[val] = maximum[val].copy()
                    
        shunt = np.zeros(len(sep_modelo.bus)).tolist()

#         shunt[4]= sep.Var(qshunt[4],0,0.40)
#         shunt[33]= sep.Var(qshunt[33],-0.2,0)
#         shunt[36]= sep.Var(qshunt[36],0,0.25)
#         shunt[43]= sep.Var(qshunt[43],-0.1,0)
#         shunt[44]= sep.Var(qshunt[44],-0.1,0)
#         shunt[45]= sep.Var(qshunt[45],-0.1,0)
#         shunt[47]= sep.Var(qshunt[47],-0.15,0)
#         shunt[73]= sep.Var(qshunt[73],-0.2,0)
#         shunt[78]= sep.Var(qshunt[78],-0.2,0)
#         shunt[81]= sep.Var(qshunt[81],-0.2,0)
#         shunt[82]= sep.Var(qshunt[82],-0.2,0)
#         shunt[104]= sep.Var(qshunt[104],-0.2,0)
#         shunt[106]= sep.Var(qshunt[106],-0.2,0)
#         shunt[109]= sep.Var(qshunt[109],-0.2,0)

        auxx = n_vizinhos_shunt
        aux_ls = shunt[4]-auxx
        aux_us = shunt[4]+auxx
        
        
        if aux_ls < 0:
            
            aux_ls = 0
            
        if aux_us>0.40:
        
            aux_us = 0.40
        
        
        shunt[4]= sep.Var(qshunt[4]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[33]-auxx
        aux_us = shunt[33]+auxx
        
        
        if aux_ls < -0.14:
            
            aux_ls = -0.14
            
        if aux_us>0:
            
            aux_us = 0
        
        shunt[33]= sep.Var(qshunt[33]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[36]-auxx
        aux_us = shunt[36]+auxx
        
        
        if aux_ls < 0:
            
            aux_ls = 0
            
        if aux_us>0.25:
            
            aux_us = 0.25  
            
        shunt[36]= sep.Var(qshunt[36]*100,aux_ls*100,aux_us*100,integer=True)
        
        
        aux_ls = shunt[43]-auxx
        aux_us = shunt[43]+auxx
        
        
        if aux_ls < -0.1:
            
            aux_ls = -0.1
            
        if aux_us>0:
            
            aux_us = 0 
            
        
        shunt[43]= sep.Var(qshunt[43]*100,aux_ls*100,aux_us*100,integer=True)

        aux_ls = shunt[44]-auxx
        aux_us = shunt[44]+auxx
        
        
        if aux_ls < -0.1:
            
            aux_ls = -0.1
            
        if aux_us>0:
            
            aux_us = 0 
        
        shunt[44]= sep.Var(qshunt[44]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[45]-auxx
        aux_us = shunt[45]+auxx
        
        
        if aux_ls < -0.1:
            
            aux_ls = -0.1
            
        if aux_us>0:
            
            aux_us = 0 
            
            
        shunt[45]= sep.Var(qshunt[45]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[47]-auxx
        aux_us = shunt[47]+auxx
        
        
        if aux_ls < -0.15:
            
            aux_ls = -0.15
            
        if aux_us>0:
            
            aux_us = 0 
            
        shunt[47]= sep.Var(qshunt[47]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[73]-auxx
        aux_us = shunt[73]+auxx
        
        
        if aux_ls < -0.12:
            
            aux_ls = -0.12
            
        if aux_us>0:
            
            aux_us = 0 
            
        shunt[73]= sep.Var(qshunt[73]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[78]-auxx
        aux_us = shunt[78]+auxx
        
        
        if aux_ls < -0.2:
            
            aux_ls = -0.2
            
        if aux_us>0:
            
            aux_us = 0 
        
        
        shunt[78]= sep.Var(qshunt[78]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[81]-auxx
        aux_us = shunt[81]+auxx
        
        
        if aux_ls < -0.2:
            
            aux_ls = -0.2
            
        if aux_us>0:
            
            aux_us = 0 
        
        
        shunt[81]= sep.Var(qshunt[81]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[82]-auxx
        aux_us = shunt[82]+auxx
        
        if aux_ls < -0.1:
            
            aux_ls = -0.1
            
        if aux_us>0:
            
            aux_us = 0 
        
        shunt[82]= sep.Var(qshunt[82]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[104]-auxx
        aux_us = shunt[104]+auxx
        
        if aux_ls < -0.2:
            
            aux_ls = -0.2
            
        if aux_us>0:
            
            aux_us = 0 
        
        shunt[104]= sep.Var(qshunt[104]*100,aux_ls*100,aux_us*100,integer=True)
        aux_ls = shunt[106]-auxx
        aux_us = shunt[106]+auxx
        
        if aux_ls < -0.06:
            
            aux_ls = -0.06
            
        if aux_us>0:
            
            aux_us = 0 
            
        shunt[106]= sep.Var(qshunt[106]*100,aux_ls*100,aux_us*100,integer=True)
        
        aux_ls = shunt[109]-auxx
        aux_us = shunt[109]+auxx
        
        if aux_ls < -0.06:
            
            aux_ls = -0.06
            
        if aux_us>0:
            
            aux_us = 0 
            
        shunt[109]= sep.Var(qshunt[109]*100,aux_ls*100,aux_us*100,integer=True)
            

#         shunt[4]= sep.sos1([0,0.4])
#         shunt[33]= sep.sos1([0,-0.06,-0.07,-0.13,-0.14,-0.2])
#         shunt[36]= sep.sos1([0,0.25])
#         shunt[43]= sep.sos1([0,-0.1])
#         shunt[44]= sep.sos1([0,-0.1])
#         shunt[45]= sep.sos1([0,-0.1])
#         shunt[47]= sep.sos1([0,-0.15])

#         shunt[73]= sep.sos1([0,-0.08,-0.12,-0.2])

#         shunt[78]= sep.sos1([0,-0.1,-0.2])

#         shunt[81]= sep.sos1([0,-0.1,-0.2])

#         shunt[82]= sep.sos1([0,-0.1,-0.2])

#         shunt[104]= sep.sos1([0,-0.1,-0.2])

#         shunt[106]= sep.sos1([0,-0.06,-0.07,-0.13,-0.14,-0.2])
#         shunt[109]= sep.sos1([0,-0.06,-0.07,-0.13,-0.14,-0.2])


        

        sep_modelo.trafo['tap_pos'][sep_modelo.trafo['tap_pos']==np.nan] = 1

        tap = sep_modelo.trafo['tap_pos'].to_numpy()

        taps = []
        
        auxx = n_vizinhos_tap
        side = sep_modelo.trafo['tap_side'].values
        i=0
        
        for valor in tap:
            

            if side[i]=='hv' or side[i]=='lv':
                
                
                aux_l = valor-auxx
                aux_u = valor+auxx
                
                if aux_l<0.9:
                    aux_l = 0.9
                    
                if aux_u>1.1:
                    aux_u = 1.1
                    

                taps.append(sep.Var(valor*100,aux_l*100,aux_u*100,integer=True))

#                 taps.append(sep.sos1([0.88,0.8875,0.895,0.9025,0.91,0.9175,0.925,0.9325,0.94,0.9475,0.955,0.9625,0.97,0.9775,0.985,0.9925,1.0,1.0075,1.015,1.0225,1.03,1.0375,1.045,1.0525,1.06,1.0675,1.075,1.0825,1.09,1.0975,1.105,1.1125,1.12]))

            else: 
                taps.append(100)
            
            i=i+1


        if travado == True:
            taps=tap
            shunt = qshunt

        
        qgs = np.zeros(len(pg))

        qg = qgs.tolist()

        pgs = np.zeros(len(pg))
        
        pgs = pgs.tolist()
        

        for bus in sep_modelo.gen['bus'].to_numpy():

            qg[bus] = sep.Var((sep_modelo.res_gen[sep_modelo.gen['bus']==bus]['q_mvar'].to_numpy()/100)[0], (sep_modelo.gen[sep_modelo.gen['bus']==bus]['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.gen[sep_modelo.gen['bus']==bus]['max_q_mvar'].to_numpy()/100)[0] )  
            pgs[bus] = sep.Var(pg[bus],pgs_min[bus],pgs_max[bus])

        qg[68] = sep.Var((sep_modelo.res_ext_grid['q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['max_q_mvar'].to_numpy()/100)[0])
        
        pgs[68] = sep.Var(pg[68],pgs_min[68],pgs_max[68])

        qg[68] = sep.Var((sep_modelo.res_ext_grid['q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['max_q_mvar'].to_numpy()/100)[0])
        
       
        for barra in range(0,len(sep_modelo.bus)):


            sep.Equation(balanco_potencia_reativa(sep,qg, qc, origem.ravel(), destino.ravel(), barra, gkm, bkm, tensoes, angulos, to, td, taps, bkmt,gkmt, bsh, bsht,shunt)==0)
        
        for barra in range(0,len(sep_modelo.bus)):

    
            sep.Equation(balanco_potencia_ativa(sep,pg_sgen,pgs,pc,origem.ravel(), destino.ravel(), barra, gkm, bkm, tensoes, angulos, to,td,taps,bkmt,gkmt,pshunt)==0)
    
        a, equations = perdas(sep,gkm, gkmt, angulos, tensoes, taps, origem.ravel(), destino.ravel(), to, td)
        
        
#         custo

        polycosts = sep_modelo.poly_cost[sep_modelo.poly_cost['et'] != 'ext_grid']

        index_list = sep_modelo.gen.index.values.tolist()

        a_k = []
        b_k = []
        c_k = []

        for val in index_list:

            a_k.append(polycosts['cp2_eur_per_mw2'][polycosts['element']==val].values[0])
            b_k.append(polycosts['cp1_eur_per_mw'][polycosts['element']==val].values[0])
            c_k.append(polycosts['cp0_eur'][polycosts['element']==val].values[0])

        pmax0 = sep_modelo.ext_grid['max_p_mw'].values
        pmin0 = sep_modelo.ext_grid['min_p_mw'].values

        a_k0=sep_modelo.poly_cost['cp2_eur_per_mw2'][sep_modelo.poly_cost['et']=='ext_grid'].values
        b_k0=sep_modelo.poly_cost['cp1_eur_per_mw'][sep_modelo.poly_cost['et']=='ext_grid'].values
        c_k0=sep_modelo.poly_cost['cp0_eur'][sep_modelo.poly_cost['et']=='ext_grid'].values
        e_k0 = (5/100)*((a_k0*pmax0**2 + b_k0*pmax0 + c_k0 + a_k0*pmin0**2 + b_k0*pmin0 + c_k0 )/2)
        f_k0 = (4*np.pi/(sep_modelo.ext_grid['max_p_mw'].values-sep_modelo.ext_grid['min_p_mw'].values))

        e_k = (5/100)*((a_k*sep_modelo.gen['max_p_mw'].values**2 + b_k*sep_modelo.gen['max_p_mw'].values + c_k + a_k*sep_modelo.gen['min_p_mw'].values**2 + b_k*sep_modelo.gen['min_p_mw'].values + c_k )/2)
        f_k = (4*np.pi/(sep_modelo.gen['max_p_mw'].values-sep_modelo.gen['min_p_mw'].values))

        a_k = np.concatenate((a_k0, np.array(a_k)))
        b_k = np.concatenate((b_k0, np.array(b_k)))
        c_k = np.concatenate((c_k0, np.array(c_k)))
        e_k = np.concatenate((e_k0, np.array(e_k)))
        f_k = np.concatenate((f_k0, np.array(f_k)))

             
        alfa = np.zeros(len(sep_modelo.gen['bus'])+1)
        
        alfa = alfa.tolist()
        
        
        alfa[0] = sep.Var(sep.abs(e_k[0]*sep.sin(f_k[0]*(100*minimum[68]-100*pgs[68]))),0,1e15)
        
        u=1
        for bus in sep_modelo.gen['bus'].to_numpy():
            alfa[u]=sep.Var(sep.abs(e_k[u]*sep.sin(f_k[u]*(100*minimum[bus]-100*pgs[bus]))),0,1e15)
            sep.Equation(-alfa[u]<=e_k[u]*sep.sin(f_k[u]*(100*minimum[bus]-100*pgs[bus])))
            sep.Equation(e_k[u]*sep.sin(f_k[u]*(100*minimum[bus]-100*pgs[bus]))<=alfa[u])
            u=u+1
            
        sep.Equation(e_k[0]*sep.sin(f_k[0]*(100*minimum[68]-100*pgs[68]))<=alfa[0])
        
        sep.Equation(-alfa[0]<=e_k[0]*sep.sin(f_k[0]*(100*minimum[68]-100*pgs[68])))
        
        lista_custos = []
        
        u=1
        
        for bus in sep_modelo.gen['bus'].to_numpy():
            
            lista_custos.append(((100*pgs[bus])**2)*a_k[u]+pgs[bus]*b_k[u]*100+c_k[u]+alfa[u])
            u=u+1
        
        lista_custos.append(((100*pgs[68])**2)*a_k[0]+100*pgs[68]*b_k[0]+c_k[0]+alfa[0])
            
        sep.Obj(sep.sum(lista_custos))

            
            
#         sep.Obj(sep.sum(equations)) # perdas
        
#         sep.Obj(voltage_dev(tensoes)) # desvio
        
        
        sep.options.SOLVER = 1
       
        sep.solver_options =['minlp_maximum_iterations 1000',
                             'nlp_maximum_iterations 1000',
                             'minlp_integer_tol 1.0e-5',
                             'objective_convergence_tolerance 1.0e-8',
                             'constraint_convergence_tolerance 1.0e-8',
                             'minlp_branch_method 3' 
                            ]


#         sep.options.RTOL = 1e-6       
#         sep.solver_options = ['tol 1e-8',\
#                                'constr_viol_tol 1e-8',\
#                             'bound_push 1e-4',\
#                             'bound_frac 1e-4']  
        sep.solve(disp=verbose)
        
        tensao = np.zeros(len(sep_modelo.gen['bus'].to_numpy())+1)
        i=1
        tensao[0] = tensoes[68][0]


        for bus in sep_modelo.gen['bus'].to_numpy():

            tensao[i] =  tensoes[bus][0]
            i=i+1

        if travado == False:
            t = np.array([taps[0],taps[1],taps[2],taps[3],taps[4],taps[5],taps[6],taps[8],taps[10]])
        else:
            t = np.array([taps[0],taps[1],taps[2],taps[3],taps[4],taps[5],taps[6],taps[8],taps[10]])

        s = np.zeros(len(sep_modelo.shunt['bus'].to_numpy()))

        i=0
        for bus in sep_modelo.shunt['bus'].to_numpy():


            
            if travado == True:
                s[i] =  shunt[bus]
            else:
                s[i] =  shunt[bus][0]

            i=i+1

        s = s*-1
        
        sep_modelo.res_bus= sep_modelo.res_bus.sort_index()

        thetas = np.zeros(len(angulos))

        voltages = np.zeros(len(angulos))

        pot_reativas = np.zeros(len(qg))


        for i in range(len(angulos)):
            
            if i==0:
                
                thetas[i]=angulos[i]
                
            else: thetas[i]=angulos[i][0]


            voltages[i]=tensoes[i][0]



        sep_modelo.res_bus['vm_pu'] = voltages
        

        sep_modelo.res_ext_grid['p_mw'] = pgs[68][0]*100
        
        u=0
        
        for bus in sep_modelo.gen['bus'].to_numpy():
            
            sep_modelo.gen['p_mw'][int(u)] = pgs[int(bus)][0]*100
   
            u=u+1


        sep_modelo.res_ext_grid['q_mvar'] = qg[68][0]*100

        sep_modelo.shunt['q_mvar'] = s*100
        
        gbest = [tensao,t,s,np.array([0,0,0,0,0,0])]
        
        sep_modelo.trafo['tap_pos'][pd.isnull(sep_modelo.trafo['tap_step_percent'])]=np.nan
        
        sep_modelo.bus['max_vm_pu'] = 1.06

        sep_modelo.bus['min_vm_pu'] = 0.94
        sep_modelo.gen['min_vm_pu'] = 0.94

        sep_modelo.gen['max_vm_pu'] = 1.06
        
        sep_modelo.ext_grid['min_vm_pu'] = 0.94

        sep_modelo.ext_grid['max_vm_pu'] = 1.06
        
    elif len(sep_modelo.bus) == 14:
        
   
        
        tensoes = []
        angulos = []

        sep_modelo.res_bus = sep_modelo.res_bus.sort_index()


        v = sep_modelo.res_bus['vm_pu'].to_numpy()


        theta = np.radians(sep_modelo.res_bus['va_degree'].to_numpy())

        sep = GEKKO(remote=True)
        
        minimum = pgs_min.copy()
        maximum = pgs_max.copy()
                
        # if 0<= pg[1]<=1.4*(1/7):
            
            
        #     pgs_min[1] = 1.4*(0/7)
        #     pgs_max[1] = 1.4*(1/7)
            

            
        # if 1.4*(2/7)<= pg[1]<=1.4*(3/7):

        #     pgs_min[1] = 1.4*(2/7)
        #     pgs_max[1] = 1.4*(3/7)
            
            
            
        # if 1.4*(4/7)<= pg[1]<=1.4*(5/7):

        #     pgs_min[1] = 1.4*(4/7)
        #     pgs_max[1] = 1.4*(5/7)
            
            
        # if 1.4*(6/7)<= pg[1]<=1.4*(7/7):

        #     pgs_min[1] = 1.4*(6/7)
        #     pgs_max[1] = 1.4*(7/7)
            
        
        
                
        # if 0<= pg[2]<=1.*(1/7):
            
            
        #     pgs_min[2] = 1.*(0/7)
        #     pgs_max[2] = 1.*(1/7)
            

            
        # if 1.*(2/7)<= pg[2]<=1.*(3/7):

        #     pgs_min[2] = 1.*(2/7)
        #     pgs_max[2] = 1.*(3/7)
            
            
            
        # if 1.*(4/7)<= pg[2]<=1.*(5/7):

        #     pgs_min[2] = 1.*(4/7)
        #     pgs_max[2] = 1.*(5/7)
            
            
        # if 1.*(6/7)<= pg[2]<=1.*(7/7):

        #     pgs_min[2] = 1.*(6/7)
        #     pgs_max[2] = 1.*(7/7)
            
            
        # if 1.*0<= pg[5]<=1.*(1/7):
            
            
        #     pgs_min[5] = 1.*(0/7)
        #     pgs_max[5] = 1.*(1/7)
            

            
        # if 1.*(2/7)<= pg[5]<=1.*(3/7):

        #     pgs_min[5] = 1.*(2/7)
        #     pgs_max[5] = 1.*(3/7)
            
            
            
        # if 1.*(4/7)<= pg[5]<=1.*(5/7):

        #     pgs_min[5] = 1.*(4/7)
        #     pgs_max[5] = 1.*(5/7)
            
                    
        # if 1.*(6/7)<= pg[5]<=1.*(7/7):

        #     pgs_min[5] = 1.*(6/7)
        #     pgs_max[5] = 1.*(7/7)
            
            
                        
        # if 1.*0<= pg[7] <=1.*(1/7):
            
            
        #     pgs_min[7] = 1.*(0/7)
        #     pgs_max[7] = 1.*(1/7)
            

            
        # if 1.*(2/7)<= pg[7] <=1.*(3/7):

        #     pgs_min[7] = 1.*(2/7)
        #     pgs_max[7] = 1.*(3/7)
            
            
            
        # if 1.*(4/7)<= pg[7] <=1.*(5/7):

        #     pgs_min[7] = 1.*(4/7)
        #     pgs_max[7] = 1.*(5/7)
            
            
        # if 1.*(6/7)<= pg[7] <=1.*(7/7):
            
            
        #     pgs_min[7] = 1.*(6/7)
        #     pgs_max[7] = 1.*(7/7)
            
            
        # if 1.*0<= pg[0] <=3.324*(1/7):
            
            
        #     pgs_min[0] = 3.324*(0/7)
        #     pgs_max[0] = 3.324*(1/7)
            

            
        # if 3.324*(2/7)<= pg[0] <=3.324*(3/7):
            
        #     pgs_min[0] = 3.324*(2/7)
        #     pgs_max[0] = 3.324*(3/7)
            
            
            
        # if 3.324*(4/7)<= pg[0] <=3.324*(5/7):
            
        #     pgs_min[0] = 3.324*(4/7)
        #     pgs_max[0] = 3.324*(5/7)
            
        
        # if 3.324*(6/7)<= pg[0] <=3.324*(7/7):
            
        #     pgs_min[0] = 3.324*(6/7)
        #     pgs_max[0] = 3.324*(7/7)
            


        pg1a = sep.Var(0,0,1, integer=True)
        
        pg1b = sep.Var(0,0,1, integer=True)
        
        pg1c = sep.Var(0,0,1, integer=True)
        
        pg1d = sep.Var(0,0,1, integer=True)
                
        pg2a = sep.Var(0,0,1, integer=True)
        
        pg2b = sep.Var(0,0,1, integer=True)
        
        pg2c = sep.Var(0,0,1, integer=True)        
        
        pg2d = sep.Var(0,0,1, integer=True)
                
        pg5a = sep.Var(0,0,1, integer=True)
        
        pg5b = sep.Var(0,0,1, integer=True)
        
        pg5c = sep.Var(0,0,1, integer=True)        
        
        pg5d = sep.Var(0,0,1, integer=True)

        pg7a = sep.Var(0,0,1, integer=True)
        
        pg7b = sep.Var(0,0,1, integer=True)
        
        pg7c = sep.Var(0,0,1, integer=True)        
        
        pg7d = sep.Var(0,0,1, integer=True)
        
        pg0a = sep.Var(0,0,1, integer=True)
        
        pg0b = sep.Var(0,0,1, integer=True)
        
        pg0c = sep.Var(0,0,1, integer=True)        
        
        pg0d = sep.Var(0,0,1, integer=True)
        
        
        for bus in range(len(sep_modelo.bus)):
            
            tensoes.append(sep.Var(v[bus],0.95,1.05))
            
            if bus>0:
            
                angulos.append(sep.Var(theta[bus],-np.pi,np.pi))
            
            else: angulos.append(0)

        
        shunt = np.zeros(len(sep_modelo.bus)).tolist()
        
        auxx = n_vizinhos_shunt

        aux_ls = qshunt[8]-auxx
        
        aux_us = qshunt[8]+auxx
        
        if aux_ls <= -0.05:
            
            aux_ls = -0.05
            
        if aux_us >= 0:
            
            aux_us = 0
            
        
        shunt[8]= sep.Var(qshunt[8]*100,aux_ls*100,aux_us*100,integer=True)
        
        sep_modelo.trafo['tap_pos'][sep_modelo.trafo['tap_step_percent']==np.nan] = 1.000

        tap = sep_modelo.trafo['tap_pos'].to_numpy()
        
        auxx = n_vizinhos_tap
        
        taps = []
        
        i=0
        
        for valor in tap:

            if i<3:
                
                aux_l = valor-auxx
        
                aux_u = valor+auxx
                
                if aux_l<0.9:
        
                    aux_l = 0.9
                    
                if aux_u>1.1:
        
                    aux_u = 1.1

                taps.append(sep.Var(valor*100,aux_l*100,aux_u*100,integer=True))
            else: 
            
                taps.append(100)
                
                   
            i=i+1

#               taps.append(sep.sos1([0.88,0.8875,0.895,0.9025,0.91,0.9175,0.925,0.9325,0.94,0.9475,0.955,0.9625,0.97,0.9775,0.985,0.9925,1.0,1.0075,1.015,1.0225,1.03,1.0375,1.045,1.0525,1.06,1.0675,1.075,1.0825,1.09,1.0975,1.105,1.1125,1.12]))


        if travado == True:

            taps=tap

            shunt = qshunt

        qgs = np.zeros(len(pg))

        qg = qgs.tolist()

        pgs = np.zeros(len(pg))
        
        pgs = pgs.tolist()

        for bus in sep_modelo.gen['bus'].to_numpy():

            qg[bus] = sep.Var((sep_modelo.res_gen[sep_modelo.gen['bus']==bus]['q_mvar'].to_numpy()/100)[0], (sep_modelo.gen[sep_modelo.gen['bus']==bus]['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.gen[sep_modelo.gen['bus']==bus]['max_q_mvar'].to_numpy()/100)[0] )  

            pgs[bus] = sep.Var(pg[bus],pgs_min[bus],pgs_max[bus])

        qg[0] = sep.Var((sep_modelo.res_ext_grid['q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['max_q_mvar'].to_numpy()/100)[0])
        
        pgs[0] = sep.Var(pg[0],pgs_min[0],pgs_max[0])
                
        for barra in range(0,len(sep_modelo.bus)):

            sep.Equation(balanco_potencia_reativa(sep,qg, qc, origem.ravel(), destino.ravel(), barra, gkm, bkm, tensoes, angulos, to, td, taps, bkmt,gkmt, bsh, bsht,shunt)==0)
        
        for barra in range(0,len(sep_modelo.bus)):
    
            sep.Equation(balanco_potencia_ativa(sep,pg_sgen,pgs,pc,origem.ravel(), destino.ravel(), barra, gkm, bkm, tensoes, angulos, to,td,taps,bkmt,gkmt,pshunt)==0)
    
        a, equations = perdas(sep,gkm, gkmt, angulos, tensoes, taps, origem.ravel(), destino.ravel(), to, td)
        
#         custo

        polycosts = sep_modelo.poly_cost[sep_modelo.poly_cost['et'] != 'ext_grid']

        index_list = sep_modelo.gen.index.values.tolist()

        a_k = []
        b_k = []
        c_k = []

        for val in index_list:

            a_k.append(polycosts['cp2_eur_per_mw2'][polycosts['element']==val].values[0])
            b_k.append(polycosts['cp1_eur_per_mw'][polycosts['element']==val].values[0])
            c_k.append(polycosts['cp0_eur'][polycosts['element']==val].values[0])

        pmax0 = sep_modelo.ext_grid['max_p_mw'].values
        pmin0 = sep_modelo.ext_grid['min_p_mw'].values

        a_k0=sep_modelo.poly_cost['cp2_eur_per_mw2'][sep_modelo.poly_cost['et']=='ext_grid'].values
        b_k0=sep_modelo.poly_cost['cp1_eur_per_mw'][sep_modelo.poly_cost['et']=='ext_grid'].values
        c_k0=sep_modelo.poly_cost['cp0_eur'][sep_modelo.poly_cost['et']=='ext_grid'].values
        e_k0 = (5/100)*((a_k0*pmax0**2 + b_k0*pmax0 + c_k0 + a_k0*pmin0**2 + b_k0*pmin0 + c_k0 )/2)
        f_k0 = (4*np.pi/(sep_modelo.ext_grid['max_p_mw'].values-sep_modelo.ext_grid['min_p_mw'].values))
   
        e_k = (5/100)*((a_k*(sep_modelo.gen['max_p_mw'].values**2) + b_k*sep_modelo.gen['max_p_mw'].values + c_k + a_k*(sep_modelo.gen['min_p_mw'].values**2) + b_k*sep_modelo.gen['min_p_mw'].values + c_k )/2)
        
        f_k = (4*np.pi/(sep_modelo.gen['max_p_mw'].values-sep_modelo.gen['min_p_mw'].values))

        a_k = np.concatenate((a_k0, np.array(a_k)))
        b_k = np.concatenate((b_k0, np.array(b_k)))
        c_k = np.concatenate((c_k0, np.array(c_k)))
        e_k = np.concatenate((e_k0, np.array(e_k)))
        f_k = np.concatenate((f_k0, np.array(f_k)))
        

        alfa = np.zeros(len(sep_modelo.gen['bus'])+1)
        
        alfa = alfa.tolist()
        
        alfa[0] = sep.Var(sep.abs(e_k[0]*sep.sin(f_k[0]*(100*0-100*pgs[0]))),0,1e15)
        
        u=1
        
        for bus in sep_modelo.gen['bus'].to_numpy():
        
            alfa[u]=sep.Var(sep.abs(e_k[u]*sep.sin(f_k[u]*(100*0-100*pgs[bus]))),0,1e15)
            sep.Equation(-alfa[u]<=e_k[u]*sep.sin(f_k[u]*(100*0-100*pgs[bus])))
            sep.Equation(e_k[u]*sep.sin(f_k[u]*(100*0-100*pgs[bus]))<=alfa[u])
            u=u+1
            
        sep.Equation(e_k[0]*sep.sin(f_k[0]*(100*0-100*pgs[0]))<=alfa[0])
        sep.Equation(-alfa[0]<=e_k[0]*sep.sin(f_k[0]*(100*0-100*pgs[0])))
        lista_custos = []
        
        u=1
        
        for bus in sep_modelo.gen['bus'].to_numpy():
            
            lista_custos.append(((100*pgs[bus])**2)*a_k[u]+pgs[bus]*b_k[u]*100+c_k[u]+alfa[u])
            u=u+1
        
        lista_custos.append(((100*pgs[0])**2)*a_k[0]+100*pgs[0]*b_k[0]+c_k[0]+alfa[0])
            
        sep.Obj(sep.sum(lista_custos))
        

        delta1 = 1.4

        sep.Equation(pgs[1] >= 0*delta1*pg1a + (2/7)*delta1*pg1b + (4/7)*delta1*pg1c + (6/7)*delta1*pg1d)

        sep.Equation(pgs[1] <= (1/7)*delta1*pg1a + (3/7)*delta1*pg1b + (5/7)*delta1*pg1c + (7/7)*delta1*pg1d)
        
        sep.Equation(pg1a+pg1b+pg1c+pg1d==1)

        delta2 = 1.0

        sep.Equation(pgs[2] >= 0*delta2*pg2a + (2/7)*delta2*pg2b + (4/7)*delta2*pg2c + (6/7)*delta2*pg2d)

        sep.Equation(pgs[2] <= (1/7)*delta2*pg2a + (3/7)*delta2*pg2b + (5/7)*delta2*pg2c + (7/7)*delta2*pg2d)
                
        sep.Equation(pg2a+pg2b+pg2c+pg2d==1)

        delta5 = 1.0
        
        sep.Equation(pgs[5] >= 0*delta5*pg5a + (2/7)*delta5*pg5b + (4/7)*delta5*pg5c + (6/7)*delta5*pg5d)

        sep.Equation(pgs[5] <= (1/7)*delta5*pg5a + (3/7)*delta5*pg5b + (5/7)*delta5*pg5c + (7/7)*delta5*pg5d)
                
        sep.Equation(pg5a+pg5b+pg5c+pg5d==1)
        
        delta7 = 1.0
        
        sep.Equation(pgs[7] >= 0*delta7*pg7a + (2/7)*delta7*pg7b + (4/7)*delta7*pg7c + (6/7)*delta7*pg7d)

        sep.Equation(pgs[7] <= (1/7)*delta7*pg7a + (3/7)*delta7*pg7b + (5/7)*delta7*pg7c + (7/7)*delta7*pg7d)
                
        sep.Equation(pg7a+pg7b+pg7c+pg7d==1)
        
        delta0 = 3.324
        
        sep.Equation(pgs[0] >= 0*delta0*pg0a + (2/7)*delta0*pg0b + (4/7)*delta0*pg0c + (6/7)*delta0*pg0d)

        sep.Equation(pgs[0] <= (1/7)*delta0*pg0a + (3/7)*delta0*pg0b + (5/7)*delta0*pg0c + (7/7)*delta0*pg0d)
                
        sep.Equation(pg0a+pg0b+pg0c+pg0d==1)
        
        
        
#         sep.Obj(sep.sum(equations)) # perdas
        
        
#         sep.Obj(voltage_dev(tensoes)) # desvio de Tensão
        
        sep.options.SOLVER = 1
#         sep.options.RTOL = 1e-6
        
#sep.solver_options = ['tol 1e-8',\
#                                'constr_viol_tol 1e-8',\
#                             'bound_push 1e-6',\
#                             'bound_frac 1e-6']  
        
        sep.solver_options =[
                            
                            'minlp_maximum_iterations 1000',
                             'nlp_maximum_iterations 1000',
                             'minlp_integer_tol 1.0e-5',
                             'objective_convergence_tolerance 1.0e-8',
                             'constraint_convergence_tolerance 1.0e-8',
                             'minlp_branch_method 3' 
                            
                            ]
    
        sep.solve(disp=verbose)
        
        tensao = np.zeros(len(sep_modelo.gen['bus'].to_numpy())+1)
        i=1
        tensao[0] = tensoes[0][0]

        for bus in sep_modelo.gen['bus'].to_numpy():

            tensao[i] =  tensoes[bus][0]
            i=i+1
        if travado == False:
            t =np.array([taps[0],taps[1],taps[2]])
        else:
            t =np.array([taps[0],taps[1],taps[2]])

        s = np.zeros(len(sep_modelo.shunt['bus'].to_numpy()))

        i=0
        for bus in sep_modelo.shunt['bus'].to_numpy():


            
            if travado == True:
                s[i] =  shunt[bus]
            else:
                s[i] =  shunt[bus][0]

            i=i+1


        s = s*-1
        
        sep_modelo.res_bus= sep_modelo.res_bus.sort_index()

        thetas = np.zeros(len(angulos))

        voltages = np.zeros(len(angulos))

        pot_reativas = np.zeros(len(qg))

        for i in range(len(angulos)):
            if i>0:
                
                thetas[i]=angulos[i][0]
            else: thetas[i]=angulos[i]
                
            voltages[i]=tensoes[i][0]



        sep_modelo.res_bus['vm_pu'] = voltages


        sep_modelo.res_ext_grid['p_mw'] = pgs[0][0]*100
        
        u=0
        for bus in sep_modelo.gen['bus'].to_numpy():
            
            sep_modelo.gen['p_mw'][u] = pgs[bus][0]*100
            u=u+1


        sep_modelo.res_ext_grid['q_mvar'] = qg[0][0]*100

        sep_modelo.shunt['q_mvar'] = s
        
       

        gbest = [tensao,t,s,np.array([0,0,0,0,0,0])]
        
        sep_modelo.trafo['tap_pos'][pd.isnull(sep_modelo.trafo['tap_step_percent'])]=np.nan
        
        sep_modelo.bus['max_vm_pu'] = 1.05

        sep_modelo.bus['min_vm_pu'] = 0.95
        sep_modelo.gen['min_vm_pu'] = 0.95

        sep_modelo.gen['max_vm_pu'] = 1.05
        
        sep_modelo.ext_grid['min_vm_pu'] = 0.95

        sep_modelo.ext_grid['max_vm_pu'] = 1.05

    elif len(sep_modelo.bus) == 30:

        tensoes = []
        
        angulos = []

        sep_modelo.res_bus = sep_modelo.res_bus.sort_index()


        v = sep_modelo.res_bus['vm_pu'].to_numpy()


        theta = np.radians(sep_modelo.res_bus['va_degree'].to_numpy())
        
        sep = GEKKO(remote=True)

        for bus in range(len(sep_modelo.bus)):

            tensoes.append(sep.Var(v[bus],0.95,1.05))
            if bus>0:
                
                angulos.append(sep.Var(theta[bus],-np.pi,np.pi))
            else: angulos.append(0)
                
        shunt = np.zeros(len(sep_modelo.bus)).tolist()
        
        auxx = n_vizinhos_shunt
        aux_ls = qshunt[9]-auxx
        aux_us = qshunt[9]+auxx
        
        if aux_ls <= -0.05:
            
            aux_ls = -0.05
            
        if aux_us >= 0:
            
            aux_us = 0
            
        
        shunt[9]= sep.Var(qshunt[9]*100,aux_ls*100,aux_us*100,integer=True)
        
        
        aux_ls = qshunt[23]-auxx
        aux_us = qshunt[23]+auxx
        
        if aux_ls <= -0.05:
            
            aux_ls = -0.05
            
        if aux_us >= 0:
            
            aux_us = 0
            
        
        shunt[23]= sep.Var(qshunt[23]*100,aux_ls*100,aux_us*100,integer=True)
        
        sep_modelo.trafo['tap_pos'][sep_modelo.trafo['tap_pos']==np.nan] = 1
        auxx=n_vizinhos_tap
        tap = sep_modelo.trafo['tap_pos'].to_numpy()
        side = sep_modelo.trafo['tap_side'].values
        taps = []
        i=0
        for valor in tap:
            
            if side[i]=='hv' or side[i]=='lv':
                                
                aux_l = valor-auxx
                aux_u = valor+auxx
                
                if aux_l<0.9:
                    aux_l = 0.9
                    
                if aux_u>1.1:
                    aux_u = 1.1
                    

                taps.append(sep.Var(valor*100,aux_l*100,aux_u*100,integer=True))

#               taps.append(sep.sos1([0.88,0.8875,0.895,0.9025,0.91,0.9175,0.925,0.9325,0.94,0.9475,0.955,0.9625,0.97,0.9775,0.985,0.9925,1.0,1.0075,1.015,1.0225,1.03,1.0375,1.045,1.0525,1.06,1.0675,1.075,1.0825,1.09,1.0975,1.105,1.1125,1.12]))

            else: 
                taps.append(100)
            
            i=i+1

        if travado == True:
            
            taps=tap
            
            shunt = qshunt


        qgs = np.zeros(len(pg))

        qg = qgs.tolist()

        pgs = np.zeros(len(pg))
        
        pgs = pgs.tolist()
        
        maximum = pgs_max.copy()
        
        minimum = pgs_min.copy()
        

        # if 0<= pg[1]<=0.21:
            
            
        #     pgs_min[1] = 0
            
        #     pgs_max[1] = 0.21
            
            
        # if 0.24<= pg[1]<=0.45:

        #     pgs_min[1] = 0.24
            
        #     pgs_max[1] = 0.45

            
        # if 0.55<= pg[1]<=maximum[1]:

        #     pgs_min[1] = 0.55
            
        #     pgs_max[1] = maximum[1]
            
            
   

        # # ### gerador 2

        # if 0<= pg[4]<=0.3:

            
        #     pgs_min[4] = 0.0
        #     pgs_max[4] = 0.3
            
            
        # if 0.36<= pg[4]<=maximum[4]:

            
        #     pgs_min[4] = 0.36
        #     pgs_max[4] = maximum[4]
            
        # # ### gerador 3
        # # ## 0-0.25, 0.3-max
        
        # if 0<= pg[7]<=0.25:

            
        #     pgs_min[7] = 0.0
        #     pgs_max[7] = 0.25
            
            
        # if 0.3<= pg[7]<=maximum[7]:

            
        #     pgs_min[7] = 0.3
        #     pgs_max[7] = maximum[7]
            
         

        # ### gerador 4
        
        # # 0-0.25, 0.28-max

        # if 0 <= pg[10] <= 0.25:

            
        #     pgs_min[10] = 0.0
        #     pgs_max[10] = 0.25
            
            
        # if 0.28<= pg[10]<=maximum[10]:

            
        #     pgs_min[10] = 0.28
        #     pgs_max[10] = maximum[10]
            
        
                
        # # ### gerador 5
        # # # 0. - 0.24, 0.3-max

        # if 0<= pg[12]<=0.24:

            
        #     pgs_min[12] = 0.0
        #     pgs_max[12] = 0.24
            
        # if 0.3<= pg[12]<=maximum[12]:

            
        #     pgs_min[12] = 0.3
        #     pgs_max[12] = maximum[12]
            
            
            
        # # #0-0.55, 0.66-0.8 1.2-max 

        # if 0<= pg[0]<=0.55:
            
        #     pgs_min[0] = 0
        #     pgs_max[0] = 0.55
            
        # if 0.66<= pg[0]<=0.8:
            
        #     pgs_min[0] = 0.66
        #     pgs_max[0] = 0.8
            
            
        # if 1.2<= pg[0]<=maximum[0]:
            
        #     pgs_min[0] = 1.2
        #     pgs_max[0] = maximum[0]

    
        for bus in sep_modelo.gen['bus'].to_numpy():

            qg[bus] = sep.Var((sep_modelo.res_gen[sep_modelo.gen['bus']==bus]['q_mvar'].to_numpy()/100)[0], (sep_modelo.gen[sep_modelo.gen['bus']==bus]['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.gen[sep_modelo.gen['bus']==bus]['max_q_mvar'].to_numpy()/100)[0] )  
            pgs[bus] = sep.Var(pg[bus],pgs_min[bus],pgs_max[bus])

        qg[0] = sep.Var((sep_modelo.res_ext_grid['q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['max_q_mvar'].to_numpy()/100)[0])
        
        pgs[0] = sep.Var(pg[0],pgs_min[0],pgs_max[0])
        
        
         
        
        for barra in range(0,len(sep_modelo.bus)):

            sep.Equation(balanco_potencia_reativa(sep,qg, qc, origem.ravel(), destino.ravel(), barra, gkm, bkm, tensoes, angulos, to, td, taps, bkmt,gkmt, bsh, bsht,shunt)==0)
        
        for barra in range(0,len(sep_modelo.bus)):
    
            sep.Equation(balanco_potencia_ativa(sep,pg_sgen,pgs,pc,origem.ravel(), destino.ravel(), barra, gkm, bkm, tensoes, angulos, to,td,taps,bkmt,gkmt,pshunt)==0)
    
        a, equations = perdas(sep,gkm, gkmt, angulos, tensoes, taps, origem.ravel(), destino.ravel(), to, td)
        
#         custo

        polycosts = sep_modelo.poly_cost[sep_modelo.poly_cost['et'] != 'ext_grid']

        index_list = sep_modelo.gen.index.values.tolist()

        a_k = []
        b_k = []
        c_k = []

        for val in index_list:

            a_k.append(polycosts['cp2_eur_per_mw2'][polycosts['element']==val].values[0])
            b_k.append(polycosts['cp1_eur_per_mw'][polycosts['element']==val].values[0])
            c_k.append(polycosts['cp0_eur'][polycosts['element']==val].values[0])

        pmax0 = sep_modelo.ext_grid['max_p_mw'].values
        pmin0 = sep_modelo.ext_grid['min_p_mw'].values

        a_k0=sep_modelo.poly_cost['cp2_eur_per_mw2'][sep_modelo.poly_cost['et']=='ext_grid'].values
        b_k0=sep_modelo.poly_cost['cp1_eur_per_mw'][sep_modelo.poly_cost['et']=='ext_grid'].values
        c_k0=sep_modelo.poly_cost['cp0_eur'][sep_modelo.poly_cost['et']=='ext_grid'].values
        e_k0 = (5/100)*((a_k0*sep_modelo.ext_grid['max_p_mw'].values**2 + b_k0*sep_modelo.ext_grid['max_p_mw'].values + c_k0 + a_k0*sep_modelo.ext_grid['min_p_mw'].values**2 + b_k0*sep_modelo.ext_grid['min_p_mw'].values + c_k0 )/2)
        f_k0 = (4*np.pi/(sep_modelo.ext_grid['max_p_mw'].values-sep_modelo.ext_grid['min_p_mw'].values))

        e_k = (5/100)*((a_k*sep_modelo.gen['max_p_mw'].values**2 + b_k*sep_modelo.gen['max_p_mw'].values + c_k + a_k*sep_modelo.gen['min_p_mw'].values**2 + b_k*sep_modelo.gen['min_p_mw'].values + c_k )/2)
        
        f_k = (4*np.pi/(sep_modelo.gen['max_p_mw'].values-sep_modelo.gen['min_p_mw'].values))

        a_k = np.concatenate((a_k0, np.array(a_k)))
        b_k = np.concatenate((b_k0, np.array(b_k)))
        c_k = np.concatenate((c_k0, np.array(c_k)))
        e_k = np.concatenate((e_k0, np.array(e_k)))
        f_k = np.concatenate((f_k0, np.array(f_k)))
        

        alfa = np.zeros(len(sep_modelo.gen['bus'])+1)
        
        alfa = alfa.tolist()
        
        
        alfa[0] = sep.Var(np.abs(e_k[0]*sep.sin(f_k[0]*(100*minimum[0]-100*pgs[0]))),0,1e15)
        
        u=1
        for bus in sep_modelo.gen['bus'].to_numpy():
            alfa[u]=sep.Var(sep.abs(e_k[u]*sep.sin(f_k[u]*(100*minimum[bus]-100*pgs[bus]))),0,1e15)
            sep.Equation(-alfa[u]<=e_k[u]*sep.sin(f_k[u]*(100*minimum[bus]-100*pgs[bus])))
            sep.Equation(e_k[u]*sep.sin(f_k[u]*(100*minimum[bus]-100*pgs[bus]))<=alfa[u])
            u=u+1
            
        sep.Equation(e_k[0]*sep.sin(f_k[0]*(100*pmin0[0]-100*pgs[0]))<=alfa[0])
        sep.Equation(-alfa[0]<=e_k[0]*sep.sin(f_k[0]*(100*pmin0[0]-100*pgs[0])))
        
        lista_custos = []
        
        u=1
        for bus in sep_modelo.gen['bus'].to_numpy():
            
            lista_custos.append(((100*pgs[bus])**2)*a_k[u]+pgs[bus]*b_k[u]*100+c_k[u]+alfa[u])
            u=u+1
            
            
        sep.solver_options =[
                            
                            'minlp_maximum_iterations 1000',
                             'nlp_maximum_iterations 1000',
                             'minlp_integer_tol 1.0e-5',
                             'objective_convergence_tolerance 1.0e-7',
                             'constraint_convergence_tolerance 1.0e-7',
                             'minlp_branch_method 3' 
                            
                            ]
        
        
        lista_custos.append(((100*pgs[0])**2)*a_k[0]+100*pgs[0]*b_k[0]+c_k[0]+alfa[0])
        
        aa = sep.Var(0,0,1, integer=True)
        
        bb = sep.Var(0,0,1, integer=True)
        
        cc = sep.Var(0,0,1, integer=True)
        
        dd = sep.Var(0,0,1, integer=True)
        
        ee = sep.Var(0,0,1, integer=True)
        
        ff = sep.Var(0,0,1, integer=True)
        
        gg = sep.Var(0,0,1, integer=True)
        
        hh = sep.Var(0,0,1, integer=True)
        
        ii = sep.Var(0,0,1, integer=True)
        
        jj = sep.Var(0,0,1, integer=True)
        
        kk = sep.Var(0,0,1, integer=True)
        
        ll = sep.Var(0,0,1, integer=True)
        
        mm = sep.Var(0,0,1, integer=True)
        
        nn = sep.Var(0,0,1, integer=True)
        
        sep.Equation(pgs[1] >= 0*aa + 0.24*bb + 0.55*cc)
        
        sep.Equation(pgs[1] <= 0.21*aa + 0.45*bb + pgs_max[1]*cc)
        
        sep.Equation(aa+bb+cc==1)
       
        sep.Equation(pgs[0] >= 0*dd + 0.66*ee + 1.2*ff)
        
        sep.Equation(pgs[0] <= 0.55*dd + 0.8*ee + pgs_max[0]*ff)
        
        sep.Equation(dd+ee+ff==1)
            
        sep.Equation(pgs[4] >= 0*gg + 0.36*hh)
        
        sep.Equation(pgs[4] <= 0.30*gg + pgs_max[4]*hh)
        
        sep.Equation(gg+hh==1)
        
        sep.Equation(pgs[7] >= 0*ii + 0.30*jj)
        
        sep.Equation(pgs[7] <= 0.25*ii + pgs_max[7]*jj)
        
        sep.Equation(ii+jj==1)
        
        sep.Equation(pgs[10] >= 0*kk + 0.28*ll)
        
        sep.Equation(pgs[10] <= 0.25*kk + pgs_max[10]*ll)
        
        sep.Equation(kk+ll==1)
        
        sep.Equation(pgs[12] >= 0*mm + 0.30*nn)
        
        sep.Equation(pgs[12] <= 0.24*mm + pgs_max[12]*nn)
        
        sep.Equation(mm+nn==1)

        sep.Obj(sep.sum(lista_custos))        
                
        sep.options.SOLVER = 1

        sep.solve(disp=verbose)

        tensao = np.zeros(len(sep_modelo.gen['bus'].to_numpy())+1)
        i=1
        tensao[0] = tensoes[0][0]

        for bus in sep_modelo.gen['bus'].to_numpy():

            tensao[i] =  tensoes[bus][0]
            i=i+1
        
        if travado == False:
            t =np.array([taps[0],taps[1],taps[4],taps[6]])
        else:
            t =np.array([taps[0],taps[1],taps[4],taps[6]])

        s = np.zeros(len(sep_modelo.shunt['bus'].to_numpy()))

        i=0
        
        for bus in sep_modelo.shunt['bus'].to_numpy():

            if travado == True:
                s[i] =  shunt[bus]
            else:
                s[i] =  shunt[bus][0]

            i=i+1

        s = s*-1
        
        sep_modelo.res_bus= sep_modelo.res_bus.sort_index()

        thetas = np.zeros(len(angulos))

        voltages = np.zeros(len(angulos))


        for i in range(len(angulos)):
            if i==0:
                thetas[i]=angulos[i]
                
            else:thetas[i]=angulos[i][0]

            voltages[i]=tensoes[i][0]



        sep_modelo.res_bus['vm_pu'] = voltages


        sep_modelo.res_ext_grid['p_mw'] = pgs[0][0]*100
        
        u=0
        for bus in sep_modelo.gen['bus'].to_numpy():
            
            sep_modelo.gen['p_mw'][u] = pgs[bus][0]*100
            u=u+1


        sep_modelo.res_ext_grid['q_mvar'] = qg[0][0]*100

        sep_modelo.shunt['q_mvar'] = s/100
        
        sep_modelo.trafo['tap_pos'][pd.isnull(sep_modelo.trafo['tap_step_percent'])]=np.nan
        
        sep_modelo.bus['max_vm_pu'] = 1.05

        sep_modelo.bus['min_vm_pu'] = 0.95

        sep_modelo.gen['min_vm_pu'] = 0.95

        sep_modelo.gen['max_vm_pu'] = 1.05
        
        sep_modelo.ext_grid['min_vm_pu'] = 0.95

        sep_modelo.ext_grid['max_vm_pu'] = 1.05
        

    volt = []
    angl = []
    ppower = []
    qpower = []
    tr = []
    bbsh = []

    for i in range(len(tensoes)):

        try: volt.append(tensoes[i][0])
        
        except: volt.append(tensoes[i])

    
    for i in range(len(angulos)):

        try: angl.append(angulos[i][0])
        
        except: angl.append(angulos[i])

    
    for i in range(len(pgs)):

        try: ppower.append(pgs[i][0])
        
        except: ppower.append(pgs[i])    
    
    
    for i in range(len(qg)):

        try: qpower.append(qg[i][0])
        
        except: qpower.append(qg[i])  

    for i in range(len(s)):

        try: bbsh.append(s[i][0]/Sbase)
        
        except: bbsh.append(s[i]/Sbase)    

    for i in range(len(t)):

        try: tr.append(t[i][0]/100)
        
        except: tr.append(t[i]/100)    
    
    results = {'system':sep_modelo,
                'voltages':volt,
                'angles':angl,
                'trafos':tr,
                'shunts':bbsh,
                'active_power':ppower,
                'reactive_power':qpower,
                'objective_function':sep.options.OBJFCNVAL,
                'iterations':sep.options.ITERATIONS,
                'a_k':a_k,
                'b_k':b_k,
                'c_k':c_k,
                'e_k':e_k,
                'f_k':f_k
                }

    return results