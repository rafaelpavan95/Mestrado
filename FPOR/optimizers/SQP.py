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
            
#             print(linhas[(barras_destino==bd) & (barras_origem==barra_atual)],'aaaaaaaaaa')
            posi = linhas[(barras_destino==bd) & (barras_origem==barra_atual)][1]
            
#         print(bd,'vai linha')
#         print(posi)
        soma = soma + gkm_linhas[posi]*(tensoes[barra_atual]**2) - tensoes[barra_atual]*tensoes[bd]*(gkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(gkm_linhas[posi]*(tensoes[barra_atual]**2) - tensoes[barra_atual]*tensoes[bd]*(gkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    for bd in barras_origem[barras_destino==barra_atual]:
        
        posi = linhas[(barras_destino==barra_atual) & (barras_origem==bd)][0]
        
        baux.append(bd)
        
        if baux[-2]==bd:
#             print(linhas[(barras_destino==barra_atual) & (barras_origem==bd)][1],'aaaaaaaaaa')
            posi = linhas[(barras_destino==barra_atual) & (barras_origem==bd)][1]
        
#         print(bd,'volta linha')
#         print(posi)
        soma = soma + gkm_linhas[posi]*(tensoes[barra_atual]**2) - tensoes[barra_atual]*tensoes[bd]*(gkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(gkm_linhas[posi]*(tensoes[barra_atual]**2) - tensoes[barra_atual]*tensoes[bd]*(gkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
    
    linhas = np.arange(0,len(to),1)
    
    for bd in td[to==barra_atual]:
        
        posi = linhas[(td==bd) & (to==barra_atual)][0]
#         print(bd,'vai trafo')
#         print(posi)
        soma = soma + (gkmt[posi]*tensoes[barra_atual]*tensoes[barra_atual]/(tap[posi]/100)**2 - (tensoes[barra_atual]*tensoes[bd]/(tap[posi]/100))*(gkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
        fluxos.append(gkmt[posi]*tensoes[barra_atual]*tensoes[barra_atual]/(tap[posi]/100)**2 - (tensoes[barra_atual]*tensoes[bd]/(tap[posi]/100))*(gkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-bkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
    
    for bd in to[td==barra_atual]:
        
        posi = linhas[(td==barra_atual) & (to==bd)][0]
#         print(bd,'volta trafo')
#         print(posi)
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
            
#             print(linhas[(barras_destino==bd) & (barras_origem==barra_atual)],'aaaaaaaaaa')
            posi = linhas[(barras_destino==bd) & (barras_origem==barra_atual)][1]
            
#         print(bd,'vai linha')
#         print(posi)
        
        soma = soma + -(-bkm_linhas[posi]+bshl[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*tensoes[bd]*(-bkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(-(-bkm_linhas[posi]+bshl[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*tensoes[bd]*(-bkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    for bd in barras_origem[barras_destino==barra_atual]:
        
        posi = linhas[(barras_destino==barra_atual) & (barras_origem==bd)][0]
        
        baux.append(bd)
        
        if baux[-2]==bd:
#             print(linhas[(barras_destino==barra_atual) & (barras_origem==bd)][1],'aaaaaaaaaa')
            posi = linhas[(barras_destino==barra_atual) & (barras_origem==bd)][1]
        
#         print(bd,'volta linha')
#         print(posi)
        
        soma = soma + -(-bkm_linhas[posi]+bshl[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*tensoes[bd]*(-bkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(-(-bkm_linhas[posi]+bshl[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*tensoes[bd]*(-bkm_linhas[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkm_linhas[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    linhas = np.arange(0,len(to),1)
    
    for bd in td[to==barra_atual]:
        
        posi = linhas[(td==bd) & (to==barra_atual)][0]
#         print(bd,'vai trafo')
#         print(posi)
        soma = soma + -(-bkmt[posi]/((tap[posi]/100)**2)+bsht[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*(1/(tap[posi]/100))*tensoes[bd]*(-bkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(-(-bkmt[posi]/((tap[posi]/100)**2)+bsht[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*(1/(tap[posi]/100))*tensoes[bd]*(-bkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    for bd in to[td==barra_atual]:
        
        posi = linhas[(td==barra_atual) & (to==bd)][0]
#         print(bd,'volta trafo')
#         print(posi)
        soma = soma + -(-bkmt[posi]+bsht[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*(1/(tap[posi]/100))*tensoes[bd]*(-bkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd]))
        fluxos.append(-(-bkmt[posi]+bsht[posi]/2)*tensoes[barra_atual]**2+tensoes[barra_atual]*(1/(tap[posi]/100))*tensoes[bd]*(-bkmt[posi]*sep.cos(angulos[barra_atual]-angulos[bd])-gkmt[posi]*sep.sin(angulos[barra_atual]-angulos[bd])))
        
    return  qg[barra_atual] - qc[barra_atual] - soma - (qshunt[barra_atual]/100)*tensoes[barra_atual]**2




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





def Sequential_Quadratic_Programming(sep_teste, verbose=True, travado=False):
    
    
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
    
    pg_sgen = np.zeros(len(sep_modelo.bus))
    qg_sgen = np.zeros(len(sep_modelo.bus))
    
    i = 0

    sep_modelo.gen = sep_modelo.gen.sort_index()

    sep_modelo.res_gen = sep_modelo.res_gen.sort_index()

    for bus in sep_modelo.gen['bus'].to_numpy():

        pg[bus] = sep_modelo.gen['p_mw'].to_numpy()[i]/Sbase
        
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
        
    if len(sep_modelo.bus)==300:

        pg[256] = sep_modelo.res_ext_grid['p_mw'].to_numpy()/100
        qg[256] = sep_modelo.res_ext_grid['q_mvar'].to_numpy()/100
        
    if len(sep_modelo.bus)==14 or len(sep_modelo.bus)==30 :

        pg[0] = sep_modelo.res_ext_grid['p_mw'].to_numpy()/100

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

            tensoes.append(sep.Var(v[bus],0.94,1.06))

            if bus > 0:

                angulos.append(sep.Var(theta[bus],-np.pi,np.pi))

            else: angulos.append(0)

        shunt = np.zeros(len(sep_modelo.bus)).tolist()
        
        shunt[4]= sep.Var(qshunt[4]*100,0*100,0.40*100, integer=False)
        shunt[33]= sep.Var(qshunt[33]*100,-0.2*100,0*100, integer=False)
        shunt[36]= sep.Var(qshunt[36]*100,0*100,0.25*100, integer=False)
        shunt[43]= sep.Var(qshunt[43]*100,-0.1*100,0*100, integer=False)
        shunt[44]= sep.Var(qshunt[44]*100,-0.1*100,0*100, integer=False)
        shunt[45]= sep.Var(qshunt[45]*100,-0.1*100,0*100, integer=False)
        shunt[47]= sep.Var(qshunt[47]*100,-0.15*100,0*100, integer=False)
        shunt[73]= sep.Var(qshunt[73]*100,-0.2*100,0*100, integer=False)
        shunt[78]= sep.Var(qshunt[78]*100,-0.2*100,0*100, integer=False)
        shunt[81]= sep.Var(qshunt[81]*100,-0.2*100,0*100, integer=False)
        shunt[82]= sep.Var(qshunt[82]*100,-0.2*100,0*100, integer=False)
        shunt[104]= sep.Var(qshunt[104]*100,-0.2*100,0*100, integer=False)
        shunt[106]= sep.Var(qshunt[106]*100,-0.2*100,0*100, integer=False)
        shunt[109]= sep.Var(qshunt[109]*100,-0.2*100,0*100, integer=False)

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

        side = sep_modelo.trafo['tap_side'].values

        i=0
        
        for valor in tap:
    
            if side[i]=='hv' or side[i]=='lv':
                
                
                taps.append(sep.Var(valor*100,90,110,integer=False))


            else: 

                taps.append(100)
            
            i=i+1

    

        qgs = np.zeros(len(pg))

        qg = qgs.tolist()


        for bus in sep_modelo.gen['bus'].to_numpy():

            qg[bus] = sep.Var((sep_modelo.res_gen[sep_modelo.gen['bus']==bus]['q_mvar'].to_numpy()/100)[0], (sep_modelo.gen[sep_modelo.gen['bus']==bus]['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.gen[sep_modelo.gen['bus']==bus]['max_q_mvar'].to_numpy()/100)[0] )  


        qg[68] = sep.Var((sep_modelo.res_ext_grid['q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['min_q_mvar'].to_numpy()/100)[0],(sep_modelo.ext_grid['max_q_mvar'].to_numpy()/100)[0])
        
        pgl=np.copy(pg)
        pgs = pgl.tolist()

        pgs[68] = sep.Var((sep_modelo.res_ext_grid['p_mw'].to_numpy()/100)[0],(sep_modelo.ext_grid['min_p_mw'].to_numpy()/100)[0],(sep_modelo.ext_grid['max_p_mw'].to_numpy()/100)[0])
        
        for barra in range(0,len(sep_modelo.bus)):


            sep.Equation(balanco_potencia_reativa(sep,qg, qc, origem.ravel(), destino.ravel(), barra, gkm, bkm, tensoes, angulos, to, td, taps, bkmt,gkmt, bsh, bsht,shunt)==0)
        
        for barra in range(0,len(sep_modelo.bus)):

    
            sep.Equation(balanco_potencia_ativa(sep,pg_sgen,pgs,pc,origem.ravel(), destino.ravel(), barra, gkm, bkm, tensoes, angulos, to,td,taps,bkmt,gkmt,pshunt)==0)
    
        a, equations = perdas(sep,gkm, gkmt, angulos, tensoes, taps, origem.ravel(), destino.ravel(), to, td)
        
        sep.Obj(sep.sum(equations))
        
        sep.options.SOLVER = 1

        sep.solve(disp=verbose)
        
               
        sep.solver_options =['minlp_maximum_iterations 100000',
                             'nlp_maximum_iterations 100000',
                             'objective_convergence_tolerance 5.0e-6',
                             'constraint_convergence_tolerance 5.0e-6',
                             
                            ]

    trafos = []

    for val in taps:

        try: trafos.append(val[0]/100)


        except: trafos.append(val/100)


    voltages = []

    for val in tensoes:

        try: voltages.append(val[0])


        except: voltages.append(val)


    angul = []

    for val in angulos:

        try: angul.append(val[0])


        except: angul.append(val)


    rpower = []

    for val in qg:

        try: rpower.append(val[0])


        except: rpower.append(val)
    
    bsh = []

    
    for val in shunt:

        try: bsh.append(val[0])


        except: bsh.append(val)
    

    results = {'system':sep_modelo,
                'voltages':voltages,
                'angles':angul,
                'trafos':trafos,
                'shunts':bsh,
                'reactive_power':rpower,
                'objective_function':sep.options.OBJFCNVAL,
                'iterations':sep.options.ITERATIONS,
                } 

                
    return results


