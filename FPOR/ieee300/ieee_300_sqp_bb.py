
import numpy as np
import pandas as pd
from gekko import GEKKO
import time 



def calcula_perdas(sep,barra_origem_linhas,barra_destino_linhas,barra_origem_trafos,barra_destino_trafos,gkm_linhas,gkm_trafos,v_total,theta_total,taps):
    i = 0

    soma = 0
    eq = []
    for linha in zip(barra_origem_linhas, barra_destino_linhas):

        eq.append(gkm_linhas[i]*(v_total[linha[0]]**2 + v_total[linha[1]]**2 - 2*v_total[linha[0]]*v_total[linha[1]]*(sep.cos(theta_total[linha[0]]-theta_total[linha[1]]))))
        i=i+1
    i=0
    for linha_t in zip(barra_origem_trafos, barra_destino_trafos):

        eq.append(gkm_trafos[i]*(((v_total[linha_t[0]])/(taps[i]/100))**2 + v_total[linha_t[1]]**2 - 2*(1/(taps[i]/100))*v_total[linha_t[0]]*v_total[linha_t[1]]*(sep.cos(theta_total[linha_t[0]]-theta_total[linha_t[1]]))))
        i=i+1

    return eq



def calcula_fluxo_p(sep,barra,v_total,theta_total,barra_destino_trafos,barra_origem_trafos,barra_destino_linhas,barra_origem_linhas,gkm_linhas,gkm_trafos,bkm_linhas,bkm_trafos,taps,p_gerada,p_demanda,p_gs):

    soma = 0

    barra = int(barra)

    baux = []
    baux.append(33333333)
    baux.append(33333331)

    for bd in barra_destino_linhas[barra_origem_linhas==barra]:

        gkm = gkm_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][0]
        bkm = bkm_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][0]
        baux.append(bd)
        if baux[-2]==bd:
            gkm = gkm_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][1]
            bkm = bkm_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][1]


        soma = soma + gkm*v_total[barra]**2 - v_total[barra]*v_total[bd]*(gkm*sep.cos(theta_total[barra]-theta_total[bd])-bkm*sep.sin(theta_total[barra]-theta_total[bd]))


    for bd in barra_origem_linhas[barra_destino_linhas==barra]:

        gkm = gkm_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][0]
        bkm = bkm_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][0]
        baux.append(bd)
        if baux[-2]==bd:
            gkm = gkm_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][1]
            bkm = bkm_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][1]
        soma = soma + gkm*v_total[barra]**2 - v_total[barra]*v_total[bd]*(gkm*sep.cos(theta_total[barra]-theta_total[bd])-bkm*sep.sin(theta_total[barra]-theta_total[bd]))

    for bd in barra_destino_trafos[barra_origem_trafos==barra]:
        gkm = gkm_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][0]
        bkm = bkm_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][0]
        tkm = taps[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][0]/100
        baux.append(bd)
        if baux[-2]==bd:
            gkm = gkm_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][1]
            bkm = bkm_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][1]
            tkm = taps[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][1]/100

        soma = soma + gkm*v_total[barra]**2/(tkm**2) - v_total[barra]*v_total[bd]*(1/tkm)*(gkm*sep.cos(theta_total[barra]-theta_total[bd])-bkm*sep.sin(theta_total[barra]-theta_total[bd]))


    for bd in barra_origem_trafos[barra_destino_trafos==barra]:

        gkm = gkm_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][0]
        bkm = bkm_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][0]
        tkm = taps[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][0]/100
        baux.append(bd)
        if baux[-2]==bd:
            gkm = gkm_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][1]
            bkm = bkm_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][1]
            tkm = taps[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][1]/100

        soma = soma + gkm*v_total[barra]**2 - v_total[barra]*v_total[bd]*(1/tkm)*(gkm*sep.cos(theta_total[barra]-theta_total[bd])-bkm*sep.sin(theta_total[barra]-theta_total[bd]))


    soma = soma + p_demanda[int(barra)] - p_gerada[int(barra)] + p_gs[int(barra)]*v_total[int(barra)]**2
    return soma



def calcula_fluxo_q(sep,barra,v_total,theta_total,barra_destino_trafos,barra_origem_trafos,barra_destino_linhas,barra_origem_linhas,gkm_linhas,gkm_trafos,bkm_linhas,bkm_trafos,taps,q_demanda,q_gerada,q_bs,bc_linhas, bc_trafos):

    soma = 0

    barra = int(barra)

    baux = []
    baux.append(33333333)
    baux.append(33333331)


    for bd in barra_destino_linhas[barra_origem_linhas==barra]:

        gkm = gkm_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][0]
        bkm = bkm_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][0]
        bc = bc_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][0]
        baux.append(bd)

        if baux[-2]==bd:
            gkm = gkm_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][1]
            bkm = bkm_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][1]
            bc = bc_linhas[(barra_origem_linhas == barra) & (barra_destino_linhas == bd)][1]


        soma = soma + -(-bkm+bc/2)*v_total[barra]**2+v_total[barra]*v_total[bd]*(-bkm*sep.cos(theta_total[barra]-theta_total[bd])-gkm*sep.sin(theta_total[barra]-theta_total[bd]))

    for bd in barra_origem_linhas[barra_destino_linhas==barra]:

        gkm = gkm_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][0]
        bkm = bkm_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][0]
        bc = bc_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][0]
        baux.append(bd)
        if baux[-2]==bd:
            gkm = gkm_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][1]
            bkm = bkm_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][1]
            bc = bc_linhas[(barra_origem_linhas == bd) & (barra_destino_linhas == barra)][1]

        soma = soma + -(-bkm+bc/2)*v_total[barra]**2+v_total[barra]*v_total[bd]*(-bkm*sep.cos(theta_total[barra]-theta_total[bd])-gkm*sep.sin(theta_total[barra]-theta_total[bd]))
    for bd in barra_destino_trafos[barra_origem_trafos==barra]:
        gkm = gkm_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][0]
        bkm = bkm_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][0]
        tkm = taps[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][0]/100
        bc = bc_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][0]
        baux.append(bd)
        if baux[-2]==bd:
            gkm = gkm_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][1]
            bkm = bkm_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][1]
            tkm = taps[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][1]/100
            bc = bc_trafos[(barra_origem_trafos == barra) & (barra_destino_trafos == bd)][1]

        soma = soma + -(-bkm+bc/2)*(v_total[barra]**2)/(tkm**2) + v_total[barra]*v_total[bd]*(1/tkm)*(-bkm*sep.cos(theta_total[barra]-theta_total[bd])-gkm*sep.sin(theta_total[barra]-theta_total[bd]))

    for bd in barra_origem_trafos[barra_destino_trafos==barra]:


        gkm = gkm_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][0]
        bkm = bkm_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][0]
        tkm = taps[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][0]/100
        bc = bc_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][0]
        baux.append(bd)
        if baux[-2]==bd:
            gkm = gkm_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][1]
            bkm = bkm_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][1]
            tkm = taps[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][1]/100
            bc = bc_trafos[(barra_origem_trafos == bd) & (barra_destino_trafos == barra)][1]

        soma = soma + -(-bkm+bc/2)*(v_total[barra]**2) + v_total[barra]*v_total[bd]*(1/tkm)*(-bkm*sep.cos(theta_total[barra]-theta_total[bd])-gkm*sep.sin(theta_total[barra]-theta_total[bd]))


    return  q_gerada[int(barra)] - q_demanda[int(barra)] + (q_bs[int(barra)]/100)*v_total[int(barra)]**2 - soma



def Sequential_Quadratic_Programming_Branch_and_Bound(path, n_shunts_vizinhos, n_taps_vizinhos):
    
    
    dados_barra = pd.read_csv(path+'dadosdebarra.csv', header=None)
    dados_barra.columns = ['bus_number','bus_type','pd','qd','gs','bs','area','vm','va','base','zone','max_vm','min_vm']
    dados_gerador = pd.read_csv(path+'dadosgerador.csv', header=None)
    dados_gerador = dados_gerador.iloc[:,[0,1,2,3,4,5,8,9]]
#     dados_custo = pd.read_csv(path+'dadoscusto.csv', header=None)
#     a_k = pd.read_csv('a_k.csv', header=None).values.reshape(-1)
#     b_k = pd.read_csv('b_k.csv', header=None).values.reshape(-1)
#     e_k = pd.read_csv('e_k.csv', header=None).values.reshape(-1)
#     f_k  = pd.read_csv('f_k.csv', header=None).values.reshape(-1)
    dados_gerador.columns = ['bus_number','pg','qg','q_max','q_min','vm','p_max','p_min']
    pg_min = dados_gerador['p_min'].values
    pg_max = dados_gerador['p_max'].values
    bus_n_g = dados_gerador['bus_number'].values
    dados_ramo =pd.read_csv(path+'dadosramo.csv', header=None)
    dados_ramo = dados_ramo.iloc[:,[0,1,2,3,4,8]]
    dados_ramo.columns = ['from','to', 'resistance', 'reactance', 'line_susceptance','tap']
    dados_ramo_linhas = dados_ramo[dados_ramo['tap']==0].copy()
    dados_ramo_trafos = dados_ramo[dados_ramo['tap']!=0].copy()
    barra_origem_linhas = dados_ramo_linhas['from'].to_numpy(dtype=int)
    barra_destino_linhas = dados_ramo_linhas['to'].to_numpy(dtype=int)
    barra_origem_trafos = dados_ramo_trafos['from'].to_numpy(dtype=int)
    barra_destino_trafos = dados_ramo_trafos['to'].to_numpy(dtype=int)
    gkm_linhas = dados_ramo_linhas['resistance'].values/(dados_ramo_linhas['resistance'].values**2+dados_ramo_linhas['reactance'].values**2)
    bkm_linhas = dados_ramo_linhas['reactance'].values/(dados_ramo_linhas['resistance'].values**2+dados_ramo_linhas['reactance'].values**2)
    gkm_trafos = dados_ramo_trafos['resistance'].values/(dados_ramo_trafos['resistance'].values**2+dados_ramo_trafos['reactance'].values**2)
    bkm_trafos = dados_ramo_trafos['reactance'].values/(dados_ramo_trafos['resistance'].values**2+dados_ramo_trafos['reactance'].values**2)
    bc_linhas =  dados_ramo_linhas['line_susceptance'].values
    bc_trafos =  dados_ramo_trafos['line_susceptance'].values
    taps = dados_ramo_trafos['tap'].values
    v_total = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    theta_total = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    p_gerado = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    p_demanda = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    q_demanda = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    p_gerada = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    max_p = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    min_p = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    q_gerada = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    q_gerada_max = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    q_gerada_min = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    p_gs = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    q_bs = np.arange(start=1,stop=dados_barra.bus_number.max()+2,step=1)*0.0
    sbase = 100
    for bus in dados_barra.bus_number:
        v_total[bus] = dados_barra['vm'][dados_barra['bus_number']==bus].to_numpy(dtype=float)
        theta_total[bus] = np.radians(dados_barra['va'][dados_barra['bus_number']==bus].to_numpy(dtype=float))
        p_demanda[bus] = dados_barra['pd'][dados_barra['bus_number']==bus].to_numpy(dtype=float)/sbase
        q_demanda[bus] = dados_barra['qd'][dados_barra['bus_number']==bus].to_numpy(dtype=float)/sbase
        p_gs[bus] = dados_barra['gs'][dados_barra['bus_number']==bus].to_numpy(dtype=float)/sbase
        q_bs[bus] = dados_barra['bs'][dados_barra['bus_number']==bus].to_numpy(dtype=float)/sbase

    for bus in dados_gerador.bus_number:

        q_gerada[bus] = dados_gerador['qg'][dados_gerador['bus_number']==bus].to_numpy(dtype=float)/sbase
        q_gerada_max[bus] = dados_gerador['q_max'][dados_gerador['bus_number']==bus].to_numpy(dtype=float)/sbase
        q_gerada_min[bus] = dados_gerador['q_min'][dados_gerador['bus_number']==bus].to_numpy(dtype=float)/sbase
        p_gerada[bus] = dados_gerador['pg'][dados_gerador['bus_number']==bus].to_numpy(dtype=float)/sbase
        max_p[bus] = dados_gerador['p_max'][dados_gerador['bus_number']==bus].to_numpy(dtype=float)/sbase
        min_p[bus] = dados_gerador['p_min'][dados_gerador['bus_number']==bus].to_numpy(dtype=float)/sbase



    sep = GEKKO(remote=False)
    
    tensoes = np.zeros(len(v_total))
    angulos = np.zeros(len(v_total))

    tensoes = tensoes.tolist()
    angulos = angulos.tolist()

    for bus in dados_barra.bus_number:
        for i in range(len(v_total)):

            if i==bus:
                tensoes[bus] = sep.Var(v_total[bus],0.9,1.1)

                if bus == 7049:

                    angulos[bus] = 0

                else: angulos[bus] = sep.Var(theta_total[bus],-np.pi,np.pi)

    
    shunt = np.zeros(len(q_bs)).tolist()
    # shunt[117]= sep.Var(q_bs[117]*100,0*100,4.5*100, integer=True)#
    # shunt[120]= sep.Var(q_bs[120]*100,0*100, 0.59*100, integer=True)#
    # shunt[154]= sep.Var(q_bs[154]*100,0*100,0.39*100, integer=True)#
    # shunt[164]= sep.Var(q_bs[164]*100,-4.5*100,0*100, integer=True)#
    # shunt[166]= sep.Var(q_bs[166]*100,-4.5*100,0*100, integer=True)#
    # shunt[173]= sep.Var(q_bs[173]*100,0*100,0.59*100, integer=True)#
    # shunt[179]= sep.Var(q_bs[179]*100,0*100,0.59*100, integer=True)#
    # shunt[190]= sep.Var(q_bs[190]*100,-2.5*100,0*100, integer=True)#
    # shunt[231]= sep.Var(q_bs[231]*100,-4.5*100,0*100, integer=True)#
    # shunt[238]= sep.Var(q_bs[238]*100,-4.5*100,0*100, integer=True)#
    # shunt[240]= sep.Var(q_bs[240]*100,-1.5*100,0*100, integer=True)#
    # shunt[248]= sep.Var(q_bs[248]*100,0*100,0.59*100, integer=True)#
    # shunt[9003]= sep.Var(q_bs[9003]*100,0*100,0.15*100, integer=True)#
    # shunt[9034]= sep.Var(q_bs[9034]*100,0*100,0.15*100, integer=True)#
    
    # shunt=q_bs.copy()*100

    a = n_shunts_vizinhos

    aux_u = int(q_bs[117]*100)+a
    aux_l = int(q_bs[117]*100)-a

    if q_bs[117]*100+a>4.5*100:

        aux_u = 4.5*100

    if q_bs[117]*100-a<0:

        aux_l=0


    shunt[117]= sep.Var(q_bs[117]*100,aux_l,aux_u, integer=True)



    #######################################




    aux_u = int(q_bs[120]*100)+a
    aux_l = int(q_bs[120]*100)-a

    if q_bs[120]*100+a>0.59*100:

        aux_u = 0.59*100

    if q_bs[120]*100-a<0:

        aux_l=0



    shunt[120]= sep.Var(q_bs[120]*100,aux_l,aux_u, integer=True)#


    #######################################

    aux_u = int(q_bs[154]*100)+a
    aux_l = int(q_bs[154]*100)-a

    if q_bs[154]*100+a>0.39*100:

        aux_u = 0.39*100

    if q_bs[154]*100-a<0:

        aux_l=0

    shunt[154]= sep.Var(q_bs[154]*100,aux_l,aux_u, integer=True)#

    #######################################

    aux_u = int(q_bs[164]*100)+a
    aux_l = int(q_bs[164]*100)-a

    if q_bs[164]*100+a>0*100:

        aux_u = 0*100

    if q_bs[164]*100-a<-4.5*100:

        aux_l=-4.5*100

    shunt[164]= sep.Var(q_bs[164]*100,aux_l,aux_u, integer=True)#

    #######################################

    aux_u = int(q_bs[166]*100)+a
    aux_l = int(q_bs[166]*100)-a

    if q_bs[166]*100+a>0*100:

        aux_u = 0*100

    if q_bs[166]*100-a<-4.5*100:

        aux_l=-4.5*100

    shunt[166]= sep.Var(q_bs[166]*100,aux_l,aux_u, integer=True)#

    #######################################

    aux_u = int(q_bs[173]*100)+a
    aux_l = int(q_bs[173]*100)-a

    if q_bs[173]*100+a>0.59*100:

        aux_u = 0.59*100

    if q_bs[173]*100-a<0:

        aux_l=0



    shunt[173]= sep.Var(q_bs[173]*100,aux_l,aux_u, integer=True)#

    #######################################

    aux_u = int(q_bs[179]*100)+a
    aux_l = int(q_bs[179]*100)-a

    if q_bs[179]*100+a>0.59*100:

        aux_u = 0.59*100

    if q_bs[179]*100-a<0:

        aux_l=0



    shunt[179]= sep.Var(q_bs[179]*100,aux_l,aux_u, integer=True)#

    #######################################

    aux_u = int(q_bs[190]*100)+a
    aux_l = int(q_bs[190]*100)-a

    if q_bs[190]*100+a>0*100:

        aux_u = 0*100

    if q_bs[190]*100-a<-2.5*100:

        aux_l=-2.5*100



    shunt[190]= sep.Var(q_bs[190]*100,aux_l,aux_u, integer=True)#

    #######################################

    aux_u = int(q_bs[231]*100)+a
    aux_l = int(q_bs[231]*100)-a

    if q_bs[231]*100+a>0*100:

        aux_u = 0*100

    if q_bs[231]*100-a<-4.5*100:

        aux_l=-4.5*100


    shunt[231]= sep.Var(q_bs[231]*100,aux_l,aux_u, integer=True)#


    #######################################

    aux_u = int(q_bs[238]*100)+a
    aux_l = int(q_bs[238]*100)-a

    if q_bs[238]*100+a>0*100:

        aux_u = 0*100

    if q_bs[238]*100-a<-4.5*100:

        aux_l=-4.5*100


    shunt[238]= sep.Var(q_bs[238]*100,aux_l,aux_u, integer=True)#

    #######################################

    aux_u = int(q_bs[240]*100)+a

    aux_l = int(q_bs[240]*100)-a

    if q_bs[240]*100+a>0*100:

        aux_u = 0*100

    if q_bs[240]*100-a<-1.5*100:

        aux_l=-1.5*100


    shunt[240]= sep.Var(q_bs[240]*100,aux_l,aux_u, integer=True)#


    #######################################

    aux_u = int(q_bs[248]*100)+a
    aux_l = int(q_bs[248]*100)-a

    if q_bs[248]*100+a>0.59*100:

        aux_u = 0.59*100

    if q_bs[248]*100-a<0:

        aux_l=0



    shunt[248]= sep.Var(q_bs[248]*100,aux_l,aux_u, integer=True)#



    #######################################

    aux_u = int(q_bs[9003]*100)+a
    aux_l = int(q_bs[9003]*100)-a

    if q_bs[9003]*100+a>0.15*100:

        aux_u = 0.15*100

    if q_bs[9003]*100-a<0:

        aux_l=0


    #######################################

    aux_u = int(q_bs[9034]*100)+a
    aux_l = int(q_bs[9034]*100)-a

    if q_bs[9034]*100+a>0.15*100:

        aux_u = 0.15*100

    if q_bs[9034]*100-a<0:

        aux_l=0


    shunt[9034]= sep.Var(q_bs[9034]*100,aux_l,aux_u, integer=True)#

    pg = p_gerada.copy()
    pg = pg.tolist()

    pg[7049] = sep.Var(p_gerada[7049],0,2399.005/100)

    # pg = np.array(pg)

    qg = q_gerada.copy()
    qg = qg.tolist()

    for bus in dados_barra.bus_number:
        
        for i in range(len(qg)):
            
            if i==bus:
                
                qg[bus] = sep.Var(q_gerada_max[bus],q_gerada_min[bus],q_gerada_max[bus])

    trafos = []

    b = n_taps_vizinhos

    for val in taps:

        aux_u = val*100+b
        aux_l = val*100-b

        if val*100+b>110:

            aux_u = 110

        if val*100-b<90:

            aux_l=90


        trafos.append(sep.Var(val*100, aux_l, aux_u, integer=True))
        
    total_bus = dados_barra.bus_number.to_numpy(dtype=float)
    trafos = np.array(trafos)
    shunt = np.array(shunt)
    
    # trafos = np.array(taps)*100
    
    eq= calcula_perdas(sep,barra_origem_linhas,barra_destino_linhas,barra_origem_trafos,barra_destino_trafos,gkm_linhas,gkm_trafos,tensoes,angulos,trafos)

    sep.Obj(sep.sum(eq))
    

    for barra in total_bus:

        sep.Equation(calcula_fluxo_p(sep,barra,tensoes,angulos,barra_destino_trafos,barra_origem_trafos,barra_destino_linhas,barra_origem_linhas,gkm_linhas,gkm_trafos,bkm_linhas,bkm_trafos,trafos,pg,p_demanda,p_gs)==0)


    for barra in total_bus:

        sep.Equation(calcula_fluxo_q(sep,barra,tensoes,angulos,barra_destino_trafos,barra_origem_trafos,barra_destino_linhas,barra_origem_linhas,gkm_linhas,gkm_trafos,bkm_linhas,bkm_trafos,trafos,q_demanda,qg,shunt,bc_linhas, bc_trafos)==0)
        

    sep.solver_options = ['minlp_gap_tol 1',\
                          'minlp_maximum_iterations 100000',\
                          'minlp_max_iter_with_int_sol 500',\
                          'minlp_branch_method 1',\
                          'nlp_maximum_iterations 1000',\
                          'minlp_integer_tol 5.0e-5']


    sep.options.SOLVER=1

    sep.solve(disp=True)
    

    qfilter = dados_gerador['bus_number'].values
    
    t = []
    s = []
    v = []
    a = []
    q = []
    
    for val in trafos:
        
        try: t.append(val[0])
            
        except: t.append(val)
    
    for val in shunt:
        
        try: s.append(val[0])
            
        except: s.append(val)

    for val in tensoes:
        
        try: v.append(val[0])
            
        except: pass

    i=0
    
    for val in angulos:
        
        try: a.append(val[0])
            
        except: 
            
            if i ==0:
                
                a.append(val)
                
            else: pass
            
        i=i+1
        
    
    for val in qg:
        
        try: q.append(val[0])
            
        except: 
              
            q.append(val)
            
            
    s=np.array(s)

    q=np.array(q)
    
    
    filt = [117,120,154,164,166,173,179,190,231,238,240,248,9003,9034]
    
    
    results=    {
                'funcao_objetivo':[sep.options.OBJFCNVAL],
                'iterations':[sep.options.ITERATIONS],
                'voltages':v,
                'angles':a,
                'shunts':s[filt],
                'trafos':t,
                'reactive_power':q[qfilter]
                }
            
    return results



if __name__ == "__main__":

    # path = 'init/'



    path = f'CLPSO/{number}/'
    dictionary = Sequential_Quadratic_Programming_Branch_and_Bound(path, 2,2)
    pd.DataFrame(dictionary['funcao_objetivo']).to_csv('obj.csv')
    pd.DataFrame(dictionary['iterations']).to_csv('iterations.csv')
    pd.DataFrame(dictionary['voltages']).to_csv('voltages.csv')
    pd.DataFrame(dictionary['angles']).to_csv('angles.csv')
    pd.DataFrame(dictionary['angles']).to_csv('shunts.csv')
    pd.DataFrame(dictionary['trafos']).to_csv('trafos.csv')
    pd.DataFrame(dictionary['angles']).to_csv('reactive_power.csv')