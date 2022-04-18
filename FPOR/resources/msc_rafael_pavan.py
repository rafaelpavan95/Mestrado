import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandapower as pp
import time
import random
from pandapower.networks import case14, case_ieee30, case118, case300, case4gs
import tabulate
import numba
from numba import njit
from gekko import GEKKO


def inicializa_sep(sep):
    
    if len(sep.bus)==14:
        
        sep.bus['min_vm_pu'] = 0.95

        sep.bus['max_vm_pu'] = 1.05              
      
        pp.runpp(sep, init = 'results', tolerance_mva = 1e-6, trafo_model='pi')
        
        
    if len(sep.bus)==30:
        
        sep.bus['min_vm_pu'] = 0.95
        
        sep.bus['max_vm_pu'] = 1.05
        
        sep.gen['max_q_mvar']=np.array([50,40,40,24,24])
        
        sep.gen['min_q_mvar']=np.array([-40,-40,-10,-6,-6])
        
        sep.ext_grid['max_q_mvar'] = 10
        
        sep.ext_grid['min_q_mvar'] = 0
            
        pp.runpp(sep, init = 'results', tolerance_mva = 1e-6)
                
    if len(sep.bus)==118:
        
        sep.bus['min_vm_pu'] = 0.94
        
        sep.bus['max_vm_pu'] = 1.06
        
        sep.ext_grid['va_degree'] = 0

        pp.runpp(sep, init = 'results', tolerance_mva = 1e-6, trafo_model='pi')
        
                    
    if len(sep.bus)==300:
        
        sep.bus['min_vm_pu'] = 0.9
        sep.bus['max_vm_pu'] = 1.1

        pp.runpp(sep, init = 'results', tolerance_mva = 1e-6, trafo_model='pi')
            
    voltages_init = sep.gen['vm_pu'].to_numpy()

    tap_pos = sep.trafo[~pd.isnull(sep.trafo['tap_pos'])]['tap_pos'].to_numpy(dtype=np.float64)

    tap_neutral = sep.trafo[~pd.isnull(sep.trafo['tap_neutral'])]['tap_neutral'].to_numpy(dtype=np.float64)

    tap_step_percent = sep.trafo[~pd.isnull(sep.trafo['tap_step_percent'])]['tap_step_percent'].to_numpy(dtype=np.float64)       

    valor_pu_tap = (tap_pos-tap_neutral)*(tap_step_percent/100) + 1

    valor_bshunt = (sep.shunt['q_mvar']/(-100)).to_numpy()

    zeros = np.array([0,0,0,0,0,0])

    valor_inicial = np.expand_dims(np.concatenate((voltages_init, valor_pu_tap, valor_bshunt,zeros), axis = None), 0)

    return valor_inicial



def matriz_condutancia(sep,relatorio=True):
    
    '''
    
    Calcula a matriz de condutâncias de linha, retornando apenas a parte triangular superior.
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    
    Parâmetros
    ----------   
    sep : sistema elétrico de potência carregado pelo pandapower.
    relatorio : caso relatorio = True, retorna relatório informando barras de origem e destino das linhas, resistências (pu), reatâncias (pu).
                caso relatorio = False, retorna apenas a parte triangular superior da matriz de condutâncias.
    Retorno
    -------    
    matriz_g: matriz de condutâncias entre barras com triângulo inferior zerado.
    
    Observações:
    ------------
    
    Potência Aparente de Base = 100 MVA
    
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    
    
    '''
    
    sep.line=sep.line.sort_index()
   
    sep.bus=sep.bus.sort_index()
    
    vbus = sep.bus.vn_kv.to_numpy(dtype=np.float64)
    
    zbase = np.power(np.multiply(vbus,1000), 2)/100e6
    
    # Inicializa Matriz Zerada
    
    matriz_z = np.zeros((9,len(sep.line.index.ravel())),dtype=np.float64)
    
    matriz_g = np.zeros((sep.bus.name.count(),sep.bus.name.count()), dtype=np.float64)
    
    g = np.zeros(len(sep.line.index.ravel()),dtype=np.float64)
    
    # Pega Valores de Barra Origem e Destino das Linhas
        
    matriz_z[0,:]=sep.line.from_bus
    
    matriz_z[1,:]=sep.line.to_bus
    
    
    for i in range(len(sep.line.index.ravel())):    
    
        matriz_z[2,i] = sep.line.r_ohm_per_km[i]/zbase[int(matriz_z[0,i])]
        matriz_z[3,i] = sep.line.x_ohm_per_km[i]/zbase[int(matriz_z[0,i])]
    
    # Calcula Condutâncias
    
    g = np.array(np.divide(matriz_z[2,:], np.power(matriz_z[2,:],2)+np.power(matriz_z[3],2)))
    z = np.sqrt(np.power(matriz_z[2,:],2) + np.power(matriz_z[3,:],2))
    b = np.array(np.divide(matriz_z[3,:], np.power(matriz_z[2,:],2)+np.power(matriz_z[3],2)))
    matriz_z[4,:]=g
    
    vo = []
    vd = []
    to = []
    td = []


    for bus in matriz_z[0,:]:

        vo.append(sep.res_bus['vm_pu'][sep.res_bus.index==bus].to_numpy(dtype=np.float64))
        to.append(sep.res_bus['va_degree'][sep.res_bus.index==bus].to_numpy(dtype=np.float64))


    for bus in matriz_z[1,:]:

        vd.append(sep.res_bus['vm_pu'][sep.res_bus.index==bus].to_numpy(dtype=np.float64))
        td.append(sep.res_bus['va_degree'][sep.res_bus.index==bus].to_numpy(dtype=np.float64))
    
    matriz_z[5,:] = vo
    matriz_z[6,:] = to
    matriz_z[7,:] = vd
    matriz_z[8,:] = td
    
    # Gera Matriz
    
    for posicao in range(len(sep.line.index.ravel())):
        
        matriz_g[matriz_z[0,posicao].astype(np.int),matriz_z[1,posicao].astype(np.int)] = g[posicao]
        
    
    if relatorio==True:
    
        tabela = np.zeros((len(sep.line.index.ravel()),7))
        tabela[:,0] = matriz_z[0,:]
        tabela[:,1] = matriz_z[1,:]
        tabela[:,2] = matriz_z[2,:]
        tabela[:,3] = matriz_z[3,:]
        tabela[:,4] = z
        tabela[:,5] = g
        tabela[:,6] = b

        table = tabulate.tabulate(tabela, headers = ['Barra de Origem', 'Barra de Destino','R (pu)','Xl (pu)','Z (pu)', 'G (pu)','B (pu)'], tablefmt="psql")
        print(table)
        
        if len(sep.bus)==14:

            sns.heatmap(matriz_g+matriz_g.T,annot=True,fmt='.6g',cmap='jet')
            plt.xlabel('Barra Origem')
            plt.ylabel('Barra Destino')
            plt.title('Matriz de Condutâncias de Linha Completa')

        
    if relatorio==False:
        
        return matriz_g, matriz_z



def coleta_dados_vbus(sep,relatorio=True):
    
   
    '''
    
    Coleta os Dados de Tensões e Limites Superiores e Inferiores das Barras do Sistema.
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    
    Parâmetros
    ----------
    sep : sistema elétrico de potência carregado pelo pandapower.
    relatorio : caso relatorio = True, retorna relatório informando, tensões, ângulos e limites.
                caso relatorio = False, retorna apenas as tensões, ângulos e limtes
    
    Retorno
    ----------
    vbus : vetor de tensões [pu] das barras em ordem crescente do número da barra
    theta : vetor de ângulo de tensões [°]
    v_lim_superior : 
    
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
        
    '''
    
    sep.res_bus=sep.res_bus.sort_index()
    
    sep.bus=sep.bus.sort_index()
      
    vbus = sep.res_bus['vm_pu'].to_numpy(dtype=np.float64)
    
    theta = sep.res_bus['va_degree'].to_numpy(dtype=np.float64)
    
    v_lim_superior = sep.bus["max_vm_pu"].to_numpy(dtype=np.float32)
    
    v_lim_inferior = sep.bus["min_vm_pu"].to_numpy(dtype=np.float32)

    
    
    if relatorio==True:
        
        tabela = np.zeros((len(vbus),4))
        tabela[:,0] = vbus
        tabela[:,1] = theta
        tabela[:,2] = v_lim_inferior
        tabela[:,3] = v_lim_superior

        table = tabulate.tabulate(tabela, headers = ['Tensões nas Barras (pu)', 'Ângulos das Barras (°)','Limites Inferiores','Limites Superiores'], tablefmt="psql")
        print(table)
    
        sns.scatterplot(x=np.arange(0,len(vbus),1),y=vbus,color='blue',label='Módulo da Tensão',s=75)
        sns.lineplot(x=np.arange(0,len(vbus),1),y=v_lim_superior,color='red',label='Limite Superior',alpha=0.5)
        sns.lineplot(x=np.arange(0,len(vbus),1),y=v_lim_inferior,color='orange',label='Limite Inferior',alpha=0.5)
        plt.title('Módulo da Tensão por Barra do Sistema')
        plt.xlabel('Barra do Sistema')
        plt.ylabel('Tensão [pu]')
        plt.legend(loc='best')
        plt.figure(figsize=(16,10))
        sns.scatterplot(x=np.arange(0,len(theta),1),y=theta,color='green',label='Ângulo da Tensão',s=75)
        plt.title('Ângulo da Tensão por Barra do Sistema')
        plt.xlabel('Barra do Sistema')
        plt.ylabel('Theta [°]')
        plt.legend(loc='best')
        
    
    if relatorio==False:
        
        return vbus, theta, v_lim_superior, v_lim_inferior
    
    
    ##################################################################################################################################################################################


def coleta_dados_gen(sep,relatorio=True):
       
    '''
    
    Coleta os Dados de Tensões, Potências Ativa e Reativa e Seus Respectivos Limites Superiores e Inferiores de geração.
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    
    Parâmetros
    ----------
    sep : sistema elétrico de potência carregado pelo pandapower.
    relatorio : caso relatorio = True, retorna relatório informando, limites, potências e gráficos.
                caso relatorio = False, retorna apenas as tensões, ângulos, potências e limites.
    
    Retorno
    ----------
    vgen : vetor de tensões [pu] das barras de geração
    theta : vetor de ângulo de tensões [°] das barras de geração
    p_lim_superior : Limite Superior de Potência Ativa (pu)
    p_lim_inferior : Limite Inferior de Potência Ativa (pu)
    q_lim_superior : Limite Superior de Potência Reativa (pu)
    q_lim_inferior : Limite Inferior de Potência Ativa (pu)
    
    Observações:
    - - - - - - -
    
    Potência Aparente de Base : 100 MVA
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
        
    '''
    
    sep.res_gen=sep.res_gen.sort_index()
    
    sep.gen=sep.gen.sort_index()
      
    vgen = sep.res_gen['vm_pu'].to_numpy(dtype=np.float64)
    
    barra = sep.gen['bus'].to_numpy(dtype=np.float64)
    
    thetagen = sep.res_gen['va_degree'].to_numpy(dtype=np.float64)
    
    pgen = sep.res_gen['p_mw'].to_numpy(dtype=np.float64)/100
    
    qgen = sep.res_gen['q_mvar'].to_numpy(dtype=np.float64)/100
    
    p_lim_superior = sep.gen["max_p_mw"].to_numpy(dtype=np.float32)/100
    
    p_lim_inferior = sep.gen["min_p_mw"].to_numpy(dtype=np.float32)/100
    
    q_lim_superior = sep.gen["max_q_mvar"].to_numpy(dtype=np.float32)/100
    
    q_lim_inferior = sep.gen["min_q_mvar"].to_numpy(dtype=np.float32)/100

    
    if relatorio==True:
        
        tabela = np.zeros((len(vgen),6))
        tabela[:,0] = pgen
        tabela[:,1] = p_lim_superior
        tabela[:,2] = p_lim_inferior
        tabela[:,3] = qgen
        tabela[:,4] = q_lim_superior
        tabela[:,5] = q_lim_inferior


        table = tabulate.tabulate(tabela, headers = ['P (pu)','P Lim. Sup. (pu)','P Lim. Inf. (pu)','Q (pu)','Q Lim. Sup. (pu)','Q Lim. Inf. (pu)'], tablefmt="psql")
        print(table)
    

        sns.scatterplot(x=barra,y=qgen,color='blue',label='Potência Gerada',s=75)
        sns.lineplot(x=barra,y=q_lim_superior,color='red',label='Limite Superior',alpha=0.5)
        sns.lineplot(x=barra,y=q_lim_inferior,color='orange',label='Limite Inferior',alpha=0.5)
        plt.title('Potência Reativa Gerada')
        plt.xlabel('Barra do Sistema')
        plt.ylabel('Potência Reativa (pu)')
        plt.legend(loc='best')
        
    
    if relatorio==False:
        
        return vgen, thetagen, pgen, qgen, p_lim_superior, p_lim_inferior, q_lim_superior, q_lim_inferior,barra
    
    

    ################################################################################################################################################################################################


def func_objetivo(vbarra,theta,condutancias,matriz_z,relatorio=True):
    

    '''

    Calcula as perdas nas linhas de transmissão de acordo com as tensões, ângulos das barras e condutâncias de linha.
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

    Parâmetros
    ----------
    vbarra : tensão da barra.
    theta : ângulo da barra.
    condutancias : matriz de condutâncias de linha (triângulo superior)

    caso relatorio = True, retorna relatório informando a matriz de perdas de linha e as perdas totais.
    caso relatorio = False, retorna apenas as perda em pu.
    Retorno
    ----------

    perdas : perdas de potência ativa em pu.


    Observações:
    - - - - - - -

    Potência Aparente de Base : 100 MVA
     _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

    '''
    
    matriz_v = np.zeros((len(vbarra),len(vbarra)), dtype=np.float64)
    
    matriz_theta = np.zeros((len(theta),len(theta)), dtype=np.float64)
    
    for barra in range(len(vbarra)):
        
        matriz_v[:,barra]=vbarra
        matriz_theta[:,barra]=theta
        
    
    soma_v = np.power(matriz_v,2) + np.power(matriz_v.T,2)
    
    subtrai_theta = matriz_theta - matriz_theta.T
    
    cosenotheta=np.cos(np.radians(subtrai_theta))
    
    produto = 2 * np.multiply(np.multiply(matriz_v, matriz_v.T),cosenotheta)
    
    matriz_perdas = np.multiply(condutancias,soma_v-produto) 

    perdas = np.multiply(matriz_z[4,:], np.power(matriz_z[5,:],2)+np.power(matriz_z[7,:],2)-2*np.multiply(np.multiply(matriz_z[5,:],matriz_z[7,:]),np.cos(np.radians(matriz_z[8,:]-matriz_z[6,:]))))
    perdas = np.sum(perdas)
    
    if relatorio == True:
        
        tabela = np.zeros((1,2))
        tabela[:,0] = perdas
        tabela[:,1] = perdas*100
        table = tabulate.tabulate(tabela, headers = ['Perdas Totais Nas Linhas (pu)','Perdas Totais Nas Linhas (MW)'], tablefmt="psql")
        print(table)
        
        if len(vbarra) ==14:
            plt.figure(figsize=(18,10))
            sns.heatmap(100*(matriz_perdas+matriz_perdas.T),annot=True,cmap="jet")
            plt.xlabel('Barra Origem')
            plt.ylabel('Barra Destino')
            plt.title('Matriz de Perdas de Linha Completa [MW]')

     
    else:
    
        return perdas

    

    ################################################################################################################################################################################################


def pen_tensao(vbus, limite_sup, limite_inf,relatorio=True):
    
    """    
    Calcula a parcela de penalização pura (sem constante de multiplicação) referente a violação dos limites de tensão.
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    
    Parâmetros
    ----------   
    vbus : tensões das barras do sistema elétrico.
    limite_sup : limite superior das tensões das barras do sistema elétrico.
    limite_inf : limite inferior das tensões das barras do sistema elétrico.
    
    limite_sup : tensões
    relatorio : caso relatorio = True, retorna penalização e nº de violações 
                caso relatorio = False, retorna apenas o valor de penalização.
    Retorno
    -------    
    penalizacao: somatório da diferença ao quadradado das tensões que ultrapassaram os limites inferiores ou superiores.
    
    Observações:
    ------------
    
    ...
    
    """
    
    
    inferior = vbus - limite_inf
    superior = limite_sup - vbus
    penalizacao = np.sum(np.abs(superior[superior<0]))+np.sum(np.abs(inferior[inferior<0]))
    
    if relatorio == True:
        
        print('Penalização de Tensão:\n')
        print(penalizacao,'\n')
        print('Número de Violações:\n')
        print(len(inferior[inferior<0])+len(superior[superior<0]))
     
    else:
    
        return penalizacao


 ################################################################################################################################################################################################


def pen_ger_reativo(q, limite_sup, limite_inf,sep,relatorio=True):
    
    """    
    Calcula a parcela de penalização pura (sem constante de multiplicação) referente a violação dos limites de geração de reativos.
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    
    Parâmetros
    ----------   
    q : potências reativas das barras de controle de reativo do sistema elétrico.
    limite_sup : limite superior das potências reativas das barras de controle de reativo do sistema elétrico.
    limite_inf : limite superior das potências reativas das barras de controle de reativo do sistema elétrico.
    
    limite_sup : tensões
    relatorio : caso relatorio = True, retorna penalização e nº de violações 
                caso relatorio = False, retorna apenas o valor de penalização.
    Retorno
    -------    
    penalização: somatório da diferença ao quadradado das potências reativas que ultrapassaram os limites inferiores ou superiores.
    
    Observações:
    ------------
    
    ...
    
    """
    
    inferior = limite_inf - q
    superior = limite_sup - q
    
    ext_sup = sep.ext_grid['max_q_mvar'].to_numpy()
    ext_inf = sep.ext_grid['min_q_mvar'].to_numpy()
    
    qext = sep.res_ext_grid['q_mvar'].to_numpy()
    
    inferiorext = ext_inf - qext
    superiorext =  ext_sup - qext
    
    penalizacaoext = np.sum(np.abs(superiorext[superiorext<0]))+np.sum(np.abs(inferiorext[inferiorext>0]))
    penalizacao = np.sum(np.abs(superior[superior<0]))+np.sum(np.abs(inferior[inferior>0]))
    
    
    if relatorio == True:
        
        print('Penalização de Geração de Reativos:\n')
        print(penalizacao+penalizacaoext,'\n')
        print('Número de Violações:\n')
        print(len(inferior[inferior<0])+len(superior[superior<0]))
        
    else:
    
        return penalizacao+penalizacaoext


 ################################################################################################################################################################################################


def coleta_dados_trafo(sep, relatorio=True):
    

    sep.trafo.sort_index()
   
    sep.res_trafo.sort_index()
    
    sep.trafo['tap_pos']=sep.trafo['tap_pos']
    
    n_trafos_controlados = sep.trafo['tap_pos'].count()
    
    carregamento = sep.res_trafo['loading_percent'].to_numpy()/100
    
    tap_pos = sep.trafo[~pd.isnull(sep.trafo['tap_pos'])]['tap_pos'].to_numpy(dtype=np.float64)
    
    tap_neutral = sep.trafo[~pd.isnull(sep.trafo['tap_neutral'])]['tap_neutral'].to_numpy(dtype=np.float64)
    
    tap_step_percent = sep.trafo[~pd.isnull(sep.trafo['tap_step_percent'])]['tap_step_percent'].to_numpy(dtype=np.float64)
        
    
    if len(sep.bus)==14:
        
        step = 0.01
        valores_taps = np.arange(start = 0.9, stop = 1.1, step = step)
        
        
    if len(sep.bus)==30:
        
        step = 0.01
        valores_taps = np.arange(start = 0.9, stop = 1.1, step = step)
        
                
    if len(sep.bus)==118:
        
        step = 0.01
        valores_taps = np.arange(start = 0.9, stop = 1.1, step = step)
           
                    
    if len(sep.bus)==300:
        
        step = 0.01
        valores_taps = np.arange(start = 0.9, stop = 1.1, step = step)
        
        
    if relatorio == True:
        
    
        tap_pos = sep.trafo[~pd.isnull(sep.trafo['tap_pos'])]['tap_pos'].to_numpy(dtype=np.float64)

        tap_neutral = sep.trafo[~pd.isnull(sep.trafo['tap_neutral'])]['tap_neutral'].to_numpy(dtype=np.float64)

        tap_step_percent = sep.trafo[~pd.isnull(sep.trafo['tap_step_percent'])]['tap_step_percent'].to_numpy(dtype=np.float64)
        
        valor_percentual= (tap_pos-tap_neutral)*(tap_step_percent/100) + 1
        
      
        if len(sep.bus)==14:

            plt.figure(figsize=(20,10))
            sns.scatterplot(x=np.arange(0,len(tap_pos)),y=valor_percentual,label='Valor do TAP',color='b',s=75)
            sns.lineplot(x=np.arange(0,len(tap_pos)),y=np.tile([1.12], (len(tap_pos))),label='Limite Máximo do TAP',color='r')
            sns.lineplot(x=np.arange(0,len(tap_pos)),y=np.tile([0.88], (len(tap_pos))),label='Limite Mínimo do TAP',color='orange')
            plt.grid()


        if len(sep.bus)==30:
        
            plt.figure(figsize=(20,10))
            sns.scatterplot(x=np.arange(0,len(tap_pos)),y=valor_percentual,label='Valor do TAP',color='b',s=75)
            sns.lineplot(x=np.arange(0,len(tap_pos)),y=np.tile([1.12], (len(tap_pos))),label='Limite Máximo do TAP',color='r')
            sns.lineplot(x=np.arange(0,len(tap_pos)),y=np.tile([0.88], (len(tap_pos))),label='Limite Mínimo do TAP',color='orange')
            plt.grid()
        
        if len(sep.bus)==118:
        
            plt.figure(figsize=(20,10))
            sns.scatterplot(x=np.arange(0,len(tap_pos)),y=valor_percentual,label='Valor do TAP',color='b',s=75)
            sns.lineplot(x=np.arange(0,len(tap_pos)),y=np.tile([1.12], (len(tap_pos))),label='Limite Máximo do TAP',color='r')
            sns.lineplot(x=np.arange(0,len(tap_pos)),y=np.tile([0.88], (len(tap_pos))),label='Limite Mínimo do TAP',color='orange')
            plt.grid()
        
        print('Carregamento do Trafo (pu):\n')
        print(carregamento,'\n')
        print('Número de Trafos com TAP Controlado:\n')
        print(n_trafos_controlados,'\n')
        print('Valores dos TAPs:\n')
        print(valor_percentual,'\n')
        
        
    else:
        
        return tap_pos, tap_neutral, tap_step_percent,valores_taps



 ################################################################################################################################################################################################


def pen_trafo(linha,n_tap,n_vgen):
    
    '''    
    
    
    Valores dos TAPs Retirados de:
    
    - REFORMULAÇÃO DAS RESTRIÇÕESDE COMPLEMENTARIDADE EM PROBLEMAS DE FLUXO DE POTÊNCIA ÓTIMO
      Marina Valença Alencar - Dissertação de Mestrado

    - FUNÇÕES PENALIDADE PARA O TRATAMENTO DAS VARIÁVEIS DISCRETAS DO PROBLEMA DE FLUXO DE POTÊNCIA ÓTIMO REATIVO
      Daisy Paes Silva - Dissertação de Mestrado
    

    ''' 
    """    
   
    Calcula a penalização senoidal para taps não discretos.
    _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
    
    Parâmetros
    ----------   
    linha : linha da partícula.
    n_tap : número de taps.
    n_vgen : número de geradores.
    
    Retorno
    -------    
    linha : linha da partícula com os valores da penalização do trafo atualizados.
    
    
    Observações:
    ------------
    
    ...
    
    """
    
    step = 0.01

    linha[-3] = np.sum(np.square(np.sin((linha[n_vgen:n_vgen+n_tap]*np.pi/step))))
    
    return linha


 ################################################################################################################################################################################################


def coleta_dados_bshunt(sep):
    
    '''    
    
    
    Valores dos Shunt Retirados de:
    
    - REFORMULAÇÃO DAS RESTRIÇÕESDE COMPLEMENTARIDADE EM PROBLEMAS DE FLUXO DE POTÊNCIA ÓTIMO
      Marina Valença Alencar - Dissertação de Mestrado

    - FUNÇÕES PENALIDADE PARA O TRATAMENTO DAS VARIÁVEIS DISCRETAS DO PROBLEMA DE FLUXO DE POTÊNCIA ÓTIMO REATIVO
      Daisy Paes Silva - Dissertação de Mestrado
          

    ''' 
    
    ieee118 = np.arange(0.00,0.45,0.001)
    ieee30 = np.arange(0.00,0.35,0.001)
    ieee300 = np.arange(start=-3.25,stop=3.25,step=0.1)
    
    bus = sep.shunt['bus'].sort_values().to_numpy()
    sep.shunt=sep.shunt.sort_index()
  
    
    if len(sep.bus)==14:
        
        bsh = np.array([[0,0.01,0.02,0.03,0.04,0.05]],dtype=object)
        
        
    if len(sep.bus)==30:
        
        bsh = np.array([[0,0.01, 0.02, 0.03, 0.04, 0.05],[0,0.01, 0.02, 0.03, 0.04, 0.05]],dtype=object)
        
                
    if len(sep.bus)==118:
        
        bsh = np.array([np.arange(start=-0.40, stop=0,step=0.01), 
                       np.arange(start=0, stop=0.14,step=0.01), 
                       np.arange(start=-0.25, stop=0,step=0.01), 
                       np.arange(start=0, stop=0.1,step=0.01), 
                       np.arange(start=0, stop=0.1,step=0.01),
                       np.arange(start=0, stop=0.1,step=0.01),
                       np.arange(start=0, stop=0.15,step=0.01),
                       np.arange(start=0, stop=0.12,step=0.01),                       
                       np.arange(start=0, stop=0.2,step=0.01),
                       np.arange(start=0, stop=0.2,step=0.01),
                       np.arange(start=0, stop=0.1,step=0.01),
                       np.arange(start=0, stop=0.2,step=0.01),
                       np.arange(start=0, stop=0.06,step=0.01),
                       np.arange(start=0, stop=0.06,step=0.01)],dtype=object)

#         bsh = np.array([[-0.40, 0], 
#                        [0, 0.06, 0.07, 0.13, 0.14, 0.20], 
#                        [-0.25, 0], 
#                        [0, 0.10], 
#                        [0, 0.10], 
#                         [0, 0.10],
#                         [0, 0.15],
#                         [0.08, 0.12, 0.20],
#                         [0, 0.10, 0.20],
#                         [0, 0.10, 0.20],
#                         [0, 0.10, 0.20],
#                         [0, 0.10, 0.20],
#                         [0, 0.06, 0.07, 0.13, 0.14, 0.20],
#                         [0, 0.06, 0.07, 0.13, 0.14, 0.20]],dtype=object)
        
              

    if len(sep.bus)==300:
        
        bsh = np.array([[0,1.5,3,4.5], #95
                [0, 0.15, 0.30, 0.60], #98
                [0,0.15,0.30,0.45], #132
                [-4.5,-3,-1.5,0], #142
                [-4.5,-3,-1.5,0], #144
                [0, 0.15,0.30,0.45,0.60], #151
                [0, 0.15,0.30,0.45,0.60], #157
                [-4.5,-3,-1.5,0], #168
                [-4.5,-3,-1.5,0], #209
                [-4.5,-3,-1.5,0],#216
                [-4.5,-3,-1.5,0], #218
                [0, 0.15,0.30,0.45,0.60], #226
                [0, 0.05,0.01,0.15], #267
                [0], #274
                [0], #276
                [0], #277
                [0], #278
                [0], #279
                [0], #280
                [0], #281
                [0, 0.05,0.01,0.15], #282
                [0], #283
                [0], #285
                [0], #286
                [0], #287
                [0], #288
                [0], #296
                [0], #297
                [0], #299
               ],dtype=object)

#         bsh = np.array([[0,2,3.5,4.5], #95
#                 [0, 0.25, 0.44, 0.59], #98
#                 [0,0.19,0.34,0.39], #132
#                 [-4.5,0], #142
#                 [-4.5,0], #144
#                 [0, 0.25,0.44,0.59], #151
#                 [0, 0.25,0.44,0.59], #157
#                 [-2.5,0], #168
#                 [-4.5,0], #209
#                 [-4.5,0],#216
#                 [-1.5,0], #218
#                 [0, 0.25, 0.44, 0.59], #226
#                 [0, 0.15], #267
#                 [0], #274
#                 [0], #276
#                 [0], #277
#                 [0], #278
#                 [0], #279
#                 [0], #280
#                 [0], #281
#                 [0,0.15], #282
#                 [0], #283
#                 [0], #285
#                 [0], #286
#                 [0], #287
#                 [0], #288
#                 [0], #296
#                 [0], #297
#                 [0], #299
#                ],dtype=object)
    
    
    
    return bsh, bus


 ################################################################################################################################################################################################


def converte_trafo(tap_pos, tap_neutral, tap_step_percent,valores_taps):
    
    '''
    Converte TAPS conforme equação disponibilizada pelo pandapower.
    
    https://pandapower.readthedocs.io/en/v2.1.0/elements/trafo.html
    
    '''
    
    taps_convertido = tap_neutral + ((valores_taps - 1.0)*(100/tap_step_percent))
    
    
    return taps_convertido


 ################################################################################################################################################################################################


def cria_alcateia(sep,n_lobos):
    
    """"
    
    Cria a alcatéia de lobos cinzentos.
    
    linhas = partículas
    
    colunas = tensões geradores, tap transformadores, susceptâncias shunt, perdas, penalização de tensão, penalização de reativo, penalização de trafo, penalização shunt, fitness
   
    
    """
    
    
    vgen, thetagen, pgen, qgen, p_lim_superior, p_lim_inferior, q_lim_superior, q_lim_inferior,barra = coleta_dados_gen(sep,relatorio=False)
    
    n_vgen=len(vgen)
    
    vbus, theta, v_lim_superior, v_lim_inferior = coleta_dados_vbus(sep,relatorio=False)
    
    tap_pos, tap_neutral, tap_step_percent,valores_taps=coleta_dados_trafo(sep,relatorio=False)
    
    n_taps = len(tap_pos)
    
    bshunt , bus = coleta_dados_bshunt(sep)
    
    n_bshunt = len(bus)
    
    dimensao = n_taps + n_vgen + n_bshunt + 6
    
    alcateia=np.zeros((n_lobos,dimensao),dtype=np.float64)
    
    alcateia[:,0:n_vgen] = np.random.uniform(np.max(v_lim_inferior), np.max(v_lim_superior), size=(n_lobos,n_vgen))
    
    alcateia[:,n_vgen:n_vgen+n_taps]=np.random.choice(valores_taps, size =(n_lobos, n_taps))
    
    i=1
    
    i=1
    
    for bsh in bshunt:
        
        alcateia[:,n_vgen+n_taps+i-1:n_vgen+n_taps+i] = np.random.choice(bsh, size =(n_lobos, 1))
        i=i+1

    return alcateia
    
    
    


 ################################################################################################################################################################################################


def cria_enxame(sep,n_particulas):
    
    """"
    
    Cria o enxame de partículas.
    
    
    linhas = partículas
    
    colunas = tensões geradores, tap transformadores, susceptâncias shunt, perdas, penalização de tensão, penalização de reativo, penalização de trafo, penalização shunt, fitness
    
    """
    
    
    vgen, thetagen, pgen, qgen, p_lim_superior, p_lim_inferior, q_lim_superior, q_lim_inferior,barra = coleta_dados_gen(sep,relatorio=False)
    
    n_vgen=len(vgen)+1
    
    vbus, theta, v_lim_superior, v_lim_inferior = coleta_dados_vbus(sep,relatorio=False)
    
    tap_pos, tap_neutral, tap_step_percent,valores_taps=coleta_dados_trafo(sep,relatorio=False)
    
    n_taps = len(tap_pos)
    
    bshunt , bus = coleta_dados_bshunt(sep)
    
    bshunt = np.array(bshunt)
    
    n_bshunt = len(bus)
    
    dimensao = n_taps + n_vgen + n_bshunt + 6 
    
    enxame=np.zeros((n_particulas,dimensao),dtype=np.float64)
    
    enxame[:,0:n_vgen] = np.random.uniform(np.max(v_lim_inferior), np.max(v_lim_superior), size=(n_particulas,n_vgen))
    
    enxame[:,n_vgen:n_vgen+n_taps]=np.random.choice(valores_taps, size =(n_particulas, n_taps))
    
    i=1
    
    for bsh in bshunt:
        
        enxame[:,n_vgen+n_taps+i-1:n_vgen+n_taps+i] = np.random.choice(bsh, size =(n_particulas, 1))
        i=i+1
        
    return enxame


def cria_enxame_v(sep,n_particulas):
    
    """"
    
    Cria o enxame de partículas.
    
    
    linhas = partículas
    
    colunas = tensões geradores, tap transformadores, susceptâncias shunt, perdas, penalização de tensão, penalização de reativo, penalização de trafo, penalização shunt, fitness
    
    """
    
    
    vgen, thetagen, pgen, qgen, p_lim_superior, p_lim_inferior, q_lim_superior, q_lim_inferior,barra = coleta_dados_gen(sep,relatorio=False)
    
    n_vgen=len(vgen)+1
    
    vbus, theta, v_lim_superior, v_lim_inferior = coleta_dados_vbus(sep,relatorio=False)
    tap_pos, tap_neutral, tap_step_percent,valores_taps=coleta_dados_trafo(sep,relatorio=False)
    
    n_taps = len(tap_pos)
    
    bshunt , bus = coleta_dados_bshunt(sep)
    
    bshunt = np.array(bshunt)
    
    n_bshunt = len(bus)
    
    dimensao = n_taps + n_vgen + n_bshunt + 6
    
    enxame=np.zeros((n_particulas,dimensao),dtype=np.float64)
    
    enxame[:,0:n_vgen] = np.random.uniform(-1*np.ones(n_vgen)*np.max(v_lim_superior), np.ones(n_vgen)*np.max(v_lim_superior), size=(n_particulas,n_vgen))
    
    enxame[:,n_vgen:n_vgen+n_taps]= np.random.uniform(-1*np.ones(n_taps)*np.max(valores_taps), np.ones(n_taps)*np.max(valores_taps), size =(n_particulas, n_taps))
    
    i=1
    
    for bsh in bshunt:
        
        enxame[:,n_vgen+n_taps+i-1:n_vgen+n_taps+i] = np.random.uniform(-1*np.ones(1)*np.max(bsh),np.ones(1)*np.max(bsh),size =(n_particulas, 1))
        i=i+1
        
    return enxame

 ################################################################################################################################################################################################

def pen_bshunt(grupo,n_tap,n_vgen,n_bshunt,sep):
    
    b = grupo[n_tap+n_vgen:n_tap+n_vgen+n_bshunt]
    
    bsh,bus=coleta_dados_bshunt(sep)
    
    penal = 0
    
    i=0

    bs=[]

    for i in range(len(bsh)):
    
        bs.append(np.array(bsh[i]))
    
    for i in range(len(bs)):
                
        if len(bs[i][bs[i]<=b[i]])==0:
            penal=1
            return penal
        if len(bs[i][bs[i]>=b[i]])==0:
            penal=1
            return penal
            
        anterior = bs[i][bs[i]<=b[i]][-1]
        posterior = bs[i][bs[i]>=b[i]][0]
        alfa = np.pi*(np.ceil((anterior/(0.001+posterior-anterior)))-(anterior/(0.001+posterior-anterior)))
        penal = penal + np.square(np.sin((b[i]/(0.001+posterior-anterior))*np.pi+alfa))

    
    return penal    
    



 ################################################################################################################################################################################################


def fluxo_de_pot(grupo, sep):
    
    n_bshunt = len(sep.shunt)
    n_vgen = len(sep.gen)+1
    n_tap = np.abs(sep.trafo['tap_pos']).count()
    
    matrizg, matrizz = matriz_condutancia(sep,relatorio=False)
    
    for linha in range(grupo.shape[0]):
        
        sep.ext_grid['vm_pu']=grupo[linha,0]
        
        sep.gen['vm_pu']=grupo[linha,1:n_vgen]
        
        tap_pos, tap_neutral, tap_step_percent,valores_taps=coleta_dados_trafo(sep,relatorio=False)
        
        sep.trafo['tap_pos'][~pd.isnull(sep.trafo['tap_pos'])]=converte_trafo(tap_pos, tap_neutral, tap_step_percent,grupo[linha,n_vgen:n_vgen+n_tap])
        
        sep.shunt['q_mvar']=grupo[linha,n_vgen+n_tap:n_vgen+n_tap+n_bshunt]*-100
        
        if len(sep.bus)==300:
        
            pp.runpp(sep,algorithm='fdbx',numba=True, init = 'flat', tolerance_mva = 1e-5,max_iteration=10000,trafo_model='pi')
        
        else:
        
            pp.runpp(sep,algorithm='fdbx',numba=True, init = 'dc', tolerance_mva = 1e-6,max_iteration=1000,enforce_q_lims=False,trafo_model='pi')

        
        vbus, theta, v_lim_superior, v_lim_inferior=coleta_dados_vbus(sep,relatorio=False)
        
        grupo[linha,-6] = (sep.res_line['pl_mw'].sum()/100 + sep.res_trafo['pl_mw'].sum()/100) + 0*np.sum((1-sep.res_bus['vm_pu'].values)**2) 
#         grupo[linha,-6]=func_objetivo(vbus,theta,matrizg,matrizz,relatorio=False)

        grupo[linha,-5] = pen_tensao(vbus, v_lim_superior, v_lim_inferior,relatorio=False)
        
        vgen, thetagen, pgen, qgen, p_lim_superior, p_lim_inferior, q_lim_superior, q_lim_inferior,barra = coleta_dados_gen(sep,relatorio=False)
        
        grupo[linha,-4] = pen_ger_reativo(qgen, q_lim_superior, q_lim_inferior,sep,relatorio=False)
        
        grupo[linha,:] = pen_trafo(grupo[linha,:],n_tap,n_vgen)
        
        
        grupo[linha,-2] = pen_bshunt(grupo[linha,:],n_tap,n_vgen,n_bshunt,sep)
  
        
        
    
    return grupo

def fluxo_de_pot_algo(grupo, sep,algo):
    
    n_bshunt = len(sep.shunt)
    n_vgen = len(sep.gen)+1
    n_tap = np.abs(sep.trafo['tap_pos']).count()
    
    matrizg, matrizz = matriz_condutancia(sep,relatorio=False)
    
    for linha in range(grupo.shape[0]):
        
        sep.ext_grid['vm_pu']=grupo[linha,0]
        
        sep.gen['vm_pu']=grupo[linha,1:n_vgen]
        
        tap_pos, tap_neutral, tap_step_percent,valores_taps=coleta_dados_trafo(sep,relatorio=False)
        
        sep.trafo['tap_pos'][~pd.isnull(sep.trafo['tap_pos'])]=converte_trafo(tap_pos, tap_neutral, tap_step_percent,grupo[linha,n_vgen:n_vgen+n_tap])
        
        sep.shunt['q_mvar']=grupo[linha,n_vgen+n_tap:n_vgen+n_tap+n_bshunt]*-100
        
        if len(sep.bus)==300:
        
            pp.runpp(sep,algorithm=algo,numba=True, init = 'result', tolerance_mva = 1e-4,max_iteration=10000)
        
        else:
        
#             pp.runpp(sep,algorithm=algo,numba=True, init = 'flat', tolerance_mva = 1e-6,max_iteration=1000,enforce_q_lims=False,trafo_model='pi')
            pp.rundcpp(sep)
        
        vbus, theta, v_lim_superior, v_lim_inferior=coleta_dados_vbus(sep,relatorio=False)
        
        grupo[linha,-6] = sep.res_line['pl_mw'].sum()/100 + sep.res_trafo['pl_mw'].sum()/100 #func_objetivo(vbus,theta,matrizg,relatorio=False)

        grupo[linha,-5] = pen_tensao(vbus, v_lim_superior, v_lim_inferior,relatorio=False)
        
        vgen, thetagen, pgen, qgen, p_lim_superior, p_lim_inferior, q_lim_superior, q_lim_inferior,barra = coleta_dados_gen(sep,relatorio=False)
        
        grupo[linha,-4] = pen_ger_reativo(qgen, q_lim_superior, q_lim_inferior,sep,relatorio=False)
        
        grupo[linha,:] = pen_trafo(grupo[linha,:],n_tap,n_vgen)
        
        
        grupo[linha,-2] = pen_bshunt(grupo[linha,:],n_tap,n_vgen,n_bshunt,sep)
  
        
        
    
    return grupo


def fluxo_de_pot_q(grupo, sep):
    
    n_bshunt = len(sep.shunt)
    n_vgen = len(sep.gen)+1
    n_tap = np.abs(sep.trafo['tap_pos']).count()
    
    matrizg = matriz_condutancia(sep,relatorio=False)
    
    for linha in range(grupo.shape[0]):
        
        sep.ext_grid['vm_pu']=grupo[linha,0]
        
        sep.gen['vm_pu']=grupo[linha,1:n_vgen]
        
        tap_pos, tap_neutral, tap_step_percent,valores_taps=coleta_dados_trafo(sep,relatorio=False)
        
        sep.trafo['tap_pos'][~pd.isnull(sep.trafo['tap_pos'])]=converte_trafo(tap_pos, tap_neutral, tap_step_percent,grupo[linha,n_vgen:n_vgen+n_tap])
        
        sep.shunt['q_mvar']=grupo[linha,n_vgen+n_tap:n_vgen+n_tap+n_bshunt]*-100
        
        if len(sep.bus)==300:
        
            pp.runpp(sep,algorithm='nr',numba=True, init = 'results', tolerance_mva = 1e-4,max_iteration=1000, enforce_q_lims=True)
        
        else:
        
            pp.runpp(sep,algorithm='nr',numba=True, init = 'results', tolerance_mva = 1e-5,max_iteration=1000,enforce_q_lims=True,trafo_model='pi')
        
        vbus, theta, v_lim_superior, v_lim_inferior=coleta_dados_vbus(sep,relatorio=False)
        
        grupo[linha,-6] = sep.res_line['pl_mw'].sum()/100 + sep.res_trafo['pl_mw'].sum()/100 #func_objetivo(vbus,theta,matrizg,relatorio=False)

        grupo[linha,-5] = pen_tensao(vbus, v_lim_superior, v_lim_inferior,relatorio=False)
        
        vgen, thetagen, pgen, qgen, p_lim_superior, p_lim_inferior, q_lim_superior, q_lim_inferior,barra = coleta_dados_gen(sep,relatorio=False)
        
        grupo[linha,-4] = pen_ger_reativo(qgen, q_lim_superior, q_lim_inferior,sep,relatorio=False)
        
        grupo[linha,:] = pen_trafo(grupo[linha,:],n_tap,n_vgen)
        
        
        grupo[linha,-2] = pen_bshunt(grupo[linha,:],n_tap,n_vgen,n_bshunt,sep)
  
        
    
    return grupo
    



 ################################################################################################################################################################################################


def fitness (grupo,zeta,psi,sigma,omega):
    
# fitness J       perdas         pen tensão         pen q mvar          pen trafo           pen bshunt       
    grupo[:,-1]=(grupo[:,-6])+(zeta*grupo[:,-5])+(psi*grupo[:,-4])+(sigma*grupo[:,-3])+(omega*grupo[:,-2])

    return grupo


 ################################################################################################################################################################################################


def validacao (sep, best_solution,relatorio=True):
       
    valida = fluxo_de_pot(np.array([best_solution]), sep)
    
    if relatorio == True:
        print('Perdas de Potência Ativa [PU]:\n')
        print(valida[0][-6])
        print(' ')

        print('Penalização de Violação de Tensão [PU]:\n')
        print(valida[0][-5])
        print(' ')

        print('Penalização de Violação de Geração de Reativo [PU]:\n')
        print(valida[0][-4])
        print(' ')

        print('Penalização de Violação de TAP Discreto [PU]:\n')
        print(valida[0][-3])
        print(' ')

        print('Penalização de Violação de Bshunt Discreto [PU]:\n')
        print(valida[0][-2])
        print(' ')
        


def validacao_q (sep, best_solution,relatorio=True):
       
    valida = fluxo_de_pot_q(np.array([best_solution]), sep)
    
    if relatorio == True:
        print('Perdas de Potência Ativa [PU]:\n')
        print(valida[0][-6])
        print(' ')

        print('Penalização de Violação de Tensão [PU]:\n')
        print(valida[0][-5])
        print(' ')

        print('Penalização de Violação de Geração de Reativo [PU]:\n')
        print(valida[0][-4])
        print(' ')

        print('Penalização de Violação de TAP Discreto [PU]:\n')
        print(valida[0][-3])
        print(' ')

        print('Penalização de Violação de Bshunt Discreto [PU]:\n')
        print(valida[0][-2])
        print(' ')
        



 ################################################################################################################################################################################################


def otimizacao_gwo_continuo(sep, zeta, psi, sigma, omega, max_iter, n_lobos,valor_inicial,relatorio=True,inicial=True):
        
    alcateia_fit = cria_alcateia(sep,n_lobos)

    if inicial == True:
        
        alcateia_fit[0,:]=valor_inicial
    
    
    j = []
    
    perdas = []
    
    tempo = []
    
    pen_v = []
    
    pen_gq = []
    
    pen_tap = []
    
    pen_bsh = []

    
    v_lim_superior = np.repeat(sep.bus['max_vm_pu'][0], len(sep.gen))
    v_lim_inferior = np.repeat(sep.bus['min_vm_pu'][0], len(sep.gen))
    
    tap_pos, tap_neutral, tap_step_percent,valores_taps = coleta_dados_trafo(sep,relatorio=False)
    
    tap_max = np.repeat(valores_taps[-1], len(tap_pos))
    
    tap_min = np.repeat(valores_taps[0], len(tap_pos))
    
    bsh,b=coleta_dados_bshunt(sep)

    bsh_max=[]
    
    bsh_min=[]
    
    alcateias = []
    
    for bs in bsh:
        bsh_max.append([np.max(bs)])
        bsh_min.append([np.min(bs)])


    maximo = np.expand_dims(np.concatenate((v_lim_superior, tap_max, bsh_max), axis = None), 0)
    minimo = np.expand_dims(np.concatenate((v_lim_inferior, tap_min, bsh_min), axis = None), 0)
     
    
    lim_sup = np.tile(maximo, (n_lobos,1))
    lim_inf = np.tile(minimo, (n_lobos,1))
    
    
    for i in range(0,max_iter):

        start = time.time()
       
        alcateia_fit = fluxo_de_pot(alcateia_fit,sep)
        
        alcateia_fit = fitness(alcateia_fit,zeta,psi,sigma,omega)

        alcateia_fit = alcateia_fit[np.argsort(alcateia_fit[:, -1])]
        
        a = (2 - (i*(2/max_iter)))
        
        mu = 0.5
        sigma = 0.15
        
        r1 = np.random.normal(mu, sigma, size = (n_lobos,alcateia_fit[:,0:-6].shape[1]))
        r2 = np.random.normal(mu, sigma, size = (n_lobos,alcateia_fit[:,0:-6].shape[1]))
        
        #r1 = np.random.random_sample(size = (n_lobos,alcateia_fit[:,0:-6].shape[1]))
        
        #r2 = np.random.random_sample(size = (n_lobos,alcateia_fit[:,0:-6].shape[1]))

        A = (2*a*r1) - a
        
        C = 2*r2
        
        if (i == 0):
        
            lobo_alfa = alcateia_fit[0, :].copy()
            lobo_beta = alcateia_fit[1, :].copy()
            lobo_delta = alcateia_fit[2, :].copy()
            
            alfa = np.expand_dims(alcateia_fit[0,0:-6].copy(), 1)
            beta = np.expand_dims(alcateia_fit[1,0:-6].copy(), 1)
            delta = np.expand_dims(alcateia_fit[2,0:-6].copy(), 1)
            
        
        for t in range(3):

            if (alcateia_fit[t, -1] < lobo_alfa[-1]):

                lobo_alfa = alcateia_fit[0,:].copy()
                    
                alcateias.append(alcateia_fit)

                alfa = np.expand_dims(alcateia_fit[1,0:-6].copy(), 1)

            if (alcateia_fit[t,-1] > lobo_alfa[-1] and alcateia_fit[t,-1] < lobo_beta[-1]):

                lobo_beta = alcateia_fit[1,:].copy()

                beta = np.expand_dims(alcateia_fit[1,0:-6].copy(), 1)

            if (alcateia_fit[t,-1] > lobo_alfa[-1] and alcateia_fit[t,-1] > lobo_beta[-1] and alcateia_fit[t,-1] < lobo_delta[-1]):

                lobo_delta = alcateia_fit[2, :].copy()

                delta = np.expand_dims(alcateia_fit[2,0:-6].copy(), 1)         
        

        d_alfa = np.abs(np.multiply(C, alfa.T) - alcateia_fit[:, 0:-6])*0.1

        d_beta = np.abs(np.multiply(C, beta.T) - alcateia_fit[:, 0:-6])*0.1

        d_delta = np.abs(np.multiply(C, delta.T) - alcateia_fit[:, 0:-6])*0.1

        x_alfa = alfa.T - np.multiply(A, d_alfa)

        x_beta = beta.T - np.multiply(A, d_beta)

        x_delta = delta.T - np.multiply(A, d_delta)

        alcateia_fit[:,0:-6] = (x_alfa + x_beta + x_delta)/3

        alca_estat = alcateia_fit[:,-6:]

        alcateia_fit = np.concatenate(( np.clip(alcateia_fit[:,0:-6], a_min = lim_inf, a_max = lim_sup, out = alcateia_fit[:,0:-6]),alca_estat),axis=1)
        
        
        end = time.time()

        elapsed = end - start
        
        j.append(lobo_alfa[-1])

        perdas.append(lobo_alfa[-6])

        pen_v.append(lobo_alfa[-5])

        pen_gq.append(lobo_alfa[-4])

        pen_tap.append(lobo_alfa[-3])

        pen_bsh.append(lobo_alfa[-2])       
        
        
        tempo.append(elapsed)
        
        if relatorio == True:
            
            print(' ')

            print('Lobo Alfa da Iteração:',i)

            print('Perdas (pu):',lobo_alfa[-6])

            print('Penalização de Tensão:',lobo_alfa[-5])

            print('Penalização de Geração de Reativo:',lobo_alfa[-4])

            print('Penalização do Tap:',lobo_alfa[-3])

            print('Penalização do Bshunt:',lobo_alfa[-2])

            print('Fitness:',lobo_alfa[-1])
            
            print('Tempo: ',elapsed)

            print(' ')

            print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')
            
    if relatorio == True:
        
            plt.figure(figsize=(18,10))
            plt.plot(perdas)
            plt.grid()
            plt.title('Otimização Por Alcateia de Lobos Cinzentos')
            plt.ylabel('Perdas de Potência Ativa (pu)')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(j)
            plt.grid()
            plt.title('Otimização Por Alcateia de Lobos Cinzentos')
            plt.ylabel('Fitness (J)')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(pen_v)
            plt.grid()
            plt.title('Otimização Por Alcateia de Lobos Cinzentos')
            plt.ylabel('Penalização de Tensão')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(pen_gq)
            plt.grid()
            plt.title('Otimização Por Alcateia de Lobos Cinzentos')
            plt.ylabel('Penalização de Geração de Reativo')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(pen_tap)
            plt.grid()
            plt.title('Otimização Por Alcateia de Lobos Cinzentos')
            plt.ylabel('Penalização do TAP')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(pen_bsh)
            plt.grid()
            plt.title('Otimização Por Alcateia de Lobos Cinzentos')
            plt.ylabel('Penalização do BShunt')
            plt.xlabel('Número da Iteração')

    return j,perdas,pen_v,pen_gq,pen_tap,pen_bsh,alcateias,lobo_alfa, lobo_beta, lobo_delta, tempo

    


 ################################################################################################################################################################################################

def otimizacao_pso_continuo(sep, zeta, psi, sigma, omega, max_iter, n_particles,c1,c2,v_amp,valor_inicial,relatorio=True,inicial=True):
        
    enxame_fit = cria_enxame(sep,n_particles)
    
    if inicial == True:
        
        enxame_fit[0,:]=valor_inicial
        
    
    w_max=0.9
    w_min=0.4
    
    j = []
    
    tempo = []
        
    perdas = []
    
    pen_v = []
    
    pen_gq = []
    
    pen_tap = []
    
    pen_bsh = []

    
    v_lim_superior = np.repeat(sep.bus['max_vm_pu'][0], len(sep.gen))
    
    v_lim_inferior = np.repeat(sep.bus['min_vm_pu'][0], len(sep.gen))
    
    tap_pos, tap_neutral, tap_step_percent,valores_taps = coleta_dados_trafo(sep,relatorio=False)
    
    tap_max = np.repeat(valores_taps[-1], len(tap_pos))
    
    tap_min = np.repeat(valores_taps[0], len(tap_pos))
    
    bsh,b=coleta_dados_bshunt(sep)

    bsh_max=[]
    
    bsh_min=[]
    

    for bs in bsh:
        bsh_max.append([np.max(bs)])
        bsh_min.append([np.min(bs)])


    maximo = np.expand_dims(np.concatenate((v_lim_superior, tap_max, bsh_max), axis = None), 0)
    minimo = np.expand_dims(np.concatenate((v_lim_inferior, tap_min, bsh_min), axis = None), 0)
     
    
    lim_sup = np.tile(maximo, (n_particles,1))
    lim_inf = np.tile(minimo, (n_particles,1))
    
    v_anterior = v_amp*cria_enxame(sep,n_particles)

    for i in range(0,max_iter):
        
        start = time.time()
        
        mu, sigma = 0.5, 0.152 # mean and standard deviation

        #r1 = np.random.normal(mu, sigma, size = (n_particles,enxame_fit.shape[1]))
        #r2 = np.random.normal(mu, sigma, size = (n_particles,enxame_fit.shape[1]))
        
        r1 = np.random.random_sample(size = (n_particles,1))
        
        r2 = np.random.random_sample(size = (n_particles,1))
    
       
        enxame_fit = fluxo_de_pot(enxame_fit,sep)
        


        enxame_fit = fitness(enxame_fit,zeta,psi,sigma,omega)

        if i==0:
            
            best_particles = enxame_fit.copy()

            global_best = best_particles[np.argsort(best_particles[:, -1])][0,:].copy()
            

            
           
        for t in range(0,n_particles):
                
            if enxame_fit[t,-1] < best_particles[t,-1]:
                    
                best_particles[t,:] = enxame_fit[t,:].copy()
                    

        global_best = best_particles[np.argsort(best_particles[:, -1])][0,:].copy()
            
        global_matriz = np.tile(global_best, (n_particles,1))   
            
        
        enxame_fit_anterior = enxame_fit.copy()
        
        w_novo = w_max-(w_max-w_min)*(i+1)/max_iter
        v_novo = np.multiply(w_novo,v_anterior.copy()) + c1*np.multiply(r1,(best_particles.copy()-enxame_fit.copy())) + c2*np.multiply(r2,(global_matriz.copy()-enxame_fit.copy()))
        
        enxame_fit_novo = enxame_fit_anterior  + v_novo
        
        v_anterior = v_novo.copy()
        
        enxame_estat = enxame_fit_novo[:,-6:]

        enxame_fit = np.concatenate(( np.clip(enxame_fit_novo[:,0:-6], a_min = lim_inf, a_max = lim_sup, out = enxame_fit_novo[:,0:-6]),enxame_estat),axis=1)   


        
        end = time.time()

        elapsed = end - start

        j.append(global_best[-1])

        perdas.append(global_best[-6])

        pen_v.append(global_best[-5])

        pen_gq.append(global_best[-4])

        pen_tap.append(global_best[-3])

        pen_bsh.append(global_best[-2])
        
        tempo.append(elapsed)
      

        if relatorio == True:
            
            print(' ')

            print('Melhor Global da Iteração:',i)

            print('Perdas (pu):', global_best[-6])

            print('Penalização de Tensão:', global_best[-5])

            print('Penalização de Geração de Reativo:', global_best[-4])

            print('Penalização do Tap:', global_best[-3])

            print('Penalização do Bshunt:', global_best[-2])

            print('Fitness:', global_best[-1])
            
            print('Tempo: ', elapsed)

            print(' ')

            print('_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ ')
            
            
    
    if relatorio == True:
        
            plt.figure(figsize=(18,10))
            plt.plot(perdas)
            plt.grid()
            plt.title('Otimização Por Enxame de Partículas')
            plt.ylabel('Perdas de Potência Ativa (pu)')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(j)
            plt.grid()
            plt.title('Otimização Por Enxame de Partículas')
            plt.ylabel('Fitness (J)')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(pen_v)
            plt.grid()
            plt.title('Otimização Por Enxame de Partículas')
            plt.ylabel('Penalização de Tensão')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(pen_gq)
            plt.grid()
            plt.title('Otimização Por Enxame de Partículas')
            plt.ylabel('Penalização de Geração de Reativo')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(pen_tap)
            plt.grid()
            plt.title('Otimização Por Enxame de Partículas')
            plt.ylabel('Penalização do TAP')
            plt.xlabel('Número da Iteração')
            
            plt.figure(figsize=(18,10))
            plt.plot(pen_bsh)
            plt.grid()
            plt.title('Otimização Por Enxame de Partículas')
            plt.ylabel('Penalização do BShunt')
            plt.xlabel('Número da Iteração')
                       
            
    return j,perdas,pen_v,pen_gq,pen_tap,pen_bsh,global_best, tempo

      



 ################################################################################################################################################################################################


def discreto_bshunt(grupo,n_tap,n_vgen,n_bshunt,sep):
    
    b = grupo[n_tap+n_vgen:n_tap+n_vgen+n_bshunt]
    
    bsh,bus=coleta_dados_bshunt(sep)
    
    penal = 0
    
    discretiza = []
    
    i = 0

    bs = []

    for i in range(len(bsh)):
    
        bs.append(np.array(bsh[i]))
    
    i = 0
    
    for c in bs:
                
        discretiza.append(c[np.argmin(np.abs(c-b[i]))])
        
        i=i+1    
        
    
    return discretiza    
    



 ################################################################################################################################################################################################



def discreto_tap(grupo,n_tap,n_vgen,n_bshunt,sep):
    
    b = grupo[n_vgen:n_vgen+n_tap]
    
    ref = np.arange(start = 0.9, stop = 1.1, step = 0.01)
    
    discretizatap = np.zeros(len(b))
    
    i = 0

    
    for i in range(len(b)):
        
        discretizatap[i]=(ref[np.argmin(np.abs(ref-b[i]))])
                   
    return discretizatap
    
    
