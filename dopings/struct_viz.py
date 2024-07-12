from pathlib import Path
import numpy as np

from dopings.config import viz_config, atoms_data, dops_data, dirs_data
from dopings.structure import Structure

class StructViz:
    """
    Classe de métodos estáticos que geram vizualizações para dados de 
    estruturas em forma de arquivos.

    Methods
    -------

    def histogram(struct: Structure, hist_out: str|Path) -> None

        Gera um histograma dos comprimentos de ligação da estrutura, no 
        caminho informado.

    charges_map(struct: Structure, map_out: str|Path=None, 
                picking_mode=False) -> None
        
        Gera uma visualizaçãp de mapa de cargas para dada estrutura, no 
        caminho informado.
    """

    @staticmethod
    def histogram(struct: Structure, hist_out: str|Path, 
                  gen_dat: bool=False) -> None:
        """
        Gera um histograma dos comprimentos de ligação da estrutura, no 
        caminho informado.

        Depende da definição dos limites de comprimentos de ligações, 
        por par de elementos, no arquivo viz_config.json.

        Parameters
        ----------

        struct : Structure
            Estrutura para qual se deseja gerar o histograma.

        hist_out : str|Path
            Endereço de saída do histograma, com extensão.

        gen_dat : bool, default = False
            Se os dados devem ser escritos na forma de um arquivo dat, 
            ao lado das plotagens do histograma.
        """

        def bonds_lengths(struct) -> list[float]:
            """
            Gera uma lista com o comprimento de cada ligação da 
            estrutura.

            Depende da definição do comprimento máximo de uma ligação 
            entre os pares de elementos, em viz_config.json.

            Returns
            -------

            list[float]
                Lista dos comprimentos de cada ligação da estrutura.
            """

            # TRANSFORMAR EM DICIONÁRIO DE LABELS
            bond_crit = viz_config['histogram']['bond_crit']
            
            # Lista para receber comprimentos de ligação
            data = []

            # Salvando quantidade de átomos na estrutura
            tam = struct.atoms.shape[0]

            # Para o índice de cada átomo na estrutura
            for i in range(tam):

                # Para o índice de cada átomo seguinte na estrutura
                for j in range(i+1,tam):

                    # Calcular e salvar comprimento da ligação
                    bond_length = struct.atoms[i].dist_to(struct.atoms[j])

                    # Label do par
                    label = "-".join(sorted([struct.atoms[i].elem, struct.atoms[j].elem]))

                    # Se o comprimento calculado for menor que o critério utilizado
                    if bond_length < bond_crit[struct.material][label]:

                        # Listar valor do comprimento como sendo de uma ligação
                        # válida 
                        data.append(bond_length)

            # Retornar dados
            return data

        def hist_dat(lengths, dat_out) -> None:
            """
            Função que gera um arquivo .dat com a contagem de ligações em cada
            intervalo, no formato do xmgrace.

            Parameters
            ----------

            lengths : numpy.array
                Array com todos os comprimentos de ligação da estrutura.

            dat_out : str|Path
                Endereço do arquivo de arquivo .dat de saída.  
            """

            # Importando bibliotecas necessárias
            import pandas as pd
            from pathlib import Path

            ##################################################################

            # Criando lista de bins 
            bins = minor_ticks

            # Fatiando os dados de comprimentos de ligação em bins, 
            # e armazenando contagens de cada bin
            cortes = pd.cut(x = lengths, 
                            bins = bins, 
                            ordered = True, 
                            include_lowest = True).value_counts()

            # Armazenando os intervalos e as contagens em uma lista
            bins_values = []

            for index in cortes.index:

                bins_values.append([ index.left, index.right, cortes[index] ])

            ####################### Formatando para o arquivo dat

            # Lista pra armazenar strings de intervalos
            dat_data = []

            # Armazenar primeiro limite de intervalo com o valor 0
            dat_data.append(f'{bins_values[0][0]} {0}\n')

            # Para as demais, listar cada limite do intervalo duas vezes: 
            # Uma com o valor anterior, e outra com o próprio valor
            for i in range(len(bins_values)):

                dat_data.append(f'{bins_values[i][0]} {bins_values[i][2]}\n')
                dat_data.append(f'{bins_values[i][1]} {bins_values[i][2]}\n')

            # Transformando endereço em objeto Path
            dat_out = Path(dat_out).expanduser()

            # Criar diretórios necessários
            Path.mkdir(dat_out.parent, exist_ok = True, parents = True)

            # Abrir arquivo em modo de escrita, e escrever linha a linha os dados
            with open(dat_out, mode = 'w') as file:
                
                for linha in dat_data:
                    file.write(linha)

            # Relatar que arquivo dat foi gerado
            print(dat_out)

        # Visualização de dados
        import seaborn as sns
        # Refinamento de visualização
        import matplotlib.pyplot as plt
        # Biblioteca de diretórios
        from pathlib import Path
        # Numpy
        import numpy as np

        ################
        ## Atribuindo
        ################

        # Limites do eixo X
        x_limits = tuple(viz_config['histogram']["x_limits"])  
        # Comprimento individual do bin
        bins_length = viz_config['histogram']["bins_length"]      
        # Distância entre ticks
        minor_tick_delta = viz_config['histogram']["delta_tick_small"]   
        major_tick_delta = viz_config['histogram']["delta_tick_big"]  
        #delta_y_ticks = viz_config['histogram']["delta_y_ticks"]
        ###################################################################

        # Obtendo lista de comprimento e transformando em array
        lengths = np.array(bonds_lengths(struct))
        
        # Guardando quantidade de ligações
        n_ligs = lengths.shape[0]

        if n_ligs == 0:
            raise AssertionError("Zero bonds found.")

        if struct.material is not None:
            y_limits = [0, int(dops_data[struct.material]["n_ligs"][struct.base])]
        else:
            y_limits=[0, n_ligs]

        ####################################################

        ###############################
        #     Gerando histograma      # 
        ###############################
        
        # Criando figura e eixos
        fig = plt.figure()
        ax = plt.axes()

        # Plotando gráfico
        sns.histplot(lengths, binwidth = bins_length, binrange = x_limits, ax = ax, stat = 'count')

        # Transformando endereço de saída em objeto Path
        hist_out = Path(hist_out).expanduser()

        # Definindo nome dos eixos
        plt.xlabel('Comprimento de ligação', fontdict={'fontsize': 22})
        plt.ylabel('Contagem', fontdict={'fontsize': 22})

        # Limites de valor dos eixos x e y no plot
        plt.xlim(*x_limits)
        plt.ylim(*y_limits)

        ############## TICKS X

        # Valores dos ticks pequenos, calculando apenas uma vez por execução da main
        if 'minor_ticks' not in globals():

            Delta = x_limits[1] - x_limits[0]
            delta = minor_tick_delta

            global minor_ticks
            minor_ticks = [ x_limits[0] + (delta)*i for i in range(int(Delta/delta)+1) ]

        # Valores dos ticks grandes, calculando apenas uma vez por execução da main
        if 'major_ticks' not in globals():

            Delta = x_limits[1] - x_limits[0]
            delta = major_tick_delta

            global major_ticks
            major_ticks = [ x_limits[0] + (delta)*i for i in range(int(Delta/delta)+1) ] 

        # Ticks do eixo y: calculando apenas uma vez para cada tamanho de molécula
        global last_size

        # Se a lista ainda não foi declada, ou se o último tamanho tiver mudado
        if 'y_ticks' not in globals() or (last_size != n_ligs):

            ############ Substituir último com valor máximo
            
            # Atualizar último tamanho
            last_size = n_ligs

            # Calcular deltas
            Delta = y_limits[1] - y_limits[0]

            # Valor tal que haja 40 intervalos no eixo y
            delta_y_ticks = Delta // 20

            # Chamando de "delta"
            delta = delta_y_ticks

            # Calcular ticks
            global y_ticks
            y_ticks = [ y_limits[0] + (delta)*i for i in range(int(Delta/delta)+1) ]
            
            # Substituir último tick pelo valor máximo
            y_ticks[-1] = y_limits[1]

            # Atualizar último tamanho
            last_size = n_ligs

        ### Atribuindo lista de ticks aos eixos
        # eixo x, ticks menores
        ax.set_xticks(minor_ticks, minor = True)
        # eixo x, ticks maiores
        ax.set_xticks(major_ticks, minor = False)
        # eixo y
        ax.set_yticks(y_ticks, minor = False)

        ### Tamanho das labels, e dos traços dos ticks
        # eixo x, ticks menores
        ax.tick_params(axis='y', which='major', labelsize=12, length=8, width=1)
        # eixo x, ticks maiores
        ax.tick_params(axis='x', which='major', labelsize=13, length=8, width=1)
        # eixo y
        ax.tick_params(axis='x', which='minor', labelcolor='w', length=4)

        # Inclinação da label dos ticks: 45 graus, sentido horário
        plt.xticks(rotation=-45)

        ###################### FIGURA

        # Salvando tamanho que deve ter a figura
        fig_size = viz_config['histogram']['fig_size']

        # Definindo tamanho da figura
        fig.set_size_inches(fig_size[0], fig_size[1])

        # Convertendo endereço para Path
        hist_out = Path(hist_out).expanduser()

        # Criando diretório para receber resultados (se ele já não existir)
        Path.mkdir(hist_out.parent, parents = True, exist_ok = True)

        #####################################################


        # Criando e definindo título
        plt.title(hist_out.stem + ' - nº de ligações: ' + str(n_ligs), fontdict={'fontsize': 20})

        ######### Quando há informações de dopagem #########
            
        # Calculando quantidade de ligações esperadas
        if struct.material is not None:
                
            expected = int(dops_data[struct.material]["n_ligs"][struct.base])

            # Se for diferente do esperado, avisar
            if n_ligs != expected:

                # Aviso
                warning = f'[Quantidade de ligações encontradas] - [Quantidade de ligações esperadas] = {n_ligs - expected:+}'
                # Imprimir na tela
                print(warning)
                # Imprimir no gráfico
                plt.suptitle(warning, fontsize = 13, color="red")

        # Salvando figura e fechando
        plt.savefig(hist_out)
        plt.clf()
        plt.close()

        # Imprimindo endereço do arquivo gerado
        print(hist_out)

        # Gerar dat no mesmo endereço, porém com extensão diferente
        if gen_dat:
            pass
            hist_dat(lengths, hist_out.with_suffix(".dat"))

    @staticmethod
    def charges_map(struct: Structure, map_out: str|Path=None, 
                    picking_mode=False) -> None:
        
        """
        Gera uma visualizaçãp de mapa de cargas para dada estrutura, no 
        caminho informado.

        Necessita que o número de elétrons na camada de valência do 
        elemento esteja definido em "atoms_data".
        
        Parameters
        ----------

        struct : Structure
            Estrutura para qual se deseja gerar o mapa de cargas.

        map_out : str
            Endereço de saída do mapa de cargas, com extensão inclusa.

        picking_mode : bool, default = False
            Se True, invés de gerar um arquivo com a plotagem, inicia
            um plot interativo, onde é possível clicar em um átomo
            para que sua carga seja impressa no terminal com 3 casas de
            precisão.
        """

        import pandas as pd

        # Definindo escala usada
        def scale(value: float, rev: bool=False) -> float:
            """
            Converte valores de "x" para valores da escala empregada.
            
            Esta escala é usada para que 

            Parameters
            ----------

            rev : bool, default = False
                Faz o cálculo inverso, tansformando escala em valor.

            Returns
            -------

            float
                Valor na escala.
            """

            from math import log, exp

            # Limite da escala
            x_max = 0.5
            # Curvatura
            c = 3.5    # Antigo -> 5.0
            # Transladador - Função passa por (0,0)
            d = exp(-c)
            # Escalador - Função passa por (x_max, x_max)
            k = x_max / ( (log(x_max + d)) + c )

            # Se for reverso
            if rev:

                # Chamando valor de y
                y = value

                # Se positivo
                if y >= 0:

                    return exp((y/k)-c) - d
                
                # Se negativo
                elif y < 0:

                    return -exp((y/-k)-c) + d

            # Se não for reverso
            else:
                
                # Chamando valor de x
                x = value

                # Se positivo
                if x >= 0:
                    return k * (log(x+d) + c)
                
                # Se negativo
                elif x < 0:
                    return -k * (log(-x+d) + c)

        # Carga dados da estrutura como data frame
        def data_to_df(struct: Structure) -> pd.DataFrame:
            """
            Cria um pd.DataFrame com os dados da estrura. O DataFrame é 
            melhor compreendido pelo seaborn, biblioteca que gera as 
            vizualizações.

            Parameters
            ----------

            struct : Structure
                Estrutura da qual se quer extrair os dados.

            Returns
            -------

            pd.DataFrame
                Dados da estrutura em formato DataFrame.
            """

            # Lista para acumular dados de cada átomo
            data = []

            # Para cada átomo na estrutura...
            for atom in struct.atoms:

                # Ordenar dados em uma lista, e anexá-la na lista 2D "data", 
                # onde cada linha é um átomo e cada coluna um dado desse átomo
                try:
                    data.append([atom.elem, atom.coord[0], atom.coord[1], 
                                 atom.charge, 
                                 atoms_data["valency"][atom.elem]])
                except:
                    raise KeyError(f"valency not found for element {atom.elem}")

            # Transformar lista "data" em um DataFrame com as colunas 
            # devidamente nomeadas
            data = pd.DataFrame(data)

            # Mudando colunas para números de 0 a 4
            data.columns = range(data.columns.size)

            # Renomeando colunas adicionadas, e armazenando em "df"
            return data.rename(columns={0:"Elem", 1:"X", 2:"Y", 3:"Charge", 
                                        4:"Valency"})
        
        # Dados
        import pandas as pd
        # Visualização de dados
        import seaborn as sns
        # Refinamento de visualização
        import matplotlib.pyplot as plt
        # Biblioteca de diretórios
        from pathlib import Path

        ###################### Lendo e tratando dados ######################

        # Carregando dados
        df = data_to_df(struct)

        # Se qualquer carga for None, reportar
        if df["Elem"].isnull().sum() != 0:

            raise TypeError("Not all atoms have defined element values.")
        
        # Se qualquer coordenada for None, reportar
        if df["X"].isnull().sum() != 0 or df["Y"].isnull().sum() != 0:

            raise TypeError("Not all atoms have a defined coordinates value.")
        
        # Se qualquer carga for None, reportar
        if df["Charge"].isnull().sum() != 0:

            raise TypeError("Not all atoms have defined charges values.")

        ################## TRANSFORMANDO DADOS DE CARGA #########################

        ###### Transformações

        #### Mudando de "População de Elétrons" para "Cargas Elementares Ganhas".
        
        # Subtraindo número de elétrons de valência da população de elétrons.
        df['Charge'] = df['Charge'] - df['Valency']
        
        # Verificando soma
        if df["Charge"].sum() > 0.001 or df["Charge"].sum() < -0.001:
            raise AssertionError("Charges total sum is not zero.")
                
        # Arredondando em 3 casas decimais
        df['Charge'] = df['Charge'].round(3)

        # Excluindo dado não mais necessário
        del df['Valency']

        ############### Criando colunas de escala ################

        # Criando coluna de valores de carga escalonados
        df['scale'] = df['Charge'].map(scale).round(3)

        ####################################################

        # Verificando se nos dados há valores fora dos limites da escala 
        # Salvando valores limites
        values_limits = tuple(viz_config["charges_map"]["values_limits"])

        # Lista para receber cargas fora dos limites
        out_of_bounds = []

        # Para cada carga nos dados...
        for charge in df['Charge']:

            # Se abaixo do valor inferior...
            if (charge < values_limits[0] or charge > values_limits[1]):

                # Anexar à lista de cargas fora dos limites
                out_of_bounds.append(charge)

        ################### Configuraçõs do plot ##############################

        ## Estimando um bom tamanho médio dos pontos
        #tam =  800/(sqrt(struct.size) - 5) + 50
        ## Amplitude de tamanhos
        #tam_raio = 100
        ## Tamanhos máximos e mínimos
        #tam_norm = (tam-tam_raio, tam+tam_raio)
                
        tam_norm = (50, 300)

        # Paleta de cores usada --- "Spectral"
        # mais negativo / mais elétron -> mais vermelo (e maior)
        # menos positivo / menos elétrons -> mais azul (e menor)
        palette = sns.color_palette("Spectral_r", as_cmap=True)

        ######################### PLOTANDO ###########################

        # Criando eixo e figura
        figura, ax = plt.subplots(1)

        # Plotando átomos de C
        data = df.query(" Elem == 'C' ")

        sns.scatterplot(data = data,
                        x = 'X', y = 'Y', 
                        hue = 'scale', palette = palette, hue_norm = values_limits, 
                        size = 'scale', sizes = tam_norm, size_norm = values_limits, 
                        legend = True, 
                        ax=ax, picker=True)
        
        # Plotando átomos de H
        # Nesse caso os cículos terão a metade do tamanho, e contarão com um 
        # contorno para diferenciá-los dos demais
        data = df.query(" Elem == 'H' ")

        sns.scatterplot(data = data,
                        x = 'X', y = 'Y',
                        hue = 'scale', palette = palette, hue_norm = values_limits, 
                        size = 'scale', sizes = tam_norm, size_norm = values_limits, 
                        legend = False,
                        ax=ax, picker=True,
                        linestyle = '-', linewidth = 0.8, edgecolor = 'black')
        

        # Listando elementos dopantes
        dop_elems = df.query( "Elem != 'C' and Elem != 'H' " )["Elem"].unique()

        # Se tiver átomos dopantes
        if len(dop_elems) > 0:

            # Para cada elemento dopante
            for elem in dop_elems:

                # Separar dados
                data = df.query(f" Elem == '{elem}' ")

                # Plotando átomos dos dopantes 
                sns.scatterplot(data = data, 
                                x = 'X', y = 'Y', 
                                hue = 'scale', palette = palette, hue_norm = values_limits, 
                                size = 'scale', sizes = tam_norm, size_norm = values_limits, 
                                legend = False,
                                ax=ax, picker=True,
                                marker = "h", linestyle = '-', linewidth = 1, edgecolor = 'black')
                
        # Abrir scatterplot interativo
        if picking_mode:

            # Seletor de cargas
            def onpick(event):

                """
                Imprime elemento e carga do átomo selecionado.
                """

                # Pegando id
                index = event.ind[0]

                # Coordenadas
                x = event.artist.get_offsets()[index][0]
                y = event.artist.get_offsets()[index][1]

                # Imprimindo linha do dataframe que tem coordenadas coincidentes
                data = df[(df["X"] == x) & (df["Y"] == y)].values[0]

                # Imprimir cargas
                print(f"{data[0]} {data[3]:.3f}")

            ax.figure.canvas.mpl_connect("pick_event", onpick)
            plt.show()

            return

        ####################  LEGENDA DE VALORES FIXOS #########################

        import matplotlib.lines as mlines

        # Valores equidistantes da escala de cores
        legend_scale_values = sorted([-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5])

        # Transformando lista de valores de escala, para lista de valores de carga
        legend_charge_values = sorted([ round(scale(valor, rev=True), 2) for valor in legend_scale_values ])

        #########################################################
        #   OCULTANDO LINHAS DESNECESSÁRIAS 
        #
        #   Mesmo se houver valor fora da norma, 
        #   será exibido apenas até o valor máximo definido
        #########################################################

        # Achando valores extremos de cargas
        max_charge = df["Charge"].max()
        min_charge = df["Charge"].min()

        #######################################################################
        # Para a legenda não queremos exibir valores desnecessários, mas é
        # desejado que haja pelo 1 valor maior que a maior carga, e um valor
        # menor que a menor carga.
        #######################################################################

        # Valores da legenda
        values = legend_charge_values

        ############ Valor maior ##################

        # Inicializando maior valor
        max_index = None

        # Encontrando primeiro valor maior que o maior
        for i in range(len(values)):

            if values[i] >= max_charge:

                max_index = i
                break

        # Se não foi encontrado valor maior, usar o valor máximo
        if max_index is None:

            max_index = i

        ###########################################

        # Inicializando menor valor
        min_index = None

        # Encontrando primeiro valor maior que o maior
        for i in reversed(range(len(values))):

            if values[i] <= min_charge:

                min_index = i
                break

        # Se não foi encontrado valor maior, usar o valor mínimo
        if min_index is None:

            min_index = i

        ############## Criando legenda final de valores exibidos na legenda #####################
            
        temp = []

        # Para cada index na lista de valores
        for index in list(range(len(values))):

            # Se estiverem entre os index máximo e o mínimo (ou forem os próprios)
            if min_index <= index and index <= max_index:

                # Anexar o valor correspondente à lista temporária
                temp.append(values[index])

        # Sobrescrevendo lista total dos valores de legenda
        legend_charge_values = temp

        ###########################################

        #################################
        # Criando legenda
        #################################

        # Lista para receber cada linha da legenda
        handles = []
        
        values_amplitude = max(legend_scale_values) - min(legend_scale_values)

        # Normalizador dos valores
        def values_normalizer(value: float, to_size: bool=False) -> float:    
            """
            Normaliza os valores para a escala usada.
            """

            normalized_value = (value - values_limits[0]) / values_amplitude

            if to_size:
                # Tamanho pode ir de 3 até 13 (3+10), já que o valor está normalizado
                return 3 + 10*normalized_value

            else:
                return normalized_value

        # Pra cada valor, uma linha de legenda. Varrer na ordem inversa da lista, para valores positivos ficarem em cima
        for value in legend_charge_values[::-1]:

            # Modificando valor para a escala utilizada
            color_value = scale(value)

            # Criando linha
            # Encontrando cores, definindo marcador, 
            # convertendo tamanho, definindo label
            line = mlines.Line2D([], [], color=palette(values_normalizer(color_value)),
                                marker='o', linestyle='',
                                markersize = values_normalizer(color_value, to_size=True),
                                label = value)

            handles.append(line)

        # Definindo título da legenda, tamanhos das fontes e criando legenda
        plt.legend(title = "Elets. ganhos", loc = "upper left",
                        title_fontsize = 12, fontsize = 10,
                        shadow = True,
                        handles = handles)

        ######################## Título e eixos #######################

        # Se tiver sido passado endereço de saída
        if map_out is not None:
        
            # Transformando caminho em objetos Path
            map_out = Path(map_out).expanduser()

            # Criando e definindo título
            plt.title(map_out.stem + ' - nº de átomos: ' + str(struct.size), fontdict={'fontsize': 20})

        # Subtítulo avisando sobre a presença de cargas fora dos limites
        if len(out_of_bounds) > 0:
            # Imprimindo no subtítulo do gráfico
            plt.suptitle(f"Cargas fora dos limites da escala: {out_of_bounds} - Limites: {list(values_limits)}", fontsize = 13, color="red")

        # Desabilitando ticks
        plt.xticks([]) 
        plt.yticks([]) 

        # Removendo legenda dos eixos
        plt.xlabel('')
        plt.ylabel('')

        ######################## SALVANDO ##########################

        # Tamanho da figura
        fig_size = viz_config["charges_map"]["fig_size"]

        # Difinindo dimensões da figura
        figura.set_size_inches(fig_size[0], fig_size[1])

        # Se tiver sido passado endereço de saída
        if map_out is not None:

            # Criando diretório para receber resultados (se ele já não existir)
            Path.mkdir(map_out.parent, parents = True, exist_ok = True)

            #Salvando como                           
            plt.savefig(map_out)
            plt.clf()
            plt.close()

        # Relatar gráfico gerado
        print(map_out)

        # Subtítulo avisando sobre a presença de cargas fora dos limites
        if len(out_of_bounds) > 0:
            # Reportando pelo terminal
            print(f"Cargas fora dos limites da escala: {out_of_bounds} - Limites: {list(values_limits)}\n")