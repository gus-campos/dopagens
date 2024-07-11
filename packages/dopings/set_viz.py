#from dops_set import dops_set
from dopings.config import atoms_data, dirs_data

class set_viz:

    @staticmethod
    def gap_site_graph(set: "dops_set") -> None:
        """
        Gera, para um conjunto completo de dopagens, vizualizações com 
        gráficos de linha, que relacionam o sítio da dopagem com o gap
        da estrutura. É gerada uma visualizaçãp para cada par 
        base-elemento, onde são comparados lado a lado o gap dos 
        diversos sítios de dopagem, em tal base, com tal elemento.

        Parameters
        ----------

        set : dops_set
            Conjunto de dopagens para o qual se deseja gerar a 
            visualizaçãp.
        """

        from pathlib import Path
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        # Tamanho da figura
        plt.rcParams["figure.figsize"] = (15,5)
            
        # Elemento
        dop_elem = set.dop_elem

        # Material
        material = set.dops_info["material"]

        # Para cada base das dopagens
        for base in set.dops_info["bases"]:
        
            # Listas de dados
            sites = []
            homos = []
            lumos = []

            # Para cada estrutura
            for struct in set.structs:

                # Se for da base atual
                if struct.base == base:

                    # Ler dados
                    homo, lumo, gap = struct.homo_lumo()

                    # Separar em diferentes listas
                    sites.append(struct.site)
                    homos.append(homo)
                    lumos.append(lumo)

            # Criando array
            arr = list(zip(sites, homos, lumos))

            # Adicionar a um dataframe
            df = pd.DataFrame(arr, columns=["site", "homo", "lumo"])

            # Plot
            sns.lineplot(data=df, x="site", y="homo", marker="s")
            sns.lineplot(data=df, x="site", y="lumo", marker="s")

            plt.ylim(*set.dops_info["graphs_config"]["ylims"]["gap"])

            # Definindo título
            plt.title(f"HOMO LUMO | {dop_elem}-{base}", fontdict={"fontsize":18})

            # Definindo label do eixo y
            plt.ylabel("HOMO - LUMO")

            # Nome do arquivo
            filename = dirs_data["processing_output"] / material / "gap_graphs" / "site_x_gap" / f"{base}-{dop_elem}.png"

            # Criando diretórios
            Path.mkdir(filename.parent, parents=True, exist_ok=True)

            # Salvando figura
            plt.savefig(filename)

            # Reportar arquivo gerado
            print(filename)

            # Limpando o plot
            plt.clf()

    ###########################################################################
    @staticmethod
    def gap_elem_graph(sets_list: list["dops_set"]) -> None:
        """
        Gera, para um conjunto completo de dopagens, vizualizações com 
        gráficos de linha, que relacionam o elemento da dopagem com o 
        gap da estrutura. É gerada uma visualizaçãp para cada par 
        base-sítio, onde são comparados lado a lado o gap dos 
        diversos elementos de dopagem, em tal base e sítio.

        Parameters
        ----------

        set_list : list[dops_set]
            Conjunto de dopagens para o qual se deseja gerar a 
            visualizaçãp.
        """

        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        from pathlib import Path

        # Material
        material = sets_list[0].dops_info["material"]

        # Lista de todas as estruturas
        all_structs = []

        # Para cada conjunto
        for set in sets_list:
            
            # Para cada estrutura deste conjunto
            for struct in set.structs:

                # Adicionar à lista geral
                all_structs.append(struct)

        # Para cada base
        for base in set.dops_info["bases"]:
            # Para cada sítio
            for site in set.dops_info["sites"][base]:

                # Dados para a plotagem
                dop_elems = []
                homos = []
                lumos = []

                # Para todas as estruturas
                for struct in all_structs:
            
                    # Se tiver dada base e for de dado sítio
                    if struct.base == base and struct.site == site:
                        
                        # Lendo homo lumo
                        homo, lumo, gap = struct.homo_lumo()

                        # Listando dados para plotar
                        dop_elems.append(struct.dop_elem)
                        homos.append(homo)
                        lumos.append(lumo)

                # Criando array
                arr = list(zip(dop_elems, homos, lumos))

                # Adicionar a um dataframe
                df = pd.DataFrame(arr, columns=["elem", "homo", "lumo"])

                # Plot
                sns.lineplot(data=df, x="elem", y="homo", marker="s")
                sns.lineplot(data=df, x="elem", y="lumo", marker="s")

                plt.ylim(*sets_list[0].dops_info["graphs_config"]["ylims"]["gap"])

                # Definindo título
                plt.title(f"HOMO LUMO | {base}-{site}", fontdict={"fontsize":18})

                # Definindo label do eixo y
                plt.ylabel("HOMO - LUMO")

                # Nome do arquivo
                filename = dirs_data["processing_output"] / material / "gap_graphs" / "elem_x_gap" / f"{base}-{site}.png"

                # Criando diretórios
                Path.mkdir(filename.parent, parents=True, exist_ok=True)

                # Salvando figura
                plt.savefig(filename)

                # Reportar arquivo gerado
                print(filename)

                # Limpando o plot
                plt.clf()

    ###########################################################################

    @staticmethod
    def energ_site_graph(sets_list: list["dops_set"], second_var: str=None, 
                         show_cov: bool=True) -> None:
        """
        Gera uma grade de gráficos, onde os gráficos da grade variam em
        elementos no eixo x, e em bases no eixo y, e cada gráfico 
        relaciona o sítio no eixo x, com a energia de formação no eixo 
        y. A linha deste gráfico se dá na cor azul.

        Opcionalmente plota conjuntamente a relação entre alguma métrica
        e o sítio, com plot em linha vermelha. 

        Pode mostrar a covariância dos dois dados naquele gráfico, junto
        ao título.

        Parameters
        ----------
        sets_list : list[dops_set]
            Lista de donjunto de dopagens a serem plotadas.

        stat : {None, "dipole", "std"}, default = None 
            Se a estatística norma do dipolo, desvio padrão das cargas, 
            ou nenhuma, deve ser plotada junto à energia.

        show_cov : bool, default = True
            Se a covariância das duas estatísticas deve ser impressa 
            junto ao título do gráfico. Se stat for None, show_cov é
            ignorada.

        Raises
        ------

        ValueError
            Invalid second_var option.

        ValueError
            At least two sets are needed in sets_list.
        """

        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        from pathlib import Path
        import numpy as np

        def charges_std(struct):
            """
            Calcula o desvio padrão das cargas de uma
            estrutura
            """

            import statistics as stats
            delta_charges = [atom.charge - atoms_data["valency"][atom.elem] for atom in struct.atoms]
            return stats.stdev(delta_charges, xbar=0)

        # Validando opção de second_var
        if second_var not in [None, "dipole", "std"]:
            raise ValueError("Invalid second_var option.")
        
        # Validando sets_list
        if len(sets_list) < 2:
            raise ValueError("At least two sets are needed in sets_list.")

        material = sets_list[0].dops_info["material"]

        if material == "graphene":

            # Tamanho da figura
            plt.rcParams["figure.figsize"] = (320,20)
            plt.rcParams['figure.dpi'] = 100  

        else:

            # Tamanho da figura
            plt.rcParams["figure.figsize"] = (100,40)
            plt.rcParams['figure.dpi'] = 100  
        
        # Firas e eixos
        figure, grid_axis = plt.subplots(len(sets_list[0].dops_info["bases"]), len(sets_list)) 

        # Para conjunto de dopagens
        for j in range(len(sets_list)):
    
            # Conjunto
            set = sets_list[j]
        
            # Elemento
            dop_elem = set.dop_elem

            # Para cada base das dopagens
            for i, base in enumerate(set.dops_info["bases"]):

                # Listas de dados
                sites = []
                energies = []
                charges_balance = []

                # Para cada estrutura
                for struct in set.structs:

                    # Se for da base atual
                    if struct.base == base:

                        # Separar em diferentes listas
                        sites.append(struct.site)
                        energies.append(struct.formation_energy()[1])

                        if second_var == "dipole":
                            charges_balance.append(np.linalg.norm(struct.dipole))
                        
                        elif second_var == "std":

                            charges_balance.append(charges_std(struct))

                if second_var is not None:

                    arr = list(zip(sites, energies, charges_balance))
                    df = pd.DataFrame(arr, columns=["site", "energy", "balance"])

                    # Eixo gêmeo
                    ax2 = grid_axis[i,j].twinx()

                    # Segundo plot
                    sns.lineplot(data=df, x="site", y="balance", marker="s", color="red", ax=ax2)

                    # Definindo label do eixo y
                    grid_axis[i,j].set_ylabel("")
                    ax2.set_ylabel("")

                    # Definindo label do eixo x
                    grid_axis[i,j].set_xlabel("")
                    ax2.set_xlabel("")

                    import numpy as np
                    cov = np.corrcoef(energies, charges_balance)[0,1]

                    # Definindo título
                    if show_cov:
                        grid_axis[i,j].set_title(f"{dop_elem}-{base} - COV = {cov:.2f})", fontdict={"fontsize":30})

                    else:
                        grid_axis[i,j].set_title(f"{dop_elem}-{base}", fontdict={"fontsize":30})

                else:
                    
                    arr = list(zip(sites, energies))
                    df = pd.DataFrame(arr, columns=["site", "energy"])

                    grid_axis[i,j].set_title(f"{dop_elem}-{base}", fontdict={"fontsize":30})

                # Plot
                sns.lineplot(data=df, x="site", y="energy", marker="s", ax=grid_axis[i,j])

                if material == "graphine":

                    # Separador de estruturas
                    for site in ["B1", "C1", "D1"]:
                        grid_axis[i,j].axvline(x=sites.index(site) - 0.5, color="black", linestyle='--')
                        
                elif material == "graphene":
                
                    for site in ["B1", "C1"]:
                        grid_axis[i,j].axvline(x=sites.index(site) - 0.5, color="black", linestyle='--')

        # Nome do arquivo
        filename = f"{second_var if second_var is not None else 'energy'}.png"
        filepath = dirs_data["processing_output"] / material / "energy_correlations" / "sites" / filename

        # Criando diretórios
        Path.mkdir(filepath.parent, parents=True, exist_ok=True)

        # Salvando figura
        plt.savefig(filepath)

        # Reportar arquivo gerado
        print(filepath)

        # Limpando o plot
        plt.clf()
        plt.close()

    ###########################################################################

    @staticmethod
    def energ_elem_graphs(sets_list: list["dops_set"]) -> None:
        """
        Gera, para uma lista de conjuntos de dopagens, visualizações que
        para cada sítio, compara a energia de ligação das dopagens 
        feitas com diferentes elementos.

        Parameters
        ----------

        sets_list : list[dops_set]
            Lista de conjunto de dopagens para os quais se deseja gerar
            as visualizações.
        """

        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        from pathlib import Path

        # Tamanho da figura
        plt.rcParams["figure.figsize"] = (8,4)
        plt.rcParams['figure.dpi'] = 100  

        material = sets_list[0].dops_info["material"]

        for base in sets_list[0].dops_info["bases"]:

            for site in sets_list[0].dops_info["sites"][base]:

                elems = []
                energies = []

                for set in sets_list:
                    for struct in set.structs:

                        if struct.base == base and struct.site == site:

                            elems.append(struct.dop_elem)
                            energies.append(struct.formation_energy()[1] / struct.size)

                # Criando array
                arr = list(zip(elems, energies))

                # Adicionar a um dataframe
                df = pd.DataFrame(arr, columns=["dop_elem", "energy"])

                # Figura e eixo
                fig, ax = plt.subplots() 
                
                # Plot
                ax = sns.lineplot(data=df, x="dop_elem", y="energy", marker="s")

                # Definindo título
                ax.set_title(f"{base}-{site}", fontdict={"fontsize":18})

                # Definindo label do eixo y
                ax.set_ylabel("")

                # Definindo label do eixo x
                ax.set_xlabel("")

                if material == "graphine" or material == "graphene":

                    ax.set_ylim(*sets_list[0].dops_info["graphs_config"]["ylims"]["energy"][base])

                # Nome do arquivo
                filepath = dirs_data["processing_output"] / material / "energy_correlations" / "elems" / f"{base}-{site}.png"

                # Criar diretórios necessário
                Path.mkdir(filepath.parent, parents=True, exist_ok=True)

                # Salvando figura
                fig.savefig(filepath)

                # Reportar arquivo gerado
                print(filepath)

                # Limpando o plot
                plt.clf()
                plt.close()

    ###########################################################################