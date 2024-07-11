#from dops_set import dops_set
from dopings.config import dops_data, graphine_data, dirs_data
from dopings.graphine_calcs import graphine_calcs as gc

"""
viz
"""

# GRÁFICOS DE GEOMETRIA
class graphine_viz:
    
    def geometry_graph(set: "dops_set") -> None:
        """
        Gera uma visualização composta dos dados de geometria do 
        grafino, que conta, por estrutura, com a torção da superfície em
        três eixos, a dilatação nos braços interno e externo, e o ângulo
        da ligação/torção do átomo no anel, dependendo da estrutura.

        Parameters
        ----------

        set : dops_set
            Conjunto de dopagens de grafino para o qual se deseja gerar
            a visualização.
        """
        import numpy as np
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        def delta_arc(value1: float, value2: float) -> float:
            """
            Traz dois arcos, menores ou iguais a 180, um valor 
            positivo de arco, e calcula a diferença entre eles.
            """

            value1 = (value1 + 360) if (value1 < 0) else value1
            value2 = (value2 + 360) if (value2 < 0) else value2

            return value2 - value1

        # Listas para anexar todos os dados
        bases = []
        sites = []
        # Torções
        torsions0 = []
        torsions1 = []
        torsions2 = []
        # Tamanhos de braços
        int_lengths = []
        ext_lengths = []
        # Bond angle
        bonds_angle = []

        # Para cada estrutura
        for struct in set.structs:

            # Base e site da struct
            bases.append(struct.base)
            sites.append(struct.site)

            # Referência da base
            ref_inter = graphine_data["int_arm_length"][bases[-1]]
            ref_exter = graphine_data["ext_arm_length"][bases[-1]]

            # Erro em relação à base
            int_lengths.append((gc.int_arm_length(struct) - ref_inter) / ref_inter)
            ext_lengths.append((gc.ext_arm_length(struct) - ref_exter) / ref_exter)
            
            # Lista de torções desta estrutura
            torsions = []

            # Calculando torções
            for i in [0, 1, 2]:

                torsions.append(gc.axis_torsion_angle(struct, arm=i))

            # Calculando diferença da base
            for i in [0, 1, 2]:

                torsions[i] = delta_arc(graphine_data["base_torsions"][bases[-1]][i], torsions[i])

            # Passando para listas separadas
            torsions0.append(torsions[0])
            torsions1.append(torsions[1])
            torsions2.append(torsions[2])

            # Se o sítio deste está na lista
            if sites[-1] in ["A1", "C1", "C2", "C3", "C4"]:

                base_atom_torsion = graphine_data["dop_atom_torsion_on_ring"][bases[-1]][sites[-1]]
                bonds_angle.append(gc.dop_atom_torsion(struct) - base_atom_torsion)

            else:

                base_bond_angle = graphine_data["bases_bonds_angles"][bases[-1]][sites[-1]]
                bonds_angle.append(gc.dop_atom_angle(struct) - base_bond_angle)

        #################### PLOTAGEM

        # Listas para anexar todos os dados
        arr_torsions = np.array([bases, sites, torsions0, torsions1, torsions2]).transpose()
        arr_lengths = np.array([bases, sites, int_lengths, ext_lengths]).transpose()
        arr_angles = np.array([bases, sites, bonds_angle]).transpose()

        # Dataframes
        df_torsions = pd.DataFrame(data=arr_torsions, columns=["base", "site", "Eixo 0", "Eixo 1", "Eixo 2"])
        df_lengths = pd.DataFrame(data=arr_lengths, columns=["base", "site", "Braço interno", "Braço externo"])
        df_angles = pd.DataFrame(data=arr_angles, columns=["base", "site", "Ângulos"])


        # Plotando
        for base in df_torsions["base"].unique():

            # Lista de sites
            sites_names = dops_data["graphine"]["sites"][base]

            # Tamanho da figura e dpi
            plt.rcParams["figure.figsize"] = (40,25)
            plt.rcParams['figure.dpi'] = 100  

            # Initialise the subplot function using number of rows and columns 
            figure, axis = plt.subplots(3) 
            # Ajustando subplots
            figure.subplots_adjust(hspace=0.4)

            # Separando dados de certa base e removendo coluna base
            temp_df = df_torsions.query(f"base == '{base}'").copy().drop("base", axis=1)

            # Tranformando em df
            temp_df = temp_df.melt(id_vars='site')

            # Transformando valores em floats
            temp_df["value"] = temp_df["value"].astype(float)

            # Plotando
            sns.barplot(data=temp_df, x='site', y='value', hue='variable', ax=axis[0])

            # Limites eixo y
            axis[0].set_ylim([-25,25])

            #  Linha do 0
            axis[0].axhline(y=0, color='black', linestyle='--')

            # Título
            axis[0].set_title(f"{set.dop_elem}-{base} - Variação da torção", fontdict={'fontsize': 22})
            axis[0].set_ylabel("Variação no ângulo de torção (°)", fontdict={'fontsize': 18})
            # Label do eixo x
            axis[0].set_xlabel("Sítio", fontdict={'fontsize': 20})
            # Fonte dos ticks
            axis[0].tick_params(axis='both', which='major', labelsize=16)

            # Separador de estruturas
            for site_name in sites_names:
                axis[0].axvline(x=sites_names.index(site_name) + 0.5, color="black", linestyle='--')

            # Removendo título da legenda
            handles, labels = axis[0].get_legend_handles_labels()
            axis[0].legend(handles=handles, labels=labels)

            # Fonte da legenda
            plt.setp(axis[0].get_legend().get_texts(), fontsize='18') # for legend text
            plt.setp(axis[0].get_legend().get_title(), fontsize='18') # for legend title 

            #####################################################################################

            # Separando dados de certa base e removendo coluna base
            temp_df = df_lengths.query(f"base == '{base}'").copy().drop("base", axis=1)

            # Tranformando em df
            temp_df = temp_df.melt(id_vars='site')

            # Transformando valores em floats
            temp_df["value"] = temp_df["value"].astype(float)

            # Plotando
            sns.barplot(data=temp_df, x='site', y='value', hue='variable', ax=axis[1])

            # Limites eixo y
            axis[1].set_ylim([-0.25, 0.25])
            # Título
            axis[1].set_title(f"{set.dop_elem}-{base} - Variação do comprimento dos braços", fontdict={'fontsize': 22})
            # Label do eixo y
            axis[1].set_ylabel("Variação no comprimento do braço (%)", fontdict={'fontsize': 18})
            # Label do eixo x
            axis[1].set_xlabel("Sítio", fontdict={'fontsize': 20})
            # Fonte dos ticks
            axis[1].tick_params(axis='both', which='major', labelsize=16)

            #  Linha do 0
            axis[1].axhline(y=0, color='black', linestyle='--')

            # Removendo título da legenda
            handles, labels = axis[1].get_legend_handles_labels()
            axis[1].legend(handles=handles, labels=labels)

            # Separador de estruturas
            for site_name in sites_names:
                axis[1].axvline(x=sites_names.index(site_name) + 0.5, color="black", linestyle='--')

            # Fonte da legenda
            plt.setp(axis[1].get_legend().get_texts(), fontsize='18') # for legend text
            plt.setp(axis[1].get_legend().get_title(), fontsize='18') # for legend title 

            ##################################################################################

            # Separando dados de certa base e removendo coluna base
            temp_df = df_angles.query(f"base == '{base}'").copy().drop("base", axis=1)

            # Tranformando em df
            temp_df = temp_df.melt(id_vars='site')

            # Transformando valores em floats
            temp_df["value"] = temp_df["value"].astype(float)
            
            # Adicionando tipos ao df
            temp_df["types"] = ["Torsão no benzeno" if site[:1] in ["A", "C"] else "Ângulo da ligação" for site in temp_df["site"]]

            # Plotando
            sns.barplot(data=temp_df, x='site', y='value', hue="types", ax=axis[2])          

            # Removendo título da legenda
            handles, labels = axis[2].get_legend_handles_labels()
            axis[2].legend(handles=handles, labels=labels)

            # Fonte da legenda
            plt.setp(axis[2].get_legend().get_texts(), fontsize='18') # for legend text

            # Limites eixo y
            axis[2].set_ylim([-80, 50])
            # Título
            axis[2].set_title(f"{set.dop_elem}-{base} - Variação no ângulo da ligação C-{set.dop_elem}-C", fontdict={'fontsize': 22})
            # Label do eixo y
            axis[2].set_ylabel(f"Variação no ângulo da ligação C-{set.dop_elem}-C (°)", fontdict={'fontsize': 18})
            # Label do eixo x
            axis[2].set_xlabel("Sítio", fontdict={'fontsize': 20})
            # Fonte dos ticks
            axis[2].tick_params(axis='both', which='major', labelsize=16)

            # Linha do 0
            axis[2].axhline(y=0, color='black', linestyle='--')

            # Lista de sites
            sites_names = dops_data["graphine"]["sites"][base]

            # Separador de estruturas
            for site_name in sites_names:    
                axis[2].axvline(x=sites_names.index(site_name) + 0.5, color="black", linestyle='--')

            ##################################################################################

            # Salvando gráfico
            dir_out = dirs_data["processing_output"] / "graphine" /  "geometry_graphs" /  f"{set.dop_elem}-{base}.png"

            from pathlib import Path
            Path.mkdir(dir_out.parent, parents=True, exist_ok=True)

            # Salvar figura
            plt.savefig(dir_out)

            # Reportar
            print(dir_out)

            # Fechando e limpando figure
            plt.close()
            plt.clf()
