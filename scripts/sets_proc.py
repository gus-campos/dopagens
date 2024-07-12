
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.config import dops_data
from dopings.dops_set import dops_set
from dopings.graphine_viz import graphine_viz 
from dopings.set_viz import set_viz 

########################################################################################

# Para cada material
for material in ["graphine", "graphene"]:

    # Lista para acumular conjuntos de dopagens
    sets = []

    # Para todos os elementos
    for dop_elem in ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]:

        # Reportar início de leitura
        print(f"Lendo {material}-{dop_elem}...")

        # Ler conjunto e acumular na lista de conjuntos
        set = dops_set(dop_elem=dop_elem, dops_info=dops_data[material], mode="read")
        sets.append(set)

    # Para cada conjunto
    for set in sets:
        
        # Gerar gráfico de gap por sítio
        set_viz.gap_site_graph(set)

        # Se for grafino, gerar gráfico de geometria de grafino
        if material == "graphine":
            graphine_viz.geometry_graph(set)

    # Gerar gráfico de gap por elemento
    set_viz.gap_elem_graph(sets)

    # Gerar gráfico de energia por elemento
    set_viz.energ_elem_graphs(sets)

    # Gerar gráfico de energia por sítio
    set_viz.energ_site_graph(sets, second_var=None)     # Só energia
    set_viz.energ_site_graph(sets, second_var="dipole") # Correlação com dipolo
    set_viz.energ_site_graph(sets, second_var="std")    # Correlação com std das cargas