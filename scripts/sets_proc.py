
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.config import dops_data
from dopings.doping_set import DopingSet
from dopings.graphine_viz import GraphineViz 
from dopings.set_viz import SetViz 

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
        set = DopingSet(dop_elem=dop_elem, dops_info=dops_data[material], mode="read")
        sets.append(set)

    # Para cada conjunto
    for set in sets:
        
        # Gerar gráfico de gap por sítio
        SetViz.gap_site_graph(set)

        # Se for grafino, gerar gráfico de geometria de grafino
        if material == "graphine":
            GraphineViz.geometry_graph(set)

    # Gerar gráfico de gap por elemento
    SetViz.gap_elem_graph(sets)

    # Gerar gráfico de energia por elemento
    SetViz.energ_elem_graphs(sets)

    # Gerar gráfico de energia por sítio
    SetViz.energ_site_graph(sets, second_var=None)     # Só energia
    SetViz.energ_site_graph(sets, second_var="dipole") # Correlação com dipolo
    SetViz.energ_site_graph(sets, second_var="std")    # Correlação com std das cargas