
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.config import dops_data
from dopings.dops_set import dops_set
from dopings.graphine_viz import graphine_viz 
from dopings.set_viz import set_viz 

########################################################################################

for material in ["graphine", "graphene"]:

    sets = []

    # Todos os elementos
    todos = ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]

    for dop_elem in todos:

        print(f"Lendo {material}-{dop_elem}...")

        sets.append(dops_set(dop_elem=dop_elem, 
                             dops_info=dops_data[material], 
                             mode="read"))

    for set in sets:

        if material == "graphine":
            graphine_viz.geometry_graph(set)
        
        # Gap dos sítios
        set_viz.gap_site_graph(set)

    # Gap dos elementos
    set_viz.gap_elem_graph(sets)

    # Energia por sítio
    set_viz.energ_site_graph(sets, second_var=None)
    set_viz.energ_site_graph(sets, second_var="dipole")
    set_viz.energ_site_graph(sets, second_var="std")

    # Energia por elemento
    set_viz.energ_elem_graphs(sets)