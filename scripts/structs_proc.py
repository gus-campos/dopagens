
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import structure
from dopings.dops_set import dops_set
from dopings.struct_viz import struct_viz
from dopings.config import dops_data

###############################################################################

# Para cada material
for material in ["graphine", "graphene"]:

    # Para cada elemento
    for dop_elem in ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]:

        # Ler um conjunto de dopagens
        set = dops_set(dop_elem=dop_elem, dops_info=dops_data[material], mode="read")
        
        #### Gerar arquivos individuais para todas as estruturas
        set.map(method=structure.output, suffix=".out")         # Output
        set.map(method=structure.frame, suffix=".xyz")          # Frame
        set.map(method=struct_viz.histogram, suffix=".png")     # Histograma
        set.map(method=struct_viz.charges_map, suffix=".png")   # Mapa de cargas

        #### Gerar arquivos CSV com todas as estruturas do conjunto
        set.map_to_csv(method=structure.homo_lumo)              # CSV de homo-lumo
        set.map_to_csv(method=structure.formation_energy)       # CSV de energia de formação