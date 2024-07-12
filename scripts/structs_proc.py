
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import Structure
from dopings.doping_set import DopingSet
from dopings.struct_viz import StructViz
from dopings.config import dops_data

###############################################################################

# Para cada material
for material in ["graphine", "graphene"]:

    # Para cada elemento
    for dop_elem in ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]:

        # Ler um conjunto de dopagens
        set = DopingSet(dop_elem=dop_elem, dops_info=dops_data[material], mode="read")
        
        #### Gerar arquivos individuais para todas as estruturas
        set.map(method=Structure.output, suffix=".out")         # Output
        set.map(method=Structure.frame, suffix=".xyz")          # Frame
        set.map(method=StructViz.histogram, suffix=".png")      # Histograma
        set.map(method=StructViz.charges_map, suffix=".png")    # Mapa de cargas

        #### Gerar arquivos CSV com todas as estruturas do conjunto
        set.map_to_csv(method=Structure.homo_lumo)              # CSV de homo-lumo
        set.map_to_csv(method=Structure.formation_energy)       # CSV de energia de formação