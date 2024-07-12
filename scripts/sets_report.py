
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.config import dops_data
from dopings.doping_set import DopingSet

###############################################################################

# Para cada material
for material in ["graphine", "graphene"]:

    # Para cada elemento
    for dop_elem in ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]:

        # Ler um conjunto de dopagens do armazenamento
        set = DopingSet(dop_elem = dop_elem, dops_info=dops_data[material], mode="look")

        # Reportar estado das otimizações
        set.map_opt(only_report=True, inverse_order=False, verbose=True)