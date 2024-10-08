
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

        # Escrever arquivos de um conjunto de dopagens
        set = DopingSet(dop_elem = dop_elem, dops_info=dops_data[material], mode="write")

        # Otimizar em fila as estruturas
        set.map_opt(only_report=False, inverse_order=False, verbose=True)