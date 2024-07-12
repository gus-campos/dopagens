
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.config import dops_data
from dopings.dops_set import dops_set

###############################################################################

# Para cada material
for material in ["graphine", "graphene"]:

    # Para cada elemento
    for dop_elem in ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]:

        # Escrever arquivos de um conjunto de dopagens
        set = dops_set(dop_elem = dop_elem, dops_info=dops_data[material], mode="write")