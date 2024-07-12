
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.config import dops_data
from dopings.dops_set import dops_set

###############################################################################

for material in ["graphine", "graphene"]:

    for dop_elem in ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]:

        set = dops_set(dop_elem = dop_elem, dops_info=dops_data[material], mode="write")

        #set.map_opt(only_report=True, inverse_order=False, verbose=False)