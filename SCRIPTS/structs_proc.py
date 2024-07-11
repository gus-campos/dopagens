
# Incluindo endereço do diretório mãe
import sys
sys.path.append('../')

# Bibliotecas internas
from DOPINGS_PACKAGE.structure import structure
from DOPINGS_PACKAGE.dops_set import dops_set
from DOPINGS_PACKAGE.struct_viz import struct_viz
from DOPINGS_PACKAGE.global_vars import dops_data, dirs_data

###############################################################################

for material in ["graphine", "graphene"]:

    for dop_elem in ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]:

        # Lê um conjunto de dopagens
        set = dops_set(dop_elem=dop_elem, dops_info=dops_data[material], mode="read")

        # Reportar status
        set.map_opt(only_report=True, verbose=False)
        
        # Gera dados
        set.map(method=structure.output, suffix=".out")
        set.map(method=structure.frame, suffix=".xyz")
        set.map(method=struct_viz.histogram, suffix=".png")
        set.map(method=struct_viz.charges_map, suffix=".png")
        set.map_to_csv(method=structure.homo_lumo)
        set.map_to_csv(method=structure.formation_energy)