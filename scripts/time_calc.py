
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.config import dops_data
from dopings.doping_set import DopingSet

from dopings.struct_read import StructRead

########################################################################################

# Para cada material
for material in ["graphine", "graphene"]:

    soma = 0

    # Lista para acumular conjuntos de dopagens
    sets = []

    # Para todos os elementos
    for dop_elem in ["Al", "B", "Li", "Mg", "N", "Na", "O", "P", "Si", "Ti", "Zn"]:

        # Reportar início de leitura
        print(f"Lendo {material}-{dop_elem}...")

        # Ler conjunto e acumular na lista de conjuntos
        set = DopingSet(dop_elem=dop_elem, dops_info=dops_data[material], 
                        mode="read")

        for struct in set.structs:

            tempo = StructRead.read_time(struct)

            if tempo:
                soma += tempo
            else:
                print(struct.dir.parent.name, struct.dir.name, "-> Tempo não encontrado") 
        

    print("\n", material, "total time:", round(soma/3600/24,1), "dias\n")