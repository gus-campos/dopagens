
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import structure
from dopings.struct_viz import struct_viz
from dopings.config import dops_data, dirs_data

# Externas
from pathlib import Path

###############################################################################

def map_to_csv(struct_list, method: "function", material=None, **kwargs) -> None:

    # Lista para receber valores
    values = []

    # Se não tiver estruturas, reportar
    if len(struct_list) == 0:
        raise AssertionError("There is no structs on this set_struct.")

    # Para cada structure
    for struct in struct_list:

        # Executar função e guardar valor
        value = method(struct, **kwargs)
            
        # Dict de identificação
        id = { "Name" : "Bases", "Base" : struct.dir.name, "Site" : "BASE"}

        # Formatação para dados em csv
        if method.__name__ == "homo_lumo":

            value = {

                "homo" : value[0],
                "lumo" : value[1],
                "homo_lumo" : value[2]
            }
        
        # Formatação para dados em csv
        elif method.__name__ == "formation_energy":

            value = {

                "energy" : value[0],
                "per_atom" : value[1],
            }
        
        # Juntando identificação e valores 
        values.append(dict(id, **value))

    # Diretório final
    dir_out = dirs_data["processing_output"] / material / method.__name__ / f"{method.__name__}-BASES.csv"
                                
    # Criando diretório
    Path.mkdir(dir_out.parent, exist_ok=True, parents=True)

    # Escrevendo arquivo CSV
    import csv
    with open(dir_out, 'w') as file:
        # Escreve dicionário como csv, usa como nomes de coluna as chaves de um item
        writer = csv.DictWriter(file, fieldnames=values[0].keys())
        # Escreve cabeçalho, e escreve valores
        writer.writeheader()
        writer.writerows(values)

    # Relatar
    print(dir_out)

###########################################################################

for material in ["graphene", "graphine"]:

    struct_list = []

    for base in dops_data[material]["bases"]:

        param_dict = {
            "graphine" : "3ob",
            "graphene" : "matsci"
        }

        struct = structure(dir=(dirs_data["bases"] / material / base), 
                           read_from_dir=True, 
                           param=param_dict[material], material=material,
                           base=base)
        
        struct_list.append(struct)

        #######################################################################

        for method, suffix in [(structure.output, "out"),
                                (structure.frame, "xyz"),
                                (struct_viz.histogram, "png"),
                                (struct_viz.charges_map, "png")]:
        
        
            method(struct, dirs_data["processing_output"] / material 
                    / method.__name__ / "BASES" / f"{base}.{suffix}")
    
    ###########################################################################

    for method in [structure.homo_lumo, structure.formation_energy]:

        map_to_csv(struct_list, method, material)



