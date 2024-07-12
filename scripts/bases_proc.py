
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import Structure
from dopings.struct_viz import StructViz
from dopings.config import dops_data, dirs_data

# Externas
from pathlib import Path

###############################################################################

def map_bases_to_csv(struct_list: list[Structure], method: "function", 
               material=None) -> None:

    """
    Aplica um função à uma lista de estruturas, recolhe seus retornos, e
    gera um cvs com os dados.

    Parameters
    ----------

    struct_list : list[Structure]
        Lista das estruturas para quais se deseja gerar o arquivo.

    method : function
        Função que se deseja aplicar a todas as estruturas da lista.

    material : str
        Material da estrutura.
    """

    # Lista para acumular os dados de retorno
    values = []

    # Se não tiver estruturas na lista, reportar
    if len(struct_list) == 0:
        raise AssertionError("There is no structs on this set_struct.")

    # Para cada estrutura na lista
    for struct in struct_list:

        # Executar função e guardar valor
        value = method(struct)
            
        # Criar dicionário de identificação
        id = { "Name" : "Bases", "Base" : struct.dir.name, "Site" : "BASE"}

        ###### Estruturar dados de retorno de acrodo com a função #############
        if method.__name__ == "homo_lumo":

            value = {
                "homo" : value[0],
                "lumo" : value[1],
                "homo_lumo" : value[2]
            }
        
        elif method.__name__ == "formation_energy":

            value = {
                "energy" : value[0],
                "per_atom" : value[1],
            }
        
        else:
            raise AssertionError("Função inválida.")
        #######################################################################
        
        # Juntando dados identificação com valores retornados 
        values.append(dict(id, **value))

    # Diretório final
    dir_out = dirs_data["processing_output"] / material / method.__name__ / f"{method.__name__}-BASES.csv"
                                
    # Criando diretórios necessários
    Path.mkdir(dir_out.parent, exist_ok=True, parents=True)

    # Escrevendo arquivo CSV
    import csv
    with open(dir_out, 'w') as file:
        writer = csv.DictWriter(file, fieldnames=values[0].keys())
        writer.writeheader()
        writer.writerows(values)

    # Reportar aquivo gerado
    print(dir_out)

###############################################################################

# Para cada material
for material in ["graphene", "graphine"]:

    # Criar lista para acumular estruturas
    struct_list = []

    # Para cada base de tal material
    for base in dops_data[material]["bases"]:

        param_dict = { 
            "graphine" : "3ob",
            "graphene" : "matsci"
        }

        # Ler estrutura
        struct = Structure(dir=dirs_data["bases"] / material / base, 
                           read_from_dir=True, param=param_dict[material], 
                           material=material, base=base)

        #######################################################################

        # Para cada função e sufixo da lista, de arquivo individual
        for method, suffix in [(Structure.output, "out"),
                               (Structure.frame, "xyz"),
                               (StructViz.histogram, "png"),
                               (StructViz.charges_map, "png")]:
        
            # Caminho do arquivo de saída
            output_path = (dirs_data["processing_output"] / material 
                           / method.__name__ / "BASES" / f"{base}.{suffix}")

            # Gerar aquivo
            method(struct, output_path)
        
        #######################################################################

        # Acumular estrutura em lista, para geração dos arquivos csv
        struct_list.append(struct)
    
    ###########################################################################

    # Para cada função de arquivo csv da lista
    for method in [Structure.homo_lumo, Structure.formation_energy]:

        # Mapear função e gerar o cvs
        map_bases_to_csv(struct_list, method, material)



