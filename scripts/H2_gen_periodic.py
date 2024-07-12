
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import structure, atom_data
from dopings.config import dirs_data, h2_gen_data
from dopings.h2_gen import h2_gen

# Externas
import numpy as np
from pathlib import Path

###############################################################################

def gen_H2_periodic(p_struct: "periodic_structure", 
                    nH2: int, output_path: str|Path, both_sides:bool=False, 
                    vertical: bool=False, plot: bool=False):

    """
    Gera um arquivo .xyz da estrutura já adicionada de moléculas de H2.

    Parameters 
    ----------

    p_struct : periodic_structure
        Objeto que representa uma estrutura periódica.

    nH2 : int
        Número máximo de moléculas de H2 a serem adicionadas na 
        estrutura.

    output_path : str|Path
        Caminho onde arquivo .xyz deve ser salvo.

    both_sides : bool
        Se a estrutura deve receber moléculas de H2 em seus dois lados.

    vetical : bool
        Se True, as moléculas serão geradas com eixo perpendicular à 
        estrutura. Se False, serão geradas com eixo paralelo à 
        estrutura, e terão rotação aleatória.

    plot : bool
        Se True, uma visualização interativa da estrutura é exibida após
        ser gerada.
    """

    # Copiando e movendo CM pro (0,0,0)
    struct = h2_gen.to_CM(p_struct.struct.copy())

    # Gerando coordenadas de pontos distribuídos no plano XY
    R2_H2_coords = h2_gen.gen_R2_coords_periodic(p_struct, nH2)

    # Gerando coeficientes que ajustam a curva aos átomos
    struct_coords = h2_gen.gen_arrays(struct)
    C = h2_gen.lin_reg(*struct_coords, n=1)
    
    # Transformando pontos distribuídos no plano, em átomos de H distribuídos 
    # sobre a estrutura
    R3_H2_coords = h2_gen.gen_R3_H2(R2_H2_coords, C, otherside=False, 
                                    vertical=vertical)

    ########################## GERANDO PRO OUTRO LADO #########################

    # Se for para gerar H2 nos dois lados da estrutura
    if both_sides:

        # Gerar mais pontos espalhados no plano
        R2_H2_coords = h2_gen.gen_R2_coords_periodic(p_struct, nH2)
        
        # Gerar pontos espalhados sobre o outro lado da estrura
        R3_H2_coords_otherside = h2_gen.gen_R3_H2(R2_H2_coords, C, 
                                                  otherside=True, 
                                                  vertical=vertical)
        
        # Juntar listas de H dos dois lados
        R3_H2_coords += R3_H2_coords_otherside

    ###########################################################################

    # Adicionar todos os átomos de H na estrutura
    for coord in R3_H2_coords:
        struct.append(atom_data(elem="H", coord=np.round(coord, decimals=8).tolist()))

    # Gerar aquivo .xyz da estrutura
    struct.frame(output_path)
    
    # Se for para plotar, gerar visualização interativa da estrutura com curva
    if plot:
        radius = h2_gen.graphine_radius(struct)
        h2_gen.plot(*struct_coords, radius=radius, C=C)

###############################################################################

# Classe de estruturas periódicas
class periodic_structure:
    """
    Classe de objetos que representam estruturas periódicas.

    Attributes
    ----------

    name : str
        O nome da estrutura

    struct : structure
        Objeto que representa a estrutura, em si, desta estrutura 
        periódica.

    nH2 : int
        Número máximo de H2 a ser gerado para a estrutura.

    vectors : list[np.ndarray[float]]
        Lista com dois vetores no R3 que delimitam a célula unitária.
    """
    def __init__(self, name=None, nH2=None, vectors=None):
        
        "Construtor do objeto periodic structure."
    
        self.name=name
        self.struct=None
        self.nH2=nH2
        self.vectors=vectors

# Listando estruturas periódicas
p_structs = []

# Para cada estrutura
for name in ["g1_s1", "g1_s2", "g1_s3", "g1_s4"]:
    
    # Criar estrutura periódica
    p_struct = periodic_structure(
            
                            name=name, 
                            nH2=h2_gen_data[name]["nH2"],
                            vectors=[np.array(h2_gen_data[name]["vectors"][0]),
                                     np.array(h2_gen_data[name]["vectors"][1])]
                        )
        
    # Acumular na lista de estruturas periódicas
    p_structs.append(p_struct)

# Para cada estrutura e seu índice
for i, p_struct in enumerate(p_structs):

    # Ler estrutura correspondente e adicionar como sua propriedade
    p_structs[i].struct = structure(dir=dirs_data["bases"] / "g1-periodic" / p_struct.name, 
                                    read_from_dir=True)

###############################################################################

# Para cada estrutura periódica
for p_struct in p_structs:

    # De 0, até o número máximo de H2 a serem adicionados
    for nH2 in range(p_struct.nH2 + 1):
        
        # Nome do arquivo .xyz de saída
        output_name = f"{p_struct.name}-{f'{i:03d}'}{'.xyz'}"

        # Caminho de saída
        output_path = dirs_data["h2_gen_output"]/ "periodic" / output_name

        # Gerar estrutura .xyz da estrutura periódica com "i" moléculas de H2
        gen_H2_periodic(p_struct, nH2, output_path, vertical=True, 
                        both_sides=True, plot=False)

###############################################################################

