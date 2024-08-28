
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import Structure
from dopings.atom import Atom
from dopings.config import dirs_data, h2_gen_data
from dopings.h2_gen import H2Gen

# Externas
import numpy as np
from pathlib import Path

###############################################################################

def gen_H2_periodic(p_struct: "PeriodicStructure", 
                    H2_count: int, output_path: str|Path, both_sides:bool=False, 
                    vertical: bool=False, random_rot=False, plot: bool=False):

    """
    Gera um arquivo .xyz da estrutura já adicionada de moléculas de H2.

    Parameters 
    ----------

    p_struct : PeriodicStructure
        Objeto que representa uma estrutura periódica.

    H2_count : int
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

    # Copiando struct
    struct = p_struct.struct.copy()

    # Movendo célula pro (0,0,0)
    for i in range(struct.size):
        struct.atoms[i].coord -= p_struct.initial_vec

    # Se both sides, distribuir contagem dos dois lados
    H2_gen_count =  (H2_count//2) if (both_sides) else (H2_count)

    # Gerando coordenadas de pontos distribuídos no plano XY
    R2_H2_coords = H2Gen.gen_R2_coords_periodic(p_struct, H2_gen_count)

    # Gerando coeficientes que ajustam a curva aos átomos
    struct_coords = H2Gen.gen_arrays(struct)
    C = H2Gen.lin_reg(*struct_coords, n=1)
    
    # Transformando pontos distribuídos no plano, em átomos de H distribuídos 
    # sobre a estrutura
    R3_H2_coords = H2Gen.gen_R3_H2(R2_H2_coords, C, otherside=False, 
                                   vertical=vertical, random_rot=random_rot)

    ########################## GERANDO PRO OUTRO LADO #########################


    # Se for para gerar H2 nos dois lados da estrutura
    if both_sides:

        # Se both sides, distribuir contagem dos dois lados
        H2_gen_count =  (H2_count//2 + H2_count%2) if (both_sides) else (H2_count)
        
        # Gerar mais pontos espalhados no plano
        R2_H2_coords = H2Gen.gen_R2_coords_periodic(p_struct, H2_count//2 + H2_count%2)
        
        # Gerar pontos espalhados sobre o outro lado da estrura
        R3_H2_coords_otherside = H2Gen.gen_R3_H2(R2_H2_coords, C, 
                                                  otherside=True, 
                                                  vertical=vertical,
                                                  random_rot=random_rot)
        
        # Juntar listas de H dos dois lados
        R3_H2_coords += R3_H2_coords_otherside

    ###########################################################################

    # Adicionar todos os átomos de H na estrutura
    for coord in R3_H2_coords:
        struct.append(Atom(elem="H", coord=np.round(coord, decimals=8).tolist()))

    # Gerar aquivo .xyz da estrutura
    struct.frame(output_path, ignore_charges=True)
    
    # Se for para plotar, gerar visualização interativa da estrutura com curva
    if plot:
        H2Gen.plot(struct, C=C)

###############################################################################

# Classe de estruturas periódicas
class PeriodicStructure:
    """
    Classe de objetos que representam estruturas periódicas.

    Attributes
    ----------

    name : str
        O nome da estrutura

    struct : Structure
        Objeto que representa a estrutura, em si, desta estrutura 
        periódica.

    vectors : list[np.ndarray[float]]
        Lista com dois vetores no R3 que delimitam a célula unitária.
    """
    def __init__(self, name=None, struct=None, vectors=None, 
                initial_vec=None):
        
        "Construtor do objeto periodic Structure."
    
        self.name=name
        self.struct=struct
        self.vectors=vectors
        self.initial_vec = initial_vec

root_dir = Path("../input/bases/periodic")

p_structs = []

for base_name in ["g1", "g2", "g3", "g4", "g5"]:
    for cell_size in ["s2"]:

        base = f"{base_name}_{cell_size}"
        dir = root_dir / base
        struct = Structure(dir=dir, read_from_dir=True)

        # Encontrando vetores e criando p_struct
        with open(dir / "geo_end.gen") as file:

            # Linhas
            lines = file.read().splitlines()

            # Ponto inicial e vetores da célula
            initial_vec = np.array([float(coord) for coord in lines[-4].split()])
            vec0 = np.array([float(coord) for coord in lines[-3].split()])
            vec1 = np.array([float(coord) for coord in lines[-2].split()])

        # p_struct
        p_structs.append(PeriodicStructure(name=base, 
                                           struct=struct,
                                           vectors=[vec0, vec1],
                                           initial_vec=initial_vec))

###############################################################################

# Contador de iterações para logging
iterations = 0

# Para cada estrutura periódica
for p_struct in p_structs:
    
    # Caminho de saída raíz
    output_path_root = dirs_data["h2_gen_output"] / "periodic"

    # Densidade máxima de H2
    max_H2_density = 0.16

    # Número de H2 baseado na área da célula
    area = np.linalg.norm(np.cross(p_struct.vectors[0], p_struct.vectors[1]))

    # Máximo de H2, arredondando e transformando em inteiro
    max_H2_count = int(np.ceil(max_H2_density * area))

    # Gerar para cada quantidade de H2
    for H2_count in range(max_H2_count + 1):
        
        # Nome do arquivo .xyz de saída
        output_name = f"{p_struct.name}-{f'{H2_count:03d}'}{'.xyz'}"

        ############### Random Rot Mono #######################################

        output_path = output_path_root / "random-mono" / output_name
        
        gen_H2_periodic(p_struct, H2_count, output_path, random_rot=True, 
                        both_sides=False, plot=False)

        ############### Vertical Mono #########################################

        # Caminho de saída
        output_path = output_path_root / "vertical-mono" / output_name

        # Gerar estrutura .xyz da estrutura periódica com "i" moléculas de H2
        gen_H2_periodic(p_struct, H2_count, output_path, vertical=True, 
                        both_sides=False, plot=False)
        
        ############### Horizontal Mono #######################################

        output_path = output_path_root / "horizontal-mono" / output_name

        gen_H2_periodic(p_struct, H2_count, output_path, vertical=False, 
                        both_sides=False, plot=False)
        
        # Log
        iterations += 1
        print("\n\n", p_struct.name, "|", "mono: ", iterations, "\n\n")

    # Gerar até o dobro da quantidade de H2
    for H2_count in range(2*max_H2_count + 1):

        # Nome do arquivo .xyz de saída
        output_name = f"{p_struct.name}-{f'{H2_count:03d}'}{'.xyz'}"
                
        ############### Random Rot Dual #######################################

        output_path = output_path_root / "random-dual" / output_name
        
        gen_H2_periodic(p_struct, H2_count, output_path, random_rot=True, 
                        both_sides=True, plot=False)

        ############### Vertical Dual #########################################

        output_path = output_path_root / "vertical-dual" / output_name

        gen_H2_periodic(p_struct, H2_count, output_path, vertical=True, 
                        both_sides=True, plot=False)
        
        ############### Horizontal Dual #######################################

        output_path = output_path_root / "horizontal-dual" / output_name
        
        gen_H2_periodic(p_struct, H2_count, output_path, vertical=False, 
                        both_sides=True, plot=False)
        
        #######################################################################

        # Log
        iterations += 1
        print("\n\n", p_struct.name, "|", "dual: ", iterations, "\n\n")
        
###############################################################################