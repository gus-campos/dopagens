
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import Structure
from dopings.atom import Atom
from dopings.config import dops_data, dirs_data
from dopings.h2_gen import H2Gen

# Externas
import numpy as np

###############################################################################

def gen_H2_flake(struct, nH2, output_path, both_sides=False, vertical=False, 
                 plot=False):
    """
    Gera um arquivo .xyz da estrutura já adicionada de moléculas de H2.

    Parameters 
    ----------

    struct : Structure
        Objeto que representa uma estrutura.

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
    struct = H2Gen.to_CM(struct)

    # Gerando coordenadas de pontos distribuídos no plano XY
    R2_H2_coords = H2Gen.gen_R2_coord_flake(struct, nH2)

    # Gerando coeficientes que ajustam a curva aos átomos
    struct_coords = H2Gen.gen_arrays(struct)
    C = H2Gen.lin_reg(*struct_coords, n=3)
    
    # Transformando pontos distribuídos no plano, em átomos de H distribuídos 
    # sobre a estrutura
    R3_H2_coords = H2Gen.gen_R3_H2(R2_H2_coords, C, otherside=False, 
                                    vertical=vertical)

    ########################## GERANDO PRO OUTRO LADO #########################

    # Se for para gerar H2 nos dois lados da estrutura
    if both_sides:

        # Gerar mais pontos espalhados no plano
        R2_H2_coords = H2Gen.gen_R2_coord_flake(struct, nH2)
        
        # Gerar pontos espalhados sobre o outro lado da estrura
        R3_H2_coords_otherside = H2Gen.gen_R3_H2(R2_H2_coords, C, 
                                                  otherside=True, 
                                                  vertical=vertical)
        
        # Juntar listas de H dos dois lados
        R3_H2_coords += R3_H2_coords_otherside

    ###########################################################################

    # Adicionar todos os átomos de H na estrutura
    for coord in R3_H2_coords:
        struct.append(Atom(elem="H", coord=np.round(coord, decimals=8).tolist()))

    # Gerar aquivo .xyz da estrutura
    struct.frame(output_path)
    
    # Se for para plotar, gerar visualização interativa da estrutura com curva
    if plot:
        radius = H2Gen.graphine_radius(struct)
        H2Gen.plot(*struct_coords, radius=radius, C=C)

###############################################################################

# Dicionário do número máximo de H2 que deve ser gerado para cada estrutura
nH2_dict = {

    "g1" : 30,
    "g2" : 50,
    "g3" : 80,
    "g4" : 110,
    "g5" : 150
}

# Para cada base
for base in dops_data["graphine"]["bases"]:

    # Lê a estrutura que receberá os H2
    struct = Structure(dir=dirs_data["bases"] / "graphine" / base, 
                       read_from_dir=True)

    # Do 0 ao máximo de H2 que cada estrutura receberá
    for nH2 in range(nH2_dict[base]+1):
        
        # Nome do arquivo .xyz de saída
        output_name = f"{base}-{nH2:03d}.xyz"

        # Caminho de saída
        output_path = dirs_data["h2_gen_output"] / "flake" / output_name

        # Gerar floco com H2 no diretório especificado
        gen_H2_flake(struct, nH2, output_path, plot=False, vertical=True,
                     both_sides=True)
