
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import structure, atom_data
from dopings.config import dops_data, dirs_data
from dopings.h2_gen import h2_gen

# Externas
import numpy as np

###############################################################################

def gen_H2_flake(struct, nH2, output_path, both_sides=False, vertical=False, 
                 plot=False):

    # Movendo pro CM
    struct = h2_gen.to_CM(struct)

    # Gerando arrays
    xs, ys, zs = h2_gen.gen_arrays(struct)

    # Gerando coeficientes
    C = h2_gen.lin_reg(xs, ys, zs, n=3)

    # Gerar H2 - 2D
    R2_H2_coords = h2_gen.gen_R2_coord_flake(struct, nH2)
    
    # Passar pro R3
    R3_H2_coords = h2_gen.gen_R3_H2(R2_H2_coords, C, 
                                        otherside=False, 
                                        vertical=vertical)

    ################################### GERANDO PRO OUTRO LADO

    if both_sides:

        # Gerar outros H2 - 2D
        R2_H2_coords = h2_gen.gen_R2_coord_flake(struct, nH2)
        
        # Passar pro R3
        R3_H2_coords_otherside = h2_gen.gen_R3_H2(R2_H2_coords, C, 
                                                      otherside=True, 
                                                      vertical=vertical)
        
        R3_H2_coords += R3_H2_coords_otherside

    # Colocando na estrutura
    for coord in R3_H2_coords:
        struct.append(atom_data(elem="H", coord=np.round(coord, decimals=8).tolist()))

    # FRAME
    struct.frame(output_path)
    
    if plot:

        # Gerando arrays
        coords = h2_gen.gen_arrays(struct)
        h2_gen.plot(*coords, radius=h2_gen.graphine_radius(struct), C=C)

###############################################################################

# Número máximo de H2 para cada estrutura
nH2_dict = {

    "g1" : 30,
    "g2" : 50,
    "g3" : 80,
    "g4" : 110,
    "g5" : 150
}

# Dicionário de estruturas
struct_dict = {}

# Para cada base
for base in dops_data["graphine"]["bases"]:

    # Lê a estrutura que receberá os H2
    struct = structure(dir=dirs_data["bases"]/"graphine"/base, 
                       read_from_dir=True)
    
    # Lista a estrutura em dicionário
    struct_dict[base] = struct

    # Do 0 ao máximo de H2 que cada estrutura receberá
    for nH2 in range(nH2_dict[base]+1):
        
        # Gerar floco com H2 no diretório especificado
        gen_H2_flake(struct, 
                     nH2, 
                     dirs_data["h2_gen_output"] / "flake" / f"{base}-{nH2:03d}.xyz", 
                     plot=False,
                     vertical=True,
                     both_sides=True)
