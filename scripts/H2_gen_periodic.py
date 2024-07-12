
# Incluindo endereço do diretório mãe
import sys
sys.path.append("../" )

# Bibliotecas internas
from dopings.structure import structure, atom_data
from dopings.config import dirs_data, h2_gen_data
from h2_gen.h2_gen import h2_gen

# Externas
import numpy as np

###############################################################################

def gen_H2_periodic(p_struct, nH2, output_path, both_sides=False, vertical=False, 
                    plot=False):

    # Copiando structure
    struct = p_struct.struct.copy()

    # Movendo pro CM
    struct = h2_gen.to_CM(struct)

    # Gerando arrays
    xs, ys, zs = h2_gen.gen_arrays(struct)

    # Gerando coeficientes
    C = h2_gen.lin_reg(xs, ys, zs, n=1)

    # Gerar H2 - 2D
    R2_H2_coords = h2_gen.gen_R2_coords_periodic(p_struct, nH2)
    
    # Passar pro R3
    R3_H2_coords = h2_gen.gen_R3_H2(R2_H2_coords, C, 
                                    otherside=False, 
                                    vertical=vertical)

    ################################### GERANDO PRO OUTRO LADO

    if both_sides:

        # Gerar outros H2 - 2D
        R2_H2_coords = h2_gen.gen_R2_coords_periodic(p_struct, nH2)
        
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

# Classe de estruturas periódicas
class periodic_structure:

    def __init__(self, name=None, nH2=None, vectors=None):

        self.name=name
        self.struct=None
        self.nH2=nH2
        self.vectors=vectors

# Listando estruturas periódicas
p_structs = []
for name in ["g1_s1", "g1_s2", "g1_s3", "g1_s4"]:
    p_structs.append(
        
        periodic_structure(name=name, 
                           nH2=h2_gen_data[name]["nH2"],
                           vectors=[np.array(h2_gen_data[name]["vectors"][0]),
                                    np.array(h2_gen_data[name]["vectors"][1])])
        )

# Lendo estruturas do armazenamento
for i, p_struct in enumerate(p_structs):
    p_structs[i].struct = structure(dir=dirs_data["bases"] / "g1-periodic" / p_struct.name, 
                                    read_from_dir=True)

###############################################################################

# Gerando todas as versões de cada estrutura periódica
for p_struct in p_structs:
    for i in range(p_struct.nH2 + 1):
        
        output_name = f"{p_struct.name}-{f'{i:03d}'}{'.xyz'}"

        gen_H2_periodic(p_struct,
                        nH2=i,
                        output_path=dirs_data["h2_gen_output"]/ "periodic" / output_name,
                        vertical=True,
                        both_sides=True,
                        plot=False)

###############################################################################

