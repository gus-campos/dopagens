import numpy as np
from dopings.structure import Structure
from dopings.atom import Atom

class unit_cell:

    def __init__(self, basis, origin=None, dims=None, unit_grid=None, nC=None):

        if isinstance(basis, list) and all([ isinstance(base, np.ndarray) for base in basis ]):
            self.basis = basis
        else:
            raise TypeError("basis")

        if (origin is not None):
            if isinstance(origin, np.ndarray):
                self.origin = origin
            else:
                raise TypeError("origin")

        else:
            self.origin = np.zeros(2)

        if (dims is not None): 
            
            if isinstance(dims, tuple) and all([ isinstance(dim, int) for dim in dims]):
                self.dims = dims
                self.steps = (1/dims[0], 1/dims[1])
            else:
                raise TypeError("dims")
        else:
            self.dims = None

        if (unit_grid is not None):

            if isinstance(unit_grid, tuple):
                self.unit_grid = unit_grid
            else:
                raise TypeError("unit_grid")
        else:
            self.unit_grid = None
        
        if (nC is not None):

            if isinstance(nC, int):
                self.nC = nC
            else:
                raise TypeError("nC")
        else:
            self.nC = None

    def random_grid(self):

        # Malha para pontos discretos
        X,Y = np.meshgrid(np.arange(0, 1, self.steps[0]), np.arange(0, 1, self.steps[1]))

        valid = False

        while not valid:

            valid = True
            Z = np.zeros(X.shape, dtype=bool)

            # Sorteando pontos
            for k in range(self.nC):

                first = True

                # Enquanto estiver sorteando pontos repetidos, refazer
                while first or Z[i,j] == True:

                    first = False

                    # Não sei pq tive que inverter os valores de "i" e "j" - Não considerar última linha e coluna
                    i = np.random.randint(low=0, high=self.dims[1], size=1)[0]
                    j = np.random.randint(low=0, high=self.dims[0], size=1)[0]

                Z[i,j] = True

            def f(i, j, validate=False):

                valid = True

                if i == -1:
                    i = Z.shape[0]-1
                    valid = False

                elif i == Z.shape[0]:
                    i = 0
                    valid = False

                if j == -1:
                    j = Z.shape[1]-1
                    valid = False

                elif j == Z.shape[1]:
                    j = 0
                    valid = False

                if validate:
                    return valid

                return i, j

            def which_side(i,j):

                sides = []

                if i == -1:
                    sides.append(1)

                if i == Z.shape[0]:
                    sides.append(2)

                if j == -1:
                    sides.append(3)

                if j == Z.shape[1]:
                    sides.append(4)

                return sides

            sides = []

            # Verificar se todos elementos têm vizinhos
            for i in range(Z.shape[0]):
                for j in range(Z.shape[1]):

                    # Para cada átomo colocado
                    if Z[i,j]:

                        near = False

                        if f(i-1, j-1, True) and Z[i-1, j-1]:
                            sides += which_side(i-1,j-1)

                        if f(i-1, j+1, True) and Z[i-1, j+1]:
                            sides += which_side(i-1,j+1)

                        if f(i+1, j-1, True) and Z[i+1, j-1]:
                            sides += which_side(i+1,j-1)

                        if f(i-1, j+1, True) and Z[i+1, j+1]:
                            sides += which_side(i+1,j+1)
                    
                        if not near:
                            valid = False

            # Tem vizinhos em todos os lados
            if sum(list(set(sides))) != 10:
                valid = False

        self.unit_grid = (X, Y, Z)

    def gen_translate(self, i, j, basis=None):

        if basis is None:
            basis = self.basis

        new_origin = self.origin + i*basis[0] + j*basis[1]

        return unit_cell(basis=self.basis, 
                         origin=new_origin, 
                         dims=self.dims, 
                         unit_grid=self.unit_grid, 
                         nC=self.nC)

    def gen_rotate(self, theta):

        # Convertendo pra radianos
        theta = np.radians(theta) 

        # Matriz de rotação
        c, s = np.cos(theta), np.sin(theta)
        R = np.array(((c,-s), (s, c)))

        # Rotacionando
        new_basis = [None, None]
        new_basis[0] = np.dot(R, self.basis[0])
        new_basis[1] = np.dot(R, self.basis[1])

        return unit_cell(basis=new_basis, 
                         origin=self.origin, 
                         dims=self.dims, 
                         unit_grid=self.unit_grid, 
                         nC=self.nC)


class cells_set:

    def __init__(self, cells_list):

        self.cells_list = cells_list

    def append(cells):
        
        if isinstance(cells, unit_cell):
            self.cells_list.append(cells)
        
        elif isinstance(cells, list) and all([ isinstance(cell, unit_cell) for cell in cells ]):
            self.cells_list += cells
        
    def gen_flake(self):

        atoms = []

        for cell in self.cells_list:

            X, Y, Z = cell.unit_grid
            origin = cell.origin
            basis = cell.basis
            dims = cell.dims

            # Gerar átomos a partir da grade e dos vetores
            for r in range(dims[1]):
                for c in range(dims[0]):

                    # Se foi escolhida tal posição
                    if Z[r,c] == True:

                        # Gerando ponto a partir da combinação linear dos vetores da célula
                        center = origin + X[r,c] * basis[0] + Y[r,c] * basis[1]

                        # Listando centro
                        atoms.append(Atom(elem="C", coord=np.append(center, [0.0]).tolist()))

        return Structure(atoms=atoms)

########################################################################

# Densidade de carbonos no grafeno
mi = 0.3917                         # Átomos de carbono por ângstron quadrado
theta = np.radians(360 / 6)         # Ângulo entre vetores
dims = (4, 4)                       # Dimensões do lado da grade, em átomos de carbono
C_bond_length = 1.3                 # Distância considerada

# Norma dos vetores de acordo com a distância desejada entre vizinhosa
norms = (C_bond_length * dims[0], C_bond_length * dims[1])

# Matriz de rotação do primeiro pro segundo vetor
c, s = np.cos(theta), np.sin(theta)
R = np.array(((c,-s), (s, c)))

# Vetores da célula
basis = [None, None]
basis[0] = np.array([1.0, 0.0], dtype=float)
basis[1] = np.dot(R, basis[0])


# Area da célula, excluindo as bordas
area = np.linalg.norm(np.cross(basis[0] * norms[0], basis[1] * norms[1]))

# Número de carbonos, de acordo com a densidade
from math import floor
nC = floor(area * mi)

# Escalonando
basis[0] = norms[0] * basis[0]
basis[1] = norms[1] * basis[1]

########################################################################

for i in range(10):

    cell = unit_cell(basis=basis, dims=dims, nC=nC)
    cell.random_grid()

    # Célula inicial
    cells = cells_set(cells_list=[cell])
    cells.gen_flake().frame(f"TESTES/g{i:02d}a.xyz")

    cells_list = []

    if True:
        n = 6
        for k in range(n):
            for l in range(n):
                cells_list.append(cell.gen_translate(k,l))

    cells = cells_set(cells_list=cells_list)
    cells.gen_flake().frame(f"TESTES/g{i:02d}b.xyz")