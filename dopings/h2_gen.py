import numpy as np

def phi(i: int, n: int, x, y):
    """
    Retorna a imagem da n-ésima função phi, tal que:

    phi{0}               = 1
    phi{1}   ... phi{n}  = x^n
    phi{n+1} ... phi{2n} = y^n

    Parameters
    ----------

    i : int
        Índice da i-ésima função phi, ou seja phi{i}.

    x : float or np.ndarray
        Valor único de x, ou array de valores de x.
        Se array, precisa ter o mesmo tamanho de y.

    y : float or np.ndarray
        Valor único de y, ou array de valores de y.
        Se array, precisa ter o mesmo tamanho de x.

    Returns
    -------

    z : float or np.ndarray
        Imagem da função phi{i}, ou array de imagens.
    
    """

    if i == 0:
        return 1
    
    elif i <= n:
        return x**i

    elif i > n and i <= 2*n:
        return y**(i-n)
    
    else:
        raise ValueError("'i' must be less or equal to 2n.")

def g(x: float, y: float, C: np.ndarray):
    """
    Dado uma matrix coluna de coeficiente do ajuste, calcula a iamgem da
    função:

    g(x,y) = C0 * phi{0}(x,y) + ... + C(2n+1) * phi{2n}(x,y)

    tal que:

    phi{0}               = 1
    phi{1}   ... phi{n}  = x^n
    phi{n+1} ... phi{2n} = y^n

    Parameters
    ----------

    x : float
        Valor de x.

    y: float
        Valor de y.

    C : np.ndarray
        Matriz coluna com os coeficientes do ajuste.

    Returns
    -------

    z : float
        Imagem da função.
    """
    
    # Encontrando o valor n, a partir da dimensão da matriz de coeficientes
    n = (C.shape[0] - 1) // 2

    # Somando as parcelas da função ajustada
    z = 0
    for i in range(2*n+1):
        z += C[i] * phi(i, n, x, y)

    # Retornando imagem
    return z

def struct_radius(struct: "Structure"):

    """"
    Estima um raio aproximado para a estrutura de grafino.
    
    Encontra o carbono mais afastado do centro de massa, e retorna 90% 
    desta distância. Assume que o centro de massa da estrutura está 
    em (0,0,0).

    Parameters
    ----------

    struct : strucutre
        Estrutura para qual se calcular a estimativa.

    Returns
    -------

    float
        Raio estimado para a estrutura.
    """

    maior = 0
    for atom in struct.atoms:
        if atom.elem == "C":
            length = np.linalg.norm(atom.coord)
            if length > maior:
                maior = length

    # Rotrnar a estimativa           
    return maior

class H2Gen:
    """
    Classe de métodos estáticos para geração de estruturas cobertas por
    moléculas de H2.

    Methods
    -------
    
    plot(xs, ys, zs, radius: float, C: np.ndarray=None)

        Gera um plot interativo de todos os átomos, dadas as coordenadas 
        de todos eles. Se for passada uma matriz de coeficientes, plota 
        também a curva ajustada.
    
    to_CM(struct : "Structure")

        Recebe uma estrutura e retorna a mesma movida para o seu centro 
        de massa. Só considera os carbonos no cálculo do centro de 
        massa.
    
    gen_arrays(struct : "Structure")

        Gera arrays xs, ys, zs das coordenadas dos átomos de uma 
        estrutura.
    
    lin_reg(xs, ys, zs, n: int)

        Calcula  os coeficientes que ajustam, aos pontos fornecidos, uma 
        curva do tipo:

        g(x,y) = C0 * phi{0}(x,y) + ... + C(2n+1) * phi{2n}(x,y)
    
    struct_radius(struct: "Structure")

        Estima um raio aproximado para a estrutura de grafino.
        
        Encontra o carbono mais afastado do centro de massa, e retorna 
        90% desta distância. Assume que o centro de massa da estrutura 
        está em (0,0,0).
    
    gen_R2_coords_periodic(periodic_struct: "PeriodicStructure",
                           nH2: int)

        Gera, para uma estrutura periódica, vários pontos no R2, 
        dispostos em uma grade, que correspondem aos centros das 
        moléculas de H2.
    
    gen_R2_coord_flake(struct: "strucutre", nH2: int)

        Gera, para uma estrutura de formato arredondado, no R2, uma 
        lista de com posições de H2.
    
    gen_R3_H2(R2_coords, C, otherside=False, vertical=False)

        Recebe pontos espalhados sobre um plano, e retorna pontos 
        espalhados sobre uma superfície ajustada sobre os átomos da 
        estrutura.
    """

    ##################### Auxiliar ################################################
    @staticmethod
    def plot(xs, ys, zs, radius: float, C: np.ndarray=None):
        """
        Gera um plot interativo de todos os átomos, dadas as coordenadas de 
        todos eles. Se for passada uma matriz de coeficientes, plota também
        a curva ajustada.

        Parameters 
        ----------

        xs, ys, zs : np.array[float]
            Array de floats das coordenadas dos átomos.

        radius : float
            Raio estimado da estrutura.

        C : np.ndarray[float], default = None
            Matriz coluna de coeficientes do ajuste de curva.
        """

        import matplotlib.pyplot as plt
        
        ######################### ÁTOMOS ##########################################

        plt.figure()
        ax = plt.subplot(111, projection='3d')
        ax.scatter(xs, ys, zs, color='b')

        ######################### CURVA ###########################################
        if C is not None:
        
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            # Malha para o plot do plano
            X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]), np.arange(ylim[0], ylim[1]))
            Z = np.zeros(X.shape)

            # Preenchendo valores de Z da malha
            for r in range(X.shape[0]):
                for c in range(X.shape[1]):
                    Z[r,c] = g(X[r][c], Y[r][c], C)

            # Plotando a grade
            ax.plot_wireframe(X,Y,Z, color='k')

        ################################# PLOTAGEM ################################

        # Limite dos eixos
        ax.set_xlim(-radius, radius)
        ax.set_ylim(-radius, radius)
        ax.set_zlim(-radius, radius)

        # Mostrando plot
        plt.show()

        # Fechando
        plt.close()
        plt.clf()

    @staticmethod
    def to_CM(struct : "Structure"):
        """
        Recebe uma estrutura e retorna a mesma movida para o seu centro 
        de massa. Só considera os carbonos no cálculo do centro de 
        massa.

        Parameters
        ----------

        struct : Structure
            Estrutura a ser transladada.

        Returns
        -------

        Structure 
            Estrutura transladada para o centro de massa.
        """

        import numpy as np

        # Copia a Structure para não alterar a original
        struct = struct.copy()

        # Inicializando variáveis
        center_of_mass = np.array([0,0,0], dtype=float)
        n_atoms = 0

        # Para cada átomo de carbono
        for atom in struct.atoms:
            if atom.elem == "C":

                # Incrementar soma das posições e quantidade de carbonos
                center_of_mass += atom.coord
                n_atoms += 1

        # Calcular centro de massa
        center_of_mass = center_of_mass / n_atoms

        # Transladar cada átomo da estrutura
        for i in range(struct.size):
            struct.atoms[i].coord -= np.round(center_of_mass, decimals=8)

        return struct

    @staticmethod
    def gen_arrays(struct : "Structure"):
        """
        Gera arrays xs, ys, zs das coordenadas dos átomos de uma 
        estrutura.

        Parameters
        ----------

        struct : Structure
            Estrutura para a qual se deseja gerar os arrays.

        Returns
        -------

        xs, ys, zs : np.array[float]
            Arrays das coordenadas x, y e z dos átomos da estrutura.
        """

        import numpy as np

        # Listas
        xs, ys, zs = [], [], []

        # Dividindo as váriaveis em arrays diferentes
        for atom in struct.atoms:

            xs.append(atom.coord[0])
            ys.append(atom.coord[1])
            zs.append(atom.coord[2])

        # Transformando em array
        xs = np.array(xs, dtype=float)
        ys = np.array(ys, dtype=float)
        zs = np.array(zs, dtype=float)

        # Retornando
        return xs, ys, zs

    #######################  Ajuste de curva ######################################

    @staticmethod
    def lin_reg(xs, ys, zs, n: int):
        """
        Calcula  os coeficientes que ajustam, aos pontos fornecidos, uma 
        curva do tipo:

        g(x,y) = C0 * phi{0}(x,y) + ... + C(2n+1) * phi{2n}(x,y)

        Parameters
        ----------

        xs, ys, zs : np.ndarray
            Arrays com as coordenadas dos pontos/átomos. 

        Returns
        -------

        np.ndarray(shape=(2n+1))
            Matriz coluna de coeficiente do ajuste de curva, na forma:

        |   C0    | 
        |   C1    | 
        |    .    | 
        |    .    | 
        |    .    | 
        | C(2n+1) |
        """

        # Inicializando matrizes
        A = np.zeros(shape=(2*n+1, 2*n+1), dtype=float)
        B = np.zeros(shape=(2*n+1       ), dtype=float)

        # Para cada item da matriz A, calcular o produto escalar correspondente
        for i in range(2*n+1):
            for j in range(2*n+1):

                A[i,j] = np.sum(phi(i, n, xs, ys) * phi(j, n, xs, ys))

        # Para cada item da matriz B, calcular o produto escalar correspondente
        for i in range(2*n+1):
            B[i] = np.sum(phi(i, n, xs, ys) * zs)

        # Resolve o sistema linear, retornando sua solução, ou seja, retornando
        # a matriz coluna de coeficientes do ajuste
        return np.linalg.solve(A, B)

    #######################  R2 GEN  ##############################################

    @staticmethod
    def gen_R2_coords_periodic(periodic_struct: "PeriodicStructure", 
                            nH2: int):
        
        """
        Gera, para uma estrutura periódica, vários pontos no R2, 
        dispostos em uma grade, que correspondem aos centros das 
        moléculas de H2.

        Parameters
        ----------

        periodic_struct : PeriodicStructure
            Estrutura periódica, dotada de vectores de célula.

        nH2 : int
            Quantidade de moléculas de H2 que deseja gerar para tal 
            estrutura.
        
        Returns
        -------

        list[np.ndarray]
            Lista de vetores no R2, com a coordenada dos centros gerados 
            para as moléculas de H2. 

        """

        # Encontrando menor quadrado perfeito, maior que nH2
        n = 0
        while n**2 < nH2:
            n += 1

        # Se n for zero, retornar lista vazia
        if n == 0:
            return []
        
        # Lista de centros
        H2_centers = []

        # Dimensão da célula da grade usada
        step = 1/n

        # Malha de pontos discretos
        X,Y = np.meshgrid(np.arange(0, 1, step), np.arange(0, 1, step))
        Z = np.zeros(X.shape, dtype=bool)

        # Vetores da cécula, a partir da exclusão da coordenada Z do vetores
        # da célula
        vec0 = np.delete(periodic_struct.vectors[0], 2)
        vec1 = np.delete(periodic_struct.vectors[1], 2)

        # Offset dos átomos em relação à origem
        # Isso evita que as moléculas fiquem na borda da célula
        offset = (vec0 * (step/2)) + (vec1 * (step/2))

        # Criando lista com todos os pontos da grade
        points = []
        for r in range(n):
            for c in range(n):
                points.append((r,c))

        # De todos os pontos da grade sorteando apenas a quantidade necessária
        import random
        chosen_points = random.sample(points, nH2)

        # Para cada ponto escolhido, gerar as coordenadas de seu centro
        for r, c in chosen_points:

            # Gerando centro a partir da combinação linear dos vetores da célula
            center = (X[r,c] * vec0 +  Y[r,c] * vec1) + offset

            # Listando centros
            H2_centers.append(center)

        return H2_centers

    @staticmethod
    def gen_R2_coord_flake(struct: "strucutre", nH2: int, 
                           main_vector: np.ndarray=None):
        """
        Gera, para uma estrutura de formato hexagonal, no R2, uma lista 
        de com posições de H2.

        Parameters
        ----------

        struct : Structure
            Estrutura para a qual se deseja gerar os centros de H2.

        nH2 : int 
            Número de H2 desejado na molécula.

        main_vector : np.ndarray[float]
            Vetor no R2 com a direção de um dos vértices do hexágono.

        Returns
        -------
        
        list[np.ndarray]
            Lista, de vetores no R2, com as posições geradas para os 
            centros de H2. 
        """

        def tri_nums(n: int):

            """n-ésimo número triangular"""

            return int((n*(n+1))/2)

        def smallest_triangular(nH2):

            """
            O menor m-ésimo número triangular, que quando usado para 
            subdividir todos os triângulos, permite conter todas as nH2
            moléculas de H2 pretendidas.
            """

            m = 0
            while (6*tri_nums(m) < nH2):
                m += 1
            return m

        def triangle_centers(vecs, n):

            """
            Divide o triângulo formado pelos vetores 'vecs' em tri_nums(n) 
            minitriângulos e retorna uma lista com seus centros.

            Parameters
            ----------

            vecs : list[np.ndarray]
                Lista de dois vetores que geram o triângulo em questão.

            n : int
                refere o n-ésimo número triangular a ser usado na subsivisão
                do triângulo.

            Returns
            -------

            list[np.ndarray]
                Lista de todos os centros possíveis, em dada subdivisão, do
                triângulo.
            """

            centers = []

            # Vetores do triângulo unitário
            vecs_unit = [vecs[0]/n, vecs[1]/n]

            # offset
            offset = 0.03 * (vecs_unit[0] + vecs[1])

            # Centro de massas do primeiro triângulo
            cm_unit = ((vecs_unit[0] + vecs_unit[1]) / 3)

            # Para cada triângulo unitário
            for k in range(0, n):
                for l in range(0, n):

                    # Apenas se tiver dentro do triângulo maior
                    if k+l < n:
                        
                        # Translação
                        transl = k*vecs_unit[0] + l*vecs_unit[1]

                        # Centro de massa transladado
                        centers.append(offset + cm_unit + transl)

            return centers

        def rotate_60(vector) -> np.ndarray:
            """
            Rotaciona um vetor no R2 em 60 graus.

            Adaptado de: https://stackoverflow.com/questions/20840692/rotation-of-a-2d-array-over-an-angle-using-rotation-matrix 
            """
            vector = np.array(vector)

            vector = vector.T

            theta = np.radians(60)

            rotation_matrix = np.array([
                [np.cos(theta), -np.sin(theta)],
                [np.sin(theta), np.cos(theta)]
            ])

            return np.matmul(rotation_matrix, vector).T

        import random

        # Se nenhum nH2 for gerado, retornar lista vazia
        if nH2 == 0:
            return []

        # Raio estimado
        radius = struct_radius(struct) #* 0.95

        # Determinando dimensão dos triângulos
        n_triangles = smallest_triangular(nH2)

        # Rotacionando o vetor principal da estrutura
        hex_vecs = []

        # Adicionando À lista o primeiro vetor
        if main_vector is not None:
            hex_vecs.append(radius * (main_vector / np.linalg.norm(main_vector)))
        else:
            # Se não foi passado, usar (1, 0)
            hex_vecs.append(radius * np.array([1.0, 0.0]))

        for i in range(5):
            hex_vecs.append(rotate_60(hex_vecs[-1]))

        # Listando os 6 pares de vetores, que geram os 6 triângulos principais
        tri_vecs = []
        for i in range(6):
            tri_vecs.append([hex_vecs[i%6], hex_vecs[(i+1)%6]])

        # Quantidade de moléculas em cada triângulo principal
        unit_samples = [nH2 // 6] * 6
        for i in range(nH2 % 6):
            unit_samples[i] += 1

        if sum(unit_samples) != nH2:
            raise AssertionError("Sum of samples dosn't match nH2")

        # Para cada triângulo, gerar todos os centros possíveis, mas escolher 
        # apenas a quantidade necessária
        centers = []
        for i, vec_pair in enumerate(tri_vecs):
            for center in random.sample(

                    triangle_centers(vecs=vec_pair, n=n_triangles), unit_samples[i]
                ):
        
                centers.append(center)

        return centers

    #######################  Posicionador no R3  ##################################

    @staticmethod
    def gen_R3_H2(R2_coords, C, otherside=False, vertical=False):
        """
        Recebe pontos espalhados sobre um plano, e retorna pontos 
        espalhados sobre uma superfície ajustada sobre os átomos da 
        estrutura.

        Parameters
        ----------

        R2_coord : list[np.ndarray]
            Lista de vetores no R2 que dão a posição dos centros de H2 
            dispostos no plano.

        C : np.ndarray
            Matriz coluna de coeficientes da função ajustada sobre os 
            pontos da estrutura.

        otherside : bool
            Se é desejado gerar tais pontos no lado oposto da estrutura,
            ou seja, no lado onde z é negativo.

        vertical : bool
            Se é deseja que o eixo da molécula de H2 esteja orientada 
            perpendicular à estrutura. Se False, o eixo será parelelo à
            estrutura, e sua rotação será aleatória.

        Returns
        -------
        
        list[np.ndarray]
            Lista com vetores no R3 das posições de cada átomo de H, 
            com posição ajustada de acordo com a curvatura da estrutura.
        """

        import math

        plane_distance = 2.7
        H2_bond_length = 0.7
        
        if vertical:
            plane_distance += H2_bond_length / 2

        # Posições do H
        H_coords = []

        # Para cada posição na lista
        for R2_posic in R2_coords:

            # Centro no R3
            center = np.array([R2_posic[0], R2_posic[1], 0.0])

            # Componente Z  
            if otherside:
                z = g(center[0], center[1], C) - plane_distance

            else:
                z = g(center[0], center[1], C) + plane_distance
            
            # Vetor distância do plano 
            plane_dist_vec = np.array([0.0, 0.0, z], dtype=float)
            
            # Meio eixo do H2
            if vertical:
                alpha = np.radians(45)
                H2_axis = (H2_bond_length / 2) * np.array([0.0, 0.0, 1.0])

            else:
                alpha = np.random.uniform(0, 2*np.pi)
                H2_axis = (H2_bond_length / 2) * np.array([math.cos(alpha), math.sin(alpha), 0.0])

            # Posição dos dois H - centro, mais distância do plano, +/- meio eixo do H2
            posic1 = center + plane_dist_vec + H2_axis 
            posic2 = center + plane_dist_vec - H2_axis

            # Adicionar à estrutura
            H_coords.append(posic1)
            H_coords.append(posic2)

        # Retornar como array
        return H_coords

    ###############################################################################