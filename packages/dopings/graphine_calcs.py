
import numpy as np
from dopings.config import graphine_data, dops_data, atoms_data
from dopings.structure import structure
from dopings.atom_data import atom_data

# Dados de sítios extras
sites_id = dops_data["graphine"]["sites_id"]
extra_sites_id = graphine_data["sites_id"]

def dihedral_angle(atoms: list[atom_data]) -> float:
    """
    Calcula o ângulo diedral, ou ângulo de torção, entre 4 átomos.

    Código adaptado de:
    https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    
    Parameters
    ----------

    atoms : list[atom_data]
        Lista de 4 átomos sobre os quais deseja-se calcular a 
        torção.

    Returns
    -------

    float
        Ângulo diedral entre os 4 átomos.
    """    

    p0 = atoms[0].coord
    p1 = atoms[1].coord
    p2 = atoms[2].coord
    p3 = atoms[3].coord

    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

class graphine_calcs:
    """
    Classe de métodos estáticos que trás implementações de calculos 
    sobre estruturas dopadas de grafinos.

    Methods
    -------

    axis_torsion_angle(struct: structure, arm: int) -> float
        Calcula a torção de um eixo da estrutura, usando o braço "arm"
        e o seu oposto para cacular um ângulo diedral.
    
    int_arm_length(struct: structure) -> float
        Calcula a soma de distâncias entre os carbonos no braço interno
        da dopagem, ou seja, no braço 0, dentre os seguintes sítios:

        [A1, B1, B2, B3, ..., B{2N}, D4]
    
    ext_arm_length(struct: structure) -> float
        Calcula a soma de distâncias entre os carbonos no braço externo
        da dopagem, ou seja, no braço 0, dentre os seguintes sítios:

        [C2, D1, D2, ..., D{2N+1}]
    
    dop_atom_angle(struct: structure) -> float|None
        Se o átomo dopante pertencer a um sítio B ou D, calcula o 
        ângulo da ligação do átomo dopado com seus dois vizinhos.
    
    dop_atom_torsion(struct: structure) -> float|None
        Quando o dopante se encontra em um anel benzênico, calcula a 
        torção de um sequência de 4 átomos, onde o dopante é o segundo.
    """

    @staticmethod
    def axis_torsion_angle(struct: structure, arm: int) -> float:
        """
        Calcula a torção de um eixo da estrutura, usando o braço "arm"
        e o seu oposto para cacular um ângulo diedral.

        Parameters
        ----------

        struct : strucutre
            Estrutura de grafino sobre a qual se deseja fazer o cálculo.

        arm : int
            Braço da estrutura usado de referência para o cálculo da 
            torção. Considera o próprio braço, e o seu oposto no 
            cálculo.

        Returns
        -------

        float
            Ângulo de torção no eixo do braço especificado.
        """

        def f(angle):
            """
            Transforma uma torção que parte do 0 para uma torção que parte do 
            180.
            """

            # Deixando todos os ângulos positivos
            if angle < 0:
                angle += 360

            # Retornando sua diferença em relação ao ângulo raso
            return (angle + 360) if (angle < 0) else (angle - 180)

        def arms_seq(initial: int, n: int):
            """
            Partindo do braço "initial", encontra o "n"-ésimo braço na
            sequência, em sentido horário.
            """

            # Trazendo pro intervalo 0, 5
            return (initial + n) % 6

        # Extra_sites
        base = struct.base

        # Próximo sítio D
        next_d = f"D{str(int(base[1:]) + 1)}"

        # Átomos
        atom1 = struct.atoms[extra_sites_id[base][arms_seq(arm, 0)][next_d] - 1]
        atom2 = struct.atoms[extra_sites_id[base][arms_seq(arm, 2)][ "A1" ] - 1]
        atom3 = struct.atoms[extra_sites_id[base][arms_seq(arm, 5)][ "A1" ] - 1]
        atom4 = struct.atoms[extra_sites_id[base][arms_seq(arm, 3)][next_d] - 1]

        angle = dihedral_angle([atom1, atom2, atom3, atom4])

        return angle

    @staticmethod
    def int_arm_length(struct: structure) -> float:
        """
        Calcula a soma de distâncias entre os carbonos no braço interno
        da dopagem, ou seja, no braço 0, dentre os seguintes sítios:

        [A1, B1, B2, B3, ..., B{2N}, D4]

        Parameters
        ----------

        struct : strucutre
            Estrutura de grafino sobre a qual se deseja fazer o cálculo.

        Returns
        -------

        float
            Soma das distâncias entre carbonos no braço interno da 
            dopagem.
        """
        # Inicializando comprimento total como 0
        total_length = 0

        # Base
        base = struct.base

        # Sites_names list
        sites_names = dops_data["graphine"]["sites"][base]

        # Separando sites B
        sites = []
        # Adicionando
        sites.append("A1")

        # Para cada site da lista, se for site B, adicionar a lista B_sites
        for site_name in sites_names:
            if site_name[0:1] == "B":

                # Apenas se for carbono
                if struct.atoms[sites_id[base][site_name] - 1].elem == "C":

                    # Adicionar à lista
                    sites.append(site_name)

        # Adicionar C4
        sites.append("C4")

        # Para cada sítio B, exceto o último
        for i in range(len(sites) - 1):

            # Átomo do i, e do i+1 na lista de sítios
            atom1 = struct.atoms[sites_id[base][sites[i  ]] - 1]
            atom2 = struct.atoms[sites_id[base][sites[i+1]] - 1]

            # Incrementar o comprimento total com a distância entre eles
            total_length += atom1.dist_to(atom2)

        return total_length

    @staticmethod
    def ext_arm_length(struct: structure) -> float:
        """
        Calcula a soma de distâncias entre os carbonos no braço externo
        da dopagem, ou seja, no braço 0, dentre os seguintes sítios:

        [C2, D1, D2, ..., D{2N+1}]

        Parameters
        ----------

        struct : strucutre
            Estrutura de grafino sobre a qual se deseja fazer o cálculo.

        Returns
        -------

        float
            Soma das distâncias entre carbonos no braço interno da 
            dopagem.
        """

        # Inicializando comprimento total como 0
        total_length = 0

        # Base
        base = struct.base

        # Sites_names list
        sites_id_copy = sites_id[base].copy()
        extra_sites_id_copy = extra_sites_id[base][0].copy()

        # Juntando os dois dicionários
        sites_id_copy.update(extra_sites_id_copy)

        # Criando novo dicionário
        sites = {}

        # Para cada chave  no dict original
        for key in sites_id_copy.keys():

            # Se atender ao critério
            if key == "C2" or (key[0:1] == "D" and struct.atoms[sites_id_copy[key] - 1].elem == "C"): 

                # Usar para formar novo par chave-valor
                sites[key] = sites_id_copy[key]
        
        keys = list(sites.keys())

        # Incrementando comprimento total baseado nos sitios separados
        for i in range(len(keys) - 1):

            # Átomo do i, e do i+1 na lista de sítios
            atom1 = struct.atoms[sites[keys[i  ]] - 1]
            atom2 = struct.atoms[sites[keys[i+1]] - 1]

            # Incrementar o comprimento total com a distância entre eles
            total_length += atom1.dist_to(atom2)

        return total_length

    @staticmethod
    def dop_atom_angle(struct: structure) -> float|None:
        """
        Se o átomo dopante pertencer a um sítio B ou D, calcula o 
        ângulo da ligação do átomo dopado com seus dois vizinhos.
        
        Se não for o caso, retorna None.

        Parameters
        ----------

        struct : strucutre
            Estrutura de grafino sobre a qual se deseja fazer o cálculo.

        Returns
        -------

        float 
            O ângulo entre o dopante e os 2 átomos ligados a ele, se o
            dopante for de um sítio B ou D.

        None
            Se o dopante não pertencer a um sítio B ou D.
        """
        def bond_angle(atoms: list[atom_data]) -> float:
            """
            Calcula o ângulo da ligação de uma sequência de 3 átomos.
            """

            vec1 = atoms[0].coord - atoms[1].coord
            vec2 = atoms[2].coord - atoms[1].coord

            return (180/np.pi) * np.arccos(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)))

        # Base
        base = struct.base
        # Sítio
        site = struct.site

        if site[0:1] == "B" or site[0:1] == "D":

            # Dicionários de sítios
            sites = dops_data["graphine"]["sites_id"][base].copy()
            sites2 = graphine_data["sites_id"][base][0].copy()

            # Juntando dicionários
            sites.update(sites2)

            # Lista de nomes de sites
            sites_names = list(sites.keys())

            # Removendo sitíos não usados
            for rm_site in ["C3", "C4", "C5", "C6"]:
                sites_names.remove(rm_site)

            # Achando a posição do site do átomo dopante na lista de 
            # nomes de sites
            index = sites_names.index(site)

            # Listando átomos envolvidos na ligação
            atoms = []

            for i in [-1, 0, 1]:

                atoms.append(struct.atoms[sites[sites_names[index + i]] -1 ])

            #print(base, site, struct.dop_elem, bond_angle(atoms))
            # Retornar a diferença para o ângulo na estrutura base
            return bond_angle(atoms)

        else: 

            return None

    @staticmethod
    def dop_atom_torsion(struct: structure) -> float|None:
        """
        Se o dopante se encontra em um anel benzênico, calcula a 
        torção de um sequência de 4 átomos, onde o dopante é o segundo.
        Do contrário, retorna None.

        Essa métrica busca representar o quão pra fora do plano da 
        estrutura o átomo dopante se encontra.

        Parameters
        ----------

        struct : strucutre
            Estrutura de grafino sobre a qual se deseja fazer o cálculo.

        Returns
        -------
        
        float 
            Ângulo de torção, caso o dopantes esteja em um sítio A ou C.

        None
            Se o dopante não pertencer a um sítio A ou C.

        Raises
        ------
            Caso o dopantes não esteja em um sítio A ou C.
        """

        base = struct.base
        site = struct.site

        # Dicionários de sítios
        sites = dops_data["graphine"]["sites_id"][base].copy()
        sites2 = graphine_data["sites_id"][base].copy()

        atoms = []

        # Se sítio A1
        if site == "A1":

            for i in [2, 1, 0, 5]:
                atoms.append(struct.atoms[sites2[i]["A1"] - 1])

            return dihedral_angle(atoms)

        elif site in ["C1","C2","C3","C4"]:

            # Juntando dicionários
            sites.update(sites2[0])

            aux_sites = {

                "C1" : ["C5", "C6", "C1", "C2"],
                "C2" : ["C6", "C1", "C2", "C3"],
                "C3" : ["C1", "C2", "C3", "C4"],
                "C4" : ["C2", "C3", "C4", "C5"]
            }

            for aux_site in aux_sites[site]:

                atoms.append(struct.atoms[sites[aux_site] - 1])

            return dihedral_angle(atoms)

        else:
            return None

    ################################################################################################