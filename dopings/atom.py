
import numpy as np

"Este módulo contém a classe Atom"

class Atom:

    """
    Representa um átomo, a partir do seu símbolo de elemento, posição e 
    carga.

    Attributes
    ----------
    elem : str
        Símbolo do elemento do átomo (maiúsculas e minúsculas são
        diferenciadas).

    coord : np.ndarray[float, float, float]
        Aray de tamanho 3, para armazenar na ordem coordenadas X, Y, Z 
        do átomo.

    charge : float
        Carga do átomo, em população de elétrons.

    Methods
    -------
    
    def __init__(self, elem: str=None, coord: list[float]=None, charge: float = None) -> None

        Inicializa objeto Atom.

    def dist_to(self, atom2: Atom) -> float

        Calcula distância de si, até outro Atom.        
    """

    def __init__(self, elem: str=None, coord: list[float]=None, 
                 charge: float = None) -> None: 

        """
        Inicializa objeto Atom.

        Parameters
        ----------
        elem : str, default = None
            Símbolo do elemento do átomo.

        coord : list[float], default = None
            Lista de 3 itens com as coordenadas X, Y, Z do átomo.

        charge : float, default = None
            Carga do átomo, em população de elétrons.

        Raises
        ------
        TypeError
            "Elem type must be string."

        TypeError
            "Coord must be list of 3 floats."

        TypeError
            "Charge must be float."
        """

        # Declarando atributos
        self.elem = None
        self.coord = None
        self.charge = None

        # Checando o tipo da variável elem e atribuindo
        if elem is not None:
            
            if isinstance(elem, str):
                self.elem = elem
            else:
                raise TypeError("Elem type must be string.")
        
        # Checando tipo e dimensão da variável coord e atribuindo
        if coord is not None:

            if (isinstance(coord, list) and len(coord) == 3
                and all(isinstance(coord, float) for coord in coord)):
                    
                # Transformar em array e salvar
                self.coord = np.array(coord, dtype=float)

            else:
                raise TypeError("coord must be list of 3 floats.")

        # Checando o tipo da variável charge e atribuindo
        if charge is not None:
            
            if isinstance(charge, float):
                self.charge = charge
            else:
                raise TypeError("Charge must be float.")  

    ############# Distância até outro átomo

    def dist_to(self, atom2: "Atom") -> float:
    
        """
        Calcula distância de si, até outro átomo.
        
        Parameters
        ----------

        atom2 : Atom
            Outro objeto Atom.

        Returns
        -------

        float
            Distância até o outro átomo.

        Raises
        ------

        TypeError
            Both atoms must have defined positions.

        TypeError
            atom2 must be an Atom object.
        """

        # Se alguma coordenada for None
        if (self.coord is None) or (atom2.coord is None):
            raise TypeError("Both atoms must have defined positions.")

        # Se for passado um objeto Atom, retornar a distância
        if isinstance(atom2, Atom):
            return np.linalg.norm(self.coord - atom2.coord)
        else:
             raise TypeError("atom2 must be an Atom object.")
