
from dopings.doping_set import DopingSet
from dopings.config import dops_data

# Classe faixa
class faixa:

    def __init__(self, id, min, max):

        self.id = id
        self.min = min
        self.max = max
        self.count = 0

    def add(self):

        self.count += 1

# Lista de faixas encontradas
faixas_graphine = {

    "Al" : [ 

        faixa("alpha", 1.63, 1.68),
        faixa("beta", 1.76, 1.84),
        faixa("gamma", 1.89, 2.04)
    ],

    "B" : [ 

        faixa("alpha", 1.29, 1.34),
        faixa("beta", 1.45, 1.50)
    ], 

    "Mg" : [

        faixa("alpha", 2.00, 2.07),
        faixa("beta", 2.08, 2.19),
        faixa("gamma", 2.20, 2.25)
    ],
    
    "N" : [ 

        faixa("alpha", 1.18, 1.20),
        faixa("beta", 1.19, 1.27),
        faixa("gamma", 1.30, 1.37)
    ],

    "Na" : [ 

        faixa("alpha", 2.49, 3.91),

    ],


    "Si" : [ 

        faixa("alpha", 1.55, 1.61),
        faixa("beta", 1.67, 1.78),
        faixa("gamma", 1.80, 1.84)
    ],

    "Li" : [ 

        faixa("alpha", 1.89, 2.09)
    ],

    "P" : [ 

        faixa("alpha", 1.49, 1.58),
        faixa("beta", 1.60, 1.74),
        faixa("gamma", 1.77, 1.89)
    ],  

    "O" : [

        faixa("alpha", 1.31, 1.42),
        faixa("beta", 1.43, 1.48),
        faixa("gamma", 1.48, 1.63)
    ],

    "Ti" : [

        faixa("alpha", 1.80, 1.88),                
        faixa("beta", 2.03, 2.28),
        faixa("gamma", 2.46, 2.56)

    ],

    "Zn" : [

        faixa("alpha", 0, 2.35),                
        #faixa("beta", 2.03, 2.28),
        #faixa("gamma", 2.46, 2.56)

    ]
}

faixas_graphene = {

    "Al" : [ 

        faixa("alpha", 1.63, 1.68),
        faixa("beta", 1.76, 1.84),
        faixa("gamma", 1.89, 2.04)
    ],

    "B" : [ 

        faixa("alpha", 1.29, 1.34),
        faixa("beta", 1.45, 1.50)
    ], 

    "Mg" : [

        faixa("alpha", 2.00, 2.07),
        faixa("beta", 2.08, 2.19),
        faixa("gamma", 2.20, 2.25)
    ],
    
    "N" : [ 

        faixa("alpha", 1.18, 1.20),
        faixa("beta", 1.19, 1.27),
        faixa("gamma", 1.30, 1.37)
    ],

    "Na" : [ 

        faixa("alpha", 2.49, 3.91),

    ],


    "Si" : [ 

        faixa("alpha", 1.55, 1.61),
        faixa("beta", 1.67, 1.78),
        faixa("gamma", 1.80, 1.84)
    ],

    "Li" : [ 

        faixa("alpha", 1.89, 2.09)
    ],

    "P" : [ 

        faixa("alpha", 1.49, 1.58),
        faixa("beta", 1.60, 1.74),
        faixa("gamma", 1.77, 1.89)
    ],  

    "O" : [

        faixa("alpha", 1.31, 1.42),
        faixa("beta", 1.43, 1.48),
        faixa("gamma", 1.48, 1.63)
    ],

    "Ti" : [

        faixa("alpha", 1.80, 1.88),                
        faixa("beta", 2.03, 2.28),
        faixa("gamma", 2.46, 2.56)

    ],

    "Zn" : [

        faixa("alpha", 0, 2.35),                
        #faixa("beta", 2.03, 2.28),
        #faixa("gamma", 2.46, 2.56)

    ]
}


# Categoriza um valor de comprimento
def categorize(value, faixas, count=False):

    for faixa in faixas:
        if faixa.min <= value and value <= faixa.max:
            
            # Incrementa contadores
            if count:
                faixa.add()

            return faixa.id

    # Se não atener a nenhum outro critério
    return ""

def bonded_to_ring(base, site):

    # Se for o primeiro sítio do 
    if (site[:1] == "B" or site[:1] == "D") and (site[1:] == "1"):
        return True
    
    # Se for B e o número for o dobro do número da base
    elif (site[:1] == "B") and (site[1:] == str(2*int(base[1:]))):
        return True
    
    else:
        return False

def detect_values_jumps(all_bonds, n=5):

    # Copiando lista de ligações
    all_bonds = sorted(all_bonds.copy())

    # Lista de saltos, com valor inicial para que tenham mesmo tamanho
    jumps = [-1]

    # Para cada índice da lista de ligações, exceto o último
    for i in range(1, len(all_bonds)):
        
        # Salvar na lista de saltos a diferença entre o índice próximo e o atual
        jumps.append(all_bonds[i] - all_bonds[i-1])

    # Ordenar de forma paralela as listas
    jumps, all_bonds = zip(*sorted(zip(jumps, all_bonds), reverse=True))

    # Cabeçalho
    print("==== Maiores saltos\n")

    # Imprimir primeiros 5 itens das listas de forma paralela
    for bond, jump in zip(all_bonds[:n], jumps[:n]):

        print(f"{bond:.3f} {jump:.3f}")


for dop_elem in ["Li"]: #["Al", "B", "Li", "N", "O", "P", "Si", "Ti"]:

        # Lê um conjunto de dopagens
        set = DopSet(dop_elem = dop_elem, dops_info=dops_data["graphene"], mode="read", param="matsci")

        # FAIXAS
        faixas = faixas_graphene

        # Para armazenar lista de ligações
        bonds = []

        # total de ligações
        total_of_bonds = 0

        # Contadores
        all_bonds = []


        # Para cada estrutura
        for struct in set.structs:

            # Base
            base = struct.base
            # Sítio
            site = struct.site


            ##############################################
            #            Seleção de sítios
            ##############################################

            # Todos
            if True:

            # Todos, exceto...
            #if not (base == "g3" and site == "D3"):

            # Sítios ligados ao anel
            #if bonded_to_ring(base, site):

            # Sítios B e D não ligados ao anel
            #if (site[:1] == "B" or site[:1] == "D") and not bonded_to_ring(base, site):            

            # Ligados a 3 carbonos
            #if site in ["A1", "C1", "C2"]:
            
            # Ligados a 2 carbonos e 1 hidrogênio
            #if site in ["C3", "C4"]:

            ###############################################

                # Elemento buscado
                elem = "H"

                # Calcular comprimentos de ligação com elem
                struct_bonds = struct.shortest_distances(dop_elem, elem, 6)      # Carbonos

                # Imprimir diretório
                print(f"==== {base}-{site}\n")

                if elem == "C":

                    # Número de ligações de acordo com o sítio
                    n = 3 if (site in ["A1", "C1", "C2"]) else 2

                if elem == "H":

                    # Número de ligações de acordo com o sítio
                    n = 1 if (site in ["C3", "C4"]) else 0

                if elem == "C" and dop_elem == "Ti":

                    n = 6

                # N primeiros comprimentos, adicionar à lista menor, e imprimir com identificação
                for bond in sorted(struct_bonds[:n]):

                    # Imprimir comprimento
                    print(f"{bond:.2f}", categorize(bond, faixas[dop_elem]))
                    
                    # Adicionar à lista de comprimentos
                    bonds.append(bond)
                    
                    # Categorizar, fazendo contagem
                    categorize(bond, faixas[dop_elem], count=True)

                # Nova linha
                print("")

                # Para todos, adicionar à lista
                for bond in struct_bonds:
                    # Adicionar à lista 
                    all_bonds.append(bond)

        print("=== Contadores\n")

        # Imprimir contadores
        for elem_faixa in faixas[dop_elem]:
            print(f"{elem_faixa.id}: {elem_faixa.count}")

        # Separador
        print("\n"+"="*30)

        # Imprimindo estatísticas
        from statistics import mean
        print("\n=== Stats\n")
        print("Média:", f"{mean(bonds):.2f}")
        print("Min:", f"{min(bonds):.2f}")
        print("Max:", f"{max(bonds):.2f}")

        # Nova linha
        print("")

        # Imprimir todos os comprimentos da lista geral
        if True:

            # Separador
            print("="*30)

            # Imprimir cada ligação, de forma ordenada
            for bond in sorted(all_bonds):
                print(f"{bond:.3f}", categorize(bond, faixas[dop_elem]))
        
            # Separador
            print("="*30)

            # Chamando detector de salto
            detect_values_jumps(bonds, n=8)

            