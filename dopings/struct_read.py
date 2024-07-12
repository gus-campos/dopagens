import numpy as np

from dopings.atom_data import atom_data
from dopings.config import atoms_data, main_config


class struct_read:
    """
    Classe de métodos estáticos para a leitura de arquivos de 
    structures.

    Methods
    -------
    read_frame(struct: "structure") -> list[list[str|int]]

        Lê o último frame de um arquivo geo_end e o retorna como uma 
        lista 2D.

    read_energy(struct: "structure", extrapolated_0K=True) -> float

        Lê a energia total da estrutura em eV no último passo de 
        otimização e a retorna.

    read_bands(struct: "structure") -> tuple[float, float]

        Lê o arquivo band.out da estrutura e retorna a energia do HOMO e
        do LUMO.

    read_dipole(struct: "structure") -> np.ndarray[float]:

        Lê o vetor dipolo da estrutura, retornando-o como um array.
    """

    @staticmethod
    def read_frame(struct: "structure") -> list[list[str|int]]:
        """
        Lê o último frame de um arquivo geo_end e o retorna como uma 
        lista 2D.

        Parameters
        ----------

        struct : structure
            Estrutura para a qual se deseja ler os dados.
  
        Returns
        -------
        
        data : list[list[str|int]]
            Lista 2D onde o primeiro índice identifica uma linha,
            e o segundo índice identifica uma coluna de dados.
            Elemento, X, Y, Z, Carga

            Exemplo:

                ['H',  2.8054,  6.321, -0.357,  0.962]
                ['H',  4.9574,  5.091, -0.248,  0.913]
                ['C', -9.1570,  0.844, -0.195,  4.002]
                ['C', -9.1580, -0.353, -0.028,  4.214]
                ['C',  2.8572,  0.993,  0.223,  4.087]

                ...]

                data[3][0] -> 'C'
                data[1][2] -> 5.091

        Raises
        ------

        AssertionError
            The file has no data lines.

        AssertionError
            Unexpected number of data in geometry file.
        """

        # Usar caminho do geo_end
        coord_dir = struct.dir / "geo_end.xyz"

        # Se o arquivo não existir, retornar None
        if not coord_dir.is_file():

            return None

        # Lendo arquivo
        with open(coord_dir) as file:

            txt = file.read()
            lines = txt.splitlines()

        # Número de linhas
        number_of_lines = int(len(lines))

        # Report se vazio
        if number_of_lines < 3:

            raise AssertionError("The file has no data lines.")

        # Lendo da primeira linha o tamanho da estrutura
        struct_size = int(lines[0].strip())
        page_size = int(struct_size + 2)
        number_of_pages = int(number_of_lines/page_size)

        # Lista para guardar linhas de inicio do frame
        frames_starting_lines = []

        # Encontando linhas
        for i in range(number_of_pages):
                
            frames_starting_lines.append(i*page_size)

        # Incluindo a posição de onde seria um próximo frame após o último
        frames_starting_lines.append(number_of_lines)

        # Pegando o último frame
        frame = lines[frames_starting_lines[-2] : frames_starting_lines[-1]]

        # Removendo cabeçalho
        frame = frame[2:]

        # Quantidade de dados
        number_of_data = len(frame[0].split())

        # Transformando cada linha em uma lista de colunas
        for i in range(len(frame)):
            # Separando colunas
            frame[i] = frame[i].split()      

            # Convertendo strings numéricas para floats   
            # Se tiver 4 valores, não tem carga    
            if number_of_data == 4:
                frame[i] = [ frame[i][0], float(frame[i][1]), float(frame[i][2]), float(frame[i][3]) ]

            # Se tiver 5 valores, tem carga
            elif number_of_data == 5:
                frame[i] = [ frame[i][0], float(frame[i][1]), float(frame[i][2]), float(frame[i][3]), float(frame[i][4]) ]

            else:
                raise AssertionError("Unexpected number of data in geometry file.")

        # Retornar
        return frame
    
    @staticmethod
    def read_energy(struct: "structure", extrapolated_0K=True) -> float:
        """
        Retorna a energia total da estrutura (eV) no último passo de 
        otimização.

        Por padrão, já lê a energia extrapolada para o 0K.

        Parameters
        ----------

        struct : structure
            Estrutura para a qual se deseja ler os dados.

        extrapolated_0K : bool, default = True
            Se deve ser a energia extrapolada para o 0K, ou só a energia
            total.

        Returns
        -------

        float
            Energia total da estrutura (eV) no último passo de 
            otimização.

        Raises
        ------

        AssertionError
            Energies not found in file.
        """

        # Diretório
        dir = struct.dir / main_config["output_name"] 

        # Se não encontrar arquivo, retornar None
        if not dir.is_file():
            return None

        # No arquivo, varrer de trás pra frente em busca de ocorrência do termo.
        # Retornar o primeiro que encontrar
        with open(dir) as file:

            # Se for pedido para calcular energia extrapolada pro 0K
            if extrapolated_0K:
                    
                # Pra cada linha do arquivo, na ordem inversa
                for linha in list(file)[::-1]:

                    # Se for encontrado tal trecho
                    if linha.find("Extrapolated to 0K") != -1:

                        # Da linha, retornar apenas o valor relevante, em float
                        return float(linha.split()[5])
                    
            else:

                # Se for encontrado tal trecho
                if linha.find("Total Energy") != -1:
                        
                    # Pra cada linha do arquivo, na ordem inversa
                    for linha in list(file)[::-1]:

                        # Da linha, retornar apenas o valor relevante, em float
                        return float(linha.split()[4])

            # Se correr todo o arquivo sem encontrar o termo, declarar erro
            raise AssertionError("Energies not found in file.")
    
    @staticmethod
    def read_bands(struct: "structure") -> tuple[float, float]:
        """
        Lê o arquivo band.out da estrutura e retorna a energia do HOMO e
        do LUMO.

        Necessita que o número de elétrons na camada de valência do 
        elemento esteja definido em "atoms_data".

        Parameters
        ----------

        struct : structure
            Estrutura para a qual se deseja ler os dados.

        Returns
        -------

        tuple[float, float]
            Energia do homo, e energia do lumo.

        Raises
        ------

        KeyError
            valency of element [element] not found.
        """

        # Importando bibliotecas necessárias
        from pathlib import Path

        # Convertendo endereço para Path
        bands_dir = struct.dir / "band.out"

        # Verificando se exite, se não existir, retornar None
        if not bands_dir.is_file():
            return None, None

        ### Listando dados do band.out
        data = []

        # Salvando linhas
        with open(bands_dir) as file:

            # Listando linhas
            for linha in file:

                # Separar os itens da linha numa lista
                data.append(linha.split())

        # Removendo cabeçalho e tail
        data = data[1:-1]

        # Convertendo strings numéricas para floats
        for i in range(len(data)):

            # Lista com o número da banda, sua energia, e sua ocupação
            data[i] = [int(data[i][0]), float(data[i][1]), float(data[i][2])] 

        # Calculando quantidade de elétrons de valência
        val_elect = 0

        # Somando os elétrons de valência de cada átomo
        for atom in struct.atoms:
            try:
                val_elect += atoms_data["valency"][atom.elem]
            except:
                raise KeyError(f"valency of element {atom.elem} not found.")

        def ceil_to_pair(n: int) -> int:

            """
            Encontra o primeiro inteiro par que vem depois de "n", caso 
            ele já não seja par.
            """

            # Se for par retornar divisão inteira
            if n % 2 == 0:
                return n // 2
            
            # Se não for par, retornar divisão inteira acrescida de 1
            else:
                return (n // 2) + 1

        # Calculando a última banda ocupada
        homo_band = ceil_to_pair(val_elect)
        lumo_band = homo_band + 1

        # Encontrando energia das bandas
        homo_energy = data[homo_band-1][1]
        lumo_energy = data[lumo_band-1][1]

        return homo_energy, lumo_energy
    
    @staticmethod
    def read_dipole(struct: "structure") -> np.ndarray[float]:
        """
        Lê o vetor dipolo da estrutura, retornando-o como um array.

        Parameters
        ----------

        struct : structure
            Estrutura para a qual se deseja ler os dados.

        Returns
        -------

        np.ndarray[float]
            Array do numpy com 3 floats que corresponde ao vetor dipolo
            da estrutura.

        None
            Se o arquivo não existir.

        Raises
        ------

        AssertionError
            detailed.out file exists, but dipole moment was not found.
        """

        import numpy as np

        # Importando bibliotecas necessárias
        from pathlib import Path

        # Convertendo endereço para Path
        detailed_dir = struct.dir / "detailed.out"

        # Verificando se exite, se não existir, retornar None
        if not detailed_dir.is_file():
            return None

        # Salvando linhas
        with open(detailed_dir) as file:

            # Listando linhas
            for linha in file:

                # Separar os itens da linha numa lista
                if "au" in linha.split():
                    
                    words = linha.split()

                    return np.array([words[2], words[3], words[4]])
                
            # Se não foi encontrado, retornar None
            return None
    