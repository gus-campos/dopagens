from __future__ import annotations

import numpy as np
from pathlib import Path

from DOPINGS_PACKAGE.atom_data import atom_data
from DOPINGS_PACKAGE.struct_read import struct_read

from DOPINGS_PACKAGE.global_vars import atoms_data, dirs_data, main_config

###############################################################################

class structure:

    """
    Representa uma estrutura formada por um conjunto de átomos.

    Attributes
    ----------
    atoms : np.ndarray[atom_data]
        Array de atom objects na estrutura.
    
    size : int
        Número de átomos na estrutura.
        Se atualiza sozinho.
    
    dir : Path, default: None
        Diretório de otimização da molécula.

    material : str
        Material dopado.

    param : str
        Set of parameters used in the optimization.
        Used only in formation energy.
    
    dop_elem : str
        Elemento da dopagem.
        
    base : str
        Base da estrutura.

    site : str
        Sítio de dopagem da estrutura.

    param : str
        Conjunto de parâmetros usado na otimização. Usado apenas para
        cálculo da energia de formação.
    
    total_energy : float
        Energia total da estrutura, lida a partir do arquivo

    homo : float
        Energia do homo da estrutura, lida a partir do arquivo

    lumo : float
        Energia do lumo da estrutura, lida a partir do arquivo

    dipole : numpy.ndarray[float]
        Vetor momento dipolo da estrutura.

    Methods
    -------

    __init__(self, atoms : list[atom_data]=None, 
             dir : str|Path=None, read_from_dir : bool=False, 
             material: str=None, base: str=None, dop_elem: str=None,
             site: str=None, param : str=None) -> None

        Inicializa o objeto structure.
    
    append(self, atoms: atom_data|list[atom_data]) -> None

        Anexa um atom_data object, ou uma lista de atom_data objects 
        à lista de átomos da estrutura.
    
    redo_files_for_resume(self) -> None

        No caso de uma otimização ser encerrada antes de concluída, este
        método refaz os arquivos de otimização de forma que parta do 
        último frame escrito, guardando o progresso até então em uma 
        pasta no diretório, que recebe o nome definido em "main_config" 
        (resume_dir_name), que por padrão é "before_restart".

    data_from_file(self) -> None

        A partir do diretório da estrutura, lê os arquivos de dados e os 
        carrega como propriedades do objeto structure.
    
    opt(self, overwrite: bool=False, verbose: bool=True, 
        skip_hard_to_conv_SCC: bool=True, resume_unfineshed: bool=True
        ) -> None

        Executa uma otimização no diretório da estrutura, espera 
        a otimização terminar e relata o status ao final.
    
    report(self, return_SCC: bool=False, return_conv: bool=False, 
           verbose: bool=True, 
           return_written: bool=False) -> None|bool|list[int]

        Reporta o estado da otimização da estrutura a partir do
        arquivo de no seu diretório, 
    
    dope_out(self, dop_elem: str, id: int, 
             kwargs: dict=None) -> structure

        Gera uma versão (dopada) de uma estrutura, onde o átomo 
        especificado é substituído por um de outro elemento.
    
    formation_energy(self) -> tuple[float, float]

        Retorna a energia de formação da estrutura, além de sua média 
        por átomo.
    
    homo_lumo(self) -> tuple[float, float, float]

        Retorna as energias de HOMO, LUMO, e o módulo da diferença 
        HOMO-LUMO de uma estrutura a partir de seu arquivo band.out.
    
    shortest_distances(self, elem1: str, elem2: str, 
                       n: int) -> list[float]

        Retorna as "n" menores distâncias entre átomos de dados 
        elementos, na estrutura.
    
    frame(self, frame_out: str|Path) -> None

        Gera e escreve o arquivo xyz do último frame da estrutura.
        Cria diretórios necessários caso não existam.
    
    output(self, output_dest: str|Path) -> None

        Copia o output de otimização da estrutura para o endereço 
        passado.

    copy(self) -> structure
        
        Retorna uma estrutura que contém os mesmo átomos que a original.
    """
     
    def __init__(self, atoms : list[atom_data]=None, dir : str|Path=None, 
                 read_from_dir : bool=False, material: str=None, 
                 base: str=None, dop_elem: str=None, site: str=None, 
                 param : str=None) -> None:

        """
        Inicializa o objeto structure.

        Parameters
        ----------

        atoms : list[atom_data], default = None
            Lista de átomos da estrutura criada.

        dir : str or Path, default: None
            Diretório de otimização da molécula.
            Obrigatório se read_from_dir for True.

        material : str
            Material dopado.

        param : str
            Set of parameters used in the optimization.
            Used only in formation energy.
        
        dop_elem : str
            Elemento da dopagem.
            
        base : str
            Base da estrutura.

        site : str
            Sítio de dopagem da estrutura.

        read_from_dir : bool, default=False
            Se dados devem ser lidos do diretório ou não.
            Se True, dir deve ser passada.

        Raises
        ------

        TypeError
            "dir must be string or Path."

        ValueError
            "If dir it's not given, read_from_dir must be False."
        """

        # Diretório de otimização da molécula
        self.dir = None
        self.size = 0
        self.atoms = np.array([], dtype = atom_data)
        
        # -- Dop info --
        self.param = None
        self.material = None
        self.base = None
        self.dop_elem = None
        self.site = None
        
        # -- Proc data --
        self.homo = None
        self.lumo = None
        self.dipole = None
        self.total_energy = None

        ####### Checando tipos e guardando variáveis

        if dir is not None:
            if isinstance(dir, str|Path):
                self.dir = Path(dir).expanduser()
            else:
                raise TypeError("dir must be string or Path.")

        # Se pediu pra ler dos arquivos
        if read_from_dir:
            if self.dir is not None:
                self.data_from_file()
            else:
                raise ValueError("If read_from_dir True, dir must be given.")
        
        # Decidindo se será anexado
        if ((dir is None) or (not read_from_dir)) and (atoms is not None):

            # Adicionar átomos que foram passados
            # Internamente o tamanho é atualizado e o type é verificado
            self.append(atoms)
        
        # Verificando propriedade "param"
        if param is not None and isinstance(param, str):
            self.param = param

        if material is not None and isinstance(material, str):
            self.material = material

        if base is not None and isinstance(base, str):
            self.base = base

        if dop_elem is not None and isinstance(dop_elem, str):
            self.dop_elem = dop_elem

        if site is not None and isinstance(site, str):
            self.site = site

    ############# Adicionar novos átomos à estrutura - Sem dados obrigatórios
             
    def append(self, atoms: atom_data|list[atom_data]) -> None:

        """
        Anexa um atom_data object, ou uma lista de atom_data objects 
        à lista de átomos da estrutura.

        Parameters
        ----------

        atoms: list[atom_data] | atom_data
            Lista de átomos a ser adicionada.

        Raises
        ------

        TypeError
            "atoms must be a non-empty list of atoms objects"
        """

        # Se for atom_data object, tranformar em lista unitária
        if isinstance(atoms, atom_data):
            atoms = [atoms]
        
        # Transformando em array
        if isinstance(atoms, list):
            atoms = np.array(atoms)

        # Se for array, composto por atoms objects
        if (isinstance(atoms, np.ndarray) and atoms.size > 0 
            and all(isinstance(atom, atom_data) for atom in atoms)):
            
            # Anexar a lista passada ao array de atoms.
            self.atoms = np.append(self.atoms, atoms)
            # Acrescentar ao tamanho, a quantidade adicionada
            self.size = self.size + len(atoms)
        
        else:
            raise TypeError(("atoms must be a non-empty list of atom_data"
            + "objects, or atom_data."))

    ############# Otimização e relatório

    def redo_files_for_resume(self) -> None:

        """
        No caso de uma otimização ser encerrada antes de concluída, este
        método refaz os arquivos de otimização de forma que parta do 
        último frame escrito, guardando o progresso até então em uma 
        pasta no diretório, que recebe o nome definido em "main_config" 
        (resume_dir_name), que por padrão é "before_restart".
        """

        import shutil

        # Diretório de backup
        backup_file = self.dir / main_config["resume_dir_name"]

        # Criando diretório
        Path.mkdir(backup_file, exist_ok=True)

        # Para cada arquivo do diretório
        for file in Path(self.dir).glob("*"):
            # Se for um arquivo
            if file.is_file():
                # Se for um arquivo dinâmico
                if file.name == "geo_end.xyz" or file.name == "output":
                    
                    # Abrir sua versão no diretório "before_restart"
                    with (backup_file / file.name).open('a') as copy_file:
                        # Abrir arquivo a ser copiado
                        with file.open() as src_file:

                            # Adicionar arquivo a ser copiado no final 
                            # do arquivo backup
                            copy_file.write(src_file.read())

                # Se for qualquer outro arquivo, apenas copiar
                else:
                    # Copiar pra pasta backup
                    shutil.copy(file, backup_file / file.name)

        # Lendo dados da structure (lê os dados mesmo que esteja em look mode)
        struct = structure(dir=self.dir, read_from_dir=True)

        # Gerando novos arquivos de otimização, e sobrescrevendo antigos 
        from dops_set import dops_set
        
        # Tentar gerar arquivos, com ReadInitialCharges
        try:
            dops_set.gen_files(self, read_initial_charges=True, redo=True)
        
        # Se não for possível, restaurar estado anterior
        except:

            # Para cada arquivo do diretório
            for file in Path(backup_file).glob("*"):

                # Copiar de volta os arquivos guardados em backup (sobrescreve)
                shutil.copy(file, backup_file.parent)

            # Remover pasta de backup
            shutil.rmtree(backup_file)

            # Anunciando que dados não foram alterados
            print("OPT FILES RESTORED BECAUSE REDO COULD NOT BE DONE")

            # Relançando erro do último bloco
            raise AssertionError("It was not possible to gen files for restart.")

        # Copiando o último frame do geo_end.xyz como inp.xyz 
        # (sobrescreve o inp.xyz gerado)
        struct.frame(struct.dir / "inp.xyz")

        # Removendo todos os arquivos, exceto os que serão usados
        for file in Path(self.dir).glob("*"):
            # Se for um aquivo
            if file.is_file():
                # Se não for charges.dat
                if file.name not in ["inp.xyz", "charges.dat", "dftb_in.hsd"]:

                    # Remover
                    Path.unlink(file)

    def data_from_file(self) -> None:

        """
        A partir do diretório da estrutura, lê os arquivos de dados e os 
        carrega como propriedades do objeto structure.
        """

        # Lendo coordenadas
        frame_data = struct_read.read_frame(self)

        # Quantidade de dados
        if frame_data is not None:
            
            number_of_data = len(frame_data[0])

            # Para cada átomo da estrutura
            for i in range(len(frame_data)):
                
                # Juntando argumentos em um dict
                temp_dict = {

                    "elem"      :   frame_data[i][0] if frame_data is not None else None,

                    "coord"     :   [frame_data[i][1], 
                                    frame_data[i][2], 
                                    frame_data[i][3]] if frame_data is not None else None,

                                                                                # e se houverem 5 dados por átomo
                    "charge"    :   frame_data[i][4] if ((frame_data is not None and number_of_data==5)) else None
                }

                # Passando pares de chaves e valores do dicionário como 
                # argumentos, usando "**"
                self.append(atom_data(**temp_dict))  

        # Lendo energia
        self.total_energy = struct_read.read_energy(self)

        # Lendo homo e lumo
        self.homo, self.lumo = struct_read.read_bands(self)

        # Lendo momento dipolo
        self.dipole = struct_read.read_dipole(self)

    def opt(self, overwrite: bool=False, verbose: bool=True, 
            skip_hard_to_conv_SCC: bool=True, resume_unfineshed: bool=True
            ) -> None:

        """
        Executa uma otimização no diretório da estrutura, espera 
        a otimização terminar e relata o status ao final.
        
        Parameters
        ----------

        overwrite
        verbose
        skip_hard_to_conv_SCC
        resume_unfineshed

        overwrite : bool, default=False
            Se a otimização deve ser executada mesmo que já exista no 
            diretório output de otimização anterior. Se True, a 
            otimização será iniciada e os arquivos anteriores serão 
            sobrescritos.

        verbose : bool, default=True
            Se, como informação extra, os passos onde não houve 
            convergência do SCC devem ser impressos também.

        skip_hard_to_conv_SCC : bool, default=True
            Se otimizações com demasiados passos SCC não convergidos 
            devem ser puladas.

            Regra para determinar uma otimizão com demasiados passos:
            
            Ter mais de 30 passos não convergidos; ou ter mais de 20
            passos não convergidos, antes de chegar no passo 30 da 
            otimização.

        resume_unfineshed : bool, default=True
            Se True, ao realizar a otimização que foi anteriormente 
            interrompida, gera os arquivos necessários para recomeçar do
            último passo escrito e inicia tal otimização.

            Se False, recomeça do primeiro frame escrito.
        """

        import subprocess
        from time import sleep

        def run_opt(opt_dir):

            # Rodar otimização no diretório e guardar pid
            pid = subprocess.run(f"nohup {bin_dir} >& {main_config['output_name']} & echo $!", 
                                 shell=True, executable='/bin/bash', 
                                 cwd=opt_dir, capture_output=True
                                 ).stdout.decode("ascii").strip()
            
            return pid
        
        def kill_pid(pid):

            print(f"try to kill {pid}")

            subprocess.run(f"kill {pid}", shell=True, 
                                    executable='/bin/bash', capture_output=True)
            
        def pid_running(pid):

            running = subprocess.run(f"if $(ps -p {pid} > /dev/null); then echo True; else echo False; fi", shell=True, 
                                     executable='/bin/bash', capture_output=True
                                     ).stdout.decode("ascii").strip() 
            
            return running == "True"

        # Transformando endereço em Path      
        # Usar diretório da estrutura como o de otimização
        opt_dir = Path(self.dir).expanduser()

        # Endereço do binário
        bin_dir = dirs_data["bin_dir"][main_config["maquina"]]

        # Se o arquivo output existir
        if Path.exists(opt_dir / main_config["output_name"]):

            # Abrir o arquivo e verificar nele o status
            with open(opt_dir / main_config["output_name"]) as file:

                # Buscar no arquivo se convergiu
                converged = file.read().find("Geometry converged") != -1

        # Se o arquivo não existir, considerar como "não convergido"
        else:
            converged = False

        # Guardar quantos passos foram escritos
        written = self.report(return_written=True)

        # A exigência de 1 passo escrito é para evitar que ele refaça os 
        # arquivos em caso de já estarem refeitos
        if ((resume_unfineshed) and (not overwrite) and (not converged) 
            and (written is not None and len(written) > 1)):

            print("Refazendo diretórios...")
            # Faz backup e refaz os arquivos para continuar de onde parou
            self.redo_files_for_resume()

        # Se ainda não convergiu ou for pedido para sobrescrever
        if (not converged) or (overwrite):
            pid = run_opt(opt_dir)
    
            # Em loop infinito
            while True:

                # Se estiver configurado para pular estruturas com grande 
                # dificuldade de convergir o SCC
                if skip_hard_to_conv_SCC:

                    # Encontrar lista de SCC não convergido
                    nConv_SCC = self.report(return_SCC=True)

                    # Se tiver mais que 30 SCC não convergido,
                    # ou chegar a 20 nconv antes do passo 30
                    if (len(nConv_SCC) > 30 
                        or (len(nConv_SCC) > 0 and nConv_SCC[-1] < 30 
                        and len(nConv_SCC) > 20)):

                        # Encerrar dftb+
                        kill_pid(pid)
                        
                # Se o dftb+ não estiver rodando mais, sair do loop
                if not pid_running(pid):

                    break
                
                sleep(5)

        # Reportar status a partir do diretório da estrutura 
        self.report(verbose=verbose)

    def report(self, return_SCC: bool=False, return_conv: bool=False, 
               verbose: bool=True, 
               return_written: bool=False) -> None|bool|list[int]:

        """
        Reporta o estado da otimização da estrutura a partir do
        arquivo de no seu diretório, 
        
        Relata:

        - Ausência de arquivo de output.
        
        - Se convergiu e com quantos passos convergiu.
        - Se houve divergência de SCC, em quais passos isso
        ocorreu e que porcentagem isso representa do total de passos.
        - Se houve divergência de SCC no último passo de otimização.
        
        - Se não convergiu.
        - Se houve divergência de SCC, em quais passos isso
        ocorreu e que porcentagem isso representa do total de passos.
        
        - Se o estado está indefinido.
        - Se nenhum passo foi escrito no arquivo.
        - Qual o último passo escrito no arquivo.
        - Se houve divergência de SCC, em quais passos isso
        ocorreu e que porcentagem isso representa do total de passos.


        Parameters
        ----------

        return_SCC : bool, default=False
            Se True, retorna uma lista dos passos não convergidos do 
            SCC.
    
        return_conv : bool, default=False
            Se True, retorna se a otimização convergiu ou não.

        verbose : bool, default=True
            Se, como informação extra, os passos onde não houve 
            convergência do SCC devem ser impressos também.

        return_written : bool, default=False
            Se True, retorna uma lista dos passos escritos no arquivo.
        
        Returns
        -------

        list[int]
            No caso de return_SCC ou return_written for True.

        bool          
            No caso de return_conv for True.

        None
            No caso de nenhuma opção de retorno ser True.
        """

        import os

        # Diretório de otimização
        opt_dir = self.dir
        output_file = opt_dir / main_config["output_name"]

        # Encontra todas as ocorrências de uma substring
        def find_all(data: str, substring: str) -> list[int]:
            """
            Função que encontra a posição de todas as ocorrências de uma 
            substring em uma string maior (como de um documento).

            Returns
            -------

            list[int]
                Lista de posições onde a substring foi encontrada.
            """

            # Inicializando lista de passos
            steps = []
            # Inicializando primeiro passo condiderado como -1
            step = -1

            while True:

                # Continuar na posição seguinte à última encontrada
                step = data.find(substring, step+1)

                # Se encontrar próxima ocorrência
                if step != -1:
                    steps.append(step)

                # Do contrário, sair do loop
                else: 
                    break

            return steps

        # Função para encontrar em que passo se encontra uma posição 
        # de caracter
        def indexes_to_steps(steps_indexes: list[int], 
                             searched_indexes: list[int]) -> list[int]:
            """
            Função que identifica em que passo da otimização se encontra
            um índice qualquer buscado.

            No caso, recebe uma lista de indexes, e devolve uma lista 
            com os passos onde cada um dos indexes se encontra.

            Parameters
            ----------

            steps_indexes : list[int]
                Lista dos índeces de início dos passos.

            searched_indexes : list[int]
                Lista de indexes procurados dentro dos passos.

            Returns
            -------

            list[int]
                Lista de mesmo tamanho que searched_indexes, onde cada 
                item é o passo onde foi encontrado o index buscado
                correspondente.

            Raises
            ------
            AssertionError
                index not found in file.
            """

            # Lista para acumular passos
            found_steps = []

            # Número de passos da otimização
            steps_num = len(steps_indexes)
            # Atribundo a última posição do arquivo como sendo um passo posterior ao último
            steps_indexes = steps_indexes + [ len(data)-1 ]

            # Para cada index buscado
            for index in searched_indexes:
            
                # Avaliar cada passo:
                for step in range(0, steps_num):

                    # Se o index buscado for maior que o do passo atual,
                    # e menor que o do passo seguinte
                    if (index > steps_indexes[step] 
                        and index < steps_indexes[step+1]):

                        # Anexar o passo avaliado como sendo o de tal index
                        found_steps.append(step)
                        break
                    
                    # Se não achou, e já processou o último passo
                    elif step == steps_num-1:
                            raise AssertionError("index not found in file.")

            # Retornar dados
            return found_steps

        # Imprime status
        def print_status(exists: bool=None, 
                         converged: bool=None, 
                         not_converged: bool=None,
                         steps_written: list=None,
                         not_conv_SCC: list=None) -> None:
            
            """
            Imprime o estado dos arquivos. De acordo com os parâmetros
            passados.

            Parameters
            ----------

            exists : bool, default=None
                Se o arquivo de output existe.

            converged : bool, default=None
                Se foi relatada convergência.

            not_converged : bool, default=None
                Se foi relatada a não convergência.
            
            steps_written : list[int], default=None
                Passos escritos no arquivo de output.
            
            not_conv_SCC : list[int], default=None
                Passos onde o SCC não convergiu. 
            """

            messages = []

            # Se arquivo não existe, imprimir que não existe
            if not exists:
                messages.append(f"{self.dir} - the output file doesn't exists")

            # Se arquivo existe
            else:

                # Se convergiu
                if converged:
                    
                    # Adicionar às mensagens que convergiu
                    messages.append(f"{self.dir} - {'Converged':<14} | Last step: {steps_written[-1]}")
                
                # Se foi EXPLICITAMENTE encontrado que não convergiu
                elif not_converged:

                    # Reportat a não convergência
                    messages.append(f"{self.dir} - {'NOT converged':<14} | Last step: {steps_written[-1]}")

                # Se não convergiu, mas não foi explicitamente encontrado que não convergiu
                else:

                    # Se tiver passos escritos
                    if len(steps_written) > 0:
                    
                        # Relatar que está indefinido
                        messages.append(f"{self.dir} - {'Undefined':<14} | Last step: {steps_written[-1]}")

                    # Se nenhum passo foi escrito no arquivo, reportar   
                    else:

                        # Relatar que está indefinido
                        messages.append(f"{self.dir} - {'Undefined':<14} | No steps written to file")

                # Se houve SCC não convergente
                if len(not_conv_SCC) != 0 and verbose:

                    # Calcular porcentagem não convergente, não considerando o 
                    # último passo, que ainda está sendo escrito
                    percent = round(100*len(not_conv_SCC_steps) / (len(steps_written) - 1), 1)

                    # Caso a leitura ocorra no pequeno intervalo que onde o passo acabou se ser escrito
                    # E o pŕoximo não se iniciou, gerando mais que 100%, jogar de volta para 100
                    percent = percent if (percent <= 100) else 100

                    # Nova linha
                    messages.append("")

                    # Adicionar à mensagem os não convergentes
                    messages.append(f"\tNot converged SCC's: {not_conv_SCC_steps} => {percent}% of total")

                    # Nova linha
                    messages.append("")

                    # Se o último passo da otimização teve SCC não convergente
                    # Verificar se o SCC divergiu no último passo
                    if not_conv_SCC_steps[-1] == steps[-1]:
                        
                        # Adicionar warning às mensagens
                        messages.append("WARNING: SCC was not converged in the final step.")

            # Imprimir todas as mensagens
            for message in messages:

                print(message)

        ######################### VERIFICANDO SE ARQUIVO EXISTE ###############

        # Se o arquivo não existir
        if not os.path.exists(output_file):

            # Se tiver em modo de retorno de passos escritos, retornar None
            if return_written:
                return None

            # Se tiver em modo de retorno da convergência, retornar que não 
            # convergiu
            if return_conv:
                return False
            
            # Imprimir status
            print_status(exists=False)

        #######################################################################

        # Do contrário...
        else:

            # Ler arquivo output
            with open(output_file) as file:
                data = file.read()

            ################### LISTANDO PASSOS DO SCC ########################

            # Listando posição de todos os passos => [12, 157, 368, ...]
            steps_indexes = find_all(data, "Geometry step:")
            # Criando lista dos números dos passos => [0, 1, 2, 3, 4, ...]
            steps = [ i for i in range(len(steps_indexes))]

            # Se estiver no modo de retorno dos passos escritos no arquivo
            if return_written:
                # Retornar passos
                return steps

            # Listando posição de SCC's não convergentes
            not_conv_SCC_indexes = find_all(data, "SCC is NOT converged")

            # Encontrando passos onde o SCC não convergiu
            not_conv_SCC_steps = indexes_to_steps(steps_indexes, 
                                                  not_conv_SCC_indexes)

            # Se for modo de retorno, retornar os SCC não conv
            if return_SCC:

                return not_conv_SCC_steps

            ############## VERIFICANDO ESTADO DA OTIMIZAÇÃO ###################

            # Se achar no arquivo, que a otmização convergiu...
            if data.find("Geometry converged") != -1:
                
                # Se estiver nesse modo de retorno de convergência, retornar 
                # que convergiu
                if return_conv:
                    return True

                # Imprimir status
                print_status(exists=True, 
                             converged=True,
                             not_conv_SCC=not_conv_SCC_indexes, 
                             steps_written=steps)

            # Se achar no arquivo que a otimização não convergiu, reportar
            elif (data.find("Geometry NOT converged") != -1 
                  or data.find("Geometry did NOT converge") != -1):
                
                # Se estiver nesse modo, retornar que não convergiu
                if return_conv:
                    return False
                
                # Imprimir status
                print_status(exists=True, 
                             converged=False,
                             not_converged=True,
                             not_conv_SCC=not_conv_SCC_indexes, 
                             steps_written=steps)
                
            # Se não foi encontrado que convergiu, mas também não foi 
            # encontrado que não convergiu
            else:

                # Se estiver nesse modo, retornar que não convergiu
                if return_conv:
                    return False
                
                # Imprimir status
                print_status(exists=True, 
                            converged=False, 
                            not_conv_SCC=not_conv_SCC_indexes,
                            steps_written=steps)
                
    ############# Dopagem

    def dope_out(self, dop_elem: str, id: int, kwargs: dict=None) -> structure:
        
        """
        Gera uma versão (dopada) de uma estrutura, onde o átomo 
        especificado é substituído por um de outro elemento.

        Parameters
        ----------

        dop_elem : str
            Elemento do átomo que substituirá o antigo.

        id : int
            ID do átomo que será substituído.

        **kwargs
            Passado para o construtor de structure.

        Returns
        -------

        structure
            Estrutura dopada, onde o átomo especificado foi substituído 
            por um de outro elemento. Obs: As cargas originais são 
            substituídas por None.

        Raises
        ------

        TypeError
            dop_elem must be string
        
        TypeError
            id must be int
        """

        #######################################################################

        # Verificando validade
        if not isinstance(dop_elem, str):
            raise TypeError("dop_elem must be string")
        
        if not isinstance(id, int):
            raise TypeError("id must be int")

        #######################################################################
        
        from copy import deepcopy
        
        # Criando cópias dos valores invés de copiar o ponteiro
        atoms = list(deepcopy(self.atoms))

        # Apagando cargas
        for atom in atoms:
            atom.charge = None

        # Modificando o símbolo de elemento químico no átomo especificado
        atoms[id-1].elem = dop_elem

        # Criando outro objeto structure, agora dopado
        doped_out = structure(atoms=atoms, **kwargs)

        # Retornando
        return doped_out

    ############# 

    def formation_energy(self) -> tuple[float, float]:

        """
        Retorna a energia de formação da estrutura, além de sua média 
        por átomo.

        Returns
        -------
        
        tuple[float, float]

            Os seguintes valores:

            Energia de formação da estrutura (eV) no último passo de 
            otimização, energia média de formação por átomo.
        """

        # Lista que receberá os dados
        data = []

        # Para cada átomo, anexar sua energia à lista "data"
        for atom in self.atoms:
            data.append(atoms_data["energ_atom"][self.param][atom.elem])

        # Retornando como array ao final da iteração
        atoms_energies_sum = np.array(data).sum()

        ############# Valência e energia dos átomos

        energy = round(self.total_energy - atoms_energies_sum, 4)

        # A energia total, menos a soma da energia de cada átomo
        return energy, round(energy / self.size, 4)

    def homo_lumo(self) -> tuple[float, float, float]: 
        """
        Retorna as energias de HOMO, LUMO, e o módulo da diferença 
        HOMO-LUMO de uma estrutura a partir de seu arquivo band.out.

        Parameters
        ----------
        None

        Returns
        -------

        tuple[float, float, float]

            Os seguintes valores (todos arredondados em 3 casas 
            decimais):
                
            homo, lumo, homo-lumo
        """

        # Retornar valores
        return self.homo, self.lumo, round(self.lumo - self.homo, 3)

    ##### Análise das ligações

    def shortest_distances(self, elem1: str, elem2: str, 
                           n: int) -> list[float]:

        """
        Retorna as "n" menores distâncias entre átomos de dados 
        elementos, na estrutura.

        Parameters
        ----------

        elem1 : str
            Símbolo do elementos do átomo 1.

        elem2 : str
            Símbolo do elementos do átomo 2.

        n : int
            Número de comprimentos retornados.

        Returns
        -------

        list[float]
            Lista com as "n" menores comprimentos entre átomos da estrutura.
        """
        
        # Lista para receber comprimentos de ligação
        data = []

        # Salvando quantidade de átomos na estrutura
        tam = self.atoms.shape[0]

        # Para o índice de cada átomo na estrutura, procuar o átomo buscado
        for i in range(tam):

            # Se dado átomo é o buscado
            if self.atoms[i].elem == elem1:

                # Para o índice de cada outro átomo na estrutura
                for j in range(tam):

                    # Se dado átomo é o buscado
                    if self.atoms[j].elem == elem2:

                        # Calcular comprimento da ligação
                        bond_length = self.atoms[i].dist_to(self.atoms[j])

                        # Listar valor do comprimento 
                        data.append(bond_length)

        # Retornar primeiros "n" resultados da lista ordenada
        return sorted(data)[0:n]

    ############ Escrita/cópia de arquivos

    def frame(self, frame_out: str|Path) -> None:
        """
        Gera e escreve o arquivo xyz do último frame da estrutura.
        Cria diretórios necessários caso não existam.

        Parameters
        ----------
        frame_out: str or Path
            Endereço de saída do frame.

        Raises
        ------

        TypeError
            "frame_out must be str or Path."

        ValueError
            This structure has no atoms.
        """

        # Se frame_out não for string, nem Path, apontar erro
        if not isinstance(frame_out, str|Path):
            raise TypeError("frame_out must be str or Path.")

        # Lista para receber dados
        data = []

        # Adicionando cabeçalho - Quantidade de átomos da estrutura
        data.append(str(self.size))

        # Imprimir nome do diretório na segunda linha do cabeçalho
        if self.dir is not None:
            data.append(self.dir.name)
        else:
            data.append("frame")

        if len(self.atoms) > 0:
            # Para cada átomo, adicionar uma linha com: elemento, X, Y, Z 
            if self.atoms[0].charge is not None:           
                for atom in self.atoms:
                    data.append(" ".join( [ atom.elem ] + [ str(round(number, 8)) for number in atom.coord] + [str(atom.charge)]))
            else:
                for atom in self.atoms:
                    data.append(" ".join( [ atom.elem ] + [ str(round(number, 8)) for number in atom.coord]))
        else:
            raise ValueError("This structure has no atoms.")

        # Convertendo tudo em texto (string formatada)
        data = "\n".join(data)

        # Convertendo endereço em caminho Path
        frame_out = Path(frame_out).expanduser()
        # Criando diretórios necessários para a escrita
        Path.mkdir(frame_out.parent, exist_ok=True, parents=True)

        # Escrevendo em arquivo no endereço passado
        with open(frame_out, mode = 'w') as file:
            file.write(data)

        # Imprimir endereço do arquivo escrito
        print(frame_out)

    def output(self, output_dest: str|Path) -> None:

        """
        Copia o output de otimização da estrutura para o endereço 
        passado.

        Parameters
        ----------

        output_dest : str | Path
            Destino do arquivo output.
        """

        from shutil import copyfile

        # Convertendo endereço em objeto Path
        output_origin = self.dir / main_config["output_name"]
        output_dest = Path(output_dest).expanduser()

        if not output_origin.is_file():
            raise AssertionError("output file not found.")

        # Criando diretórios necessário
        Path.mkdir(output_dest.parent, parents=True, exist_ok=True)

        # Coopiando arquivo
        copyfile(output_origin, output_dest)

        

        # Imprimir endereço do arquivo escrito
        print(output_dest)

    def copy(self) -> structure:
        
        """
        Retorna uma estrutura que contém os mesmo átomos que a original.
        """

        return structure(atoms=self.atoms)

####################################################################################
