from pathlib import Path

from dopings.structure import structure
from dopings.config import atoms_data, dops_data, dirs_data, main_config

####################################################################################
        
class dops_set():
    """
    Classe de objetos que representam um conjunto completo de dopagens,
    com atributos e métodos para gerar seus arquivos, organizar 
    diretórios, otimizar, e processar seus resultados em massa. 

    Parameters
    ----------
    
    structs : list[structure]
        Lista de estruturas do conjunto.

    param : str
        Conjunto de dopagens usado nas otimizações.

    main_dir : Path
        Diretório raíz de resultados de processamento dos dados.
        
    bases : { str : structure }
        Dicionário que relaciona o nome da base com a structure da 
        mesma.

    dops_info : dict
        Dicionário com dados gerais do conjunto de dopagens.

        self.dops_info = {

            "material"  : list[str],
            "bases"     : list[str]
            "sites"     : { str : list[str] },
            "sites_id"  : { str : { str : int }},
            "n_atoms"   : int,
            "n_ligs"    : int
        }

    mode : { "look", "read", "write" }
        In whick mode dops_set object was initialized.

    Methods
    -------

    
    __init__(self, dop_elem: str, dops_info: dict=None, mode: str="read", param=None, force_write: bool=False) -> None

        Contrutor de objetos dops_set.
    
    map(self, method: "function", suffix: str, **kwargs) -> None

        Executa uma função para cada estrutura do conjunto de 
        estruturas.
    
    map_to_csv(self, method: "function", **kwargs) -> None

        Executa uma função para cada estrutura do conjunto, recolhendo 
        os retornos para escrever em um arquivo csv final.
    
    map_opt(self, only_report: bool=True, overwrite: bool=False, verbose: bool=True, skip_hard_to_conv_SCC: bool=True, resume_unfineshed: bool=False, inverse_order: bool=False) -> None

        Otimiza um conjunto inteiro de dopagens em fila.
    
    atom_id(self, base, site: dict) -> int

        Retorna o ID do átomo correspondente ao sítio de dopagem da 
        estrutura.
    
    infer_base_dir(self, base) -> Path

        Infere o diretório da estrutura base que origina a estrutura 
        dopada em questão.
    
    infer_dir_out(self, struct, method, suffix=None, as_csv=False) -> Path

        Faz a inferência do diretório de escrita dos arquivos, baseado
        em outros dados do conjunto. 
    
    valid_dops_info(dops_info) -> bool

        Valida um dicionário dops_info, que descreve a estrutura do
        conjunto de dopagens.
    
    infer_dops_dir(material, dop_elem, site, base) -> Path

        Faz inferência do diretório de dopagem da estrutura.
    
    gen_files(struct, read_initial_charges=False, redo=False) -> None

        Gera os arquivos de otimização de uma estrutura.
    """

    def __init__(self, dop_elem: str, dops_info: dict=None, mode: str="read",
                 param=None, force_write: bool=False) -> None:

        """
        Contrutor de objetos dops_set.

        Parameters
        ----------

        dop_elem : str
            Símbolo do elemento substituído.

        dops_info : dict
            DIct com os dados do conjunto de dopagens, com a 
            estrutura:

            self.dops_info = {

                "material"  : list[str],
                "bases"     : list[str]
                "sites"     : { str : list[str] },
                "sites_id"  : { str : { str : int }},
                "n_atoms"   : int,
                "n_ligs"    : int
            }

        mode : {"look", "read", "write"}
            Modo de acesso às estruturas. 
            
            No modo look as estruturas são apenas listadas, para que se
            possa emitir report independente do status da otimização.

            No modo read, os dados das estruturas são lidos dos arquivos
            e carregados nas propriedades estrutura. Obs: Estruturas não
            convergidas não são adicionadas ao dops_set, porém um 
            warning é impresso.

            No modo write, os arquivos de otimização do conjunto de 
            dopagens é criado, como preparação para otimização do 
            conjunto. Porém se no local existir um diretório de reinício
            de otimização, a estrutura é pulada. Esse comportamento pode
            ser sobrescrito com o parâmetro "force_write".

        param : str
            Conjunto de parâmetro usados na otimização.

        force_write : bool, default = False
            Se True, novos arquivos de otimização são gerados para todas
            as estruturas, mesmo que para elas já exista um diretório 
            de reinício de otimização. 
            
        Raises
        ------

        TypeError
            param must be str.

        ValueError
            mode must be 'read', look or 'write'.

        ValueError
            Not valid dops_info.

        ValueError
            dops_info cannot be None.

        ValueError
            Invalid dop_elem.
        
        ValueError
            dop_elem can't be None.
        """

        import os

        # Propriedades
        self.main_dir = None
        self.param = None
        self.bases = {}
        self.structs = []
        self.dops_info = None
        self.mode = None

        # Dop_elem
        if dop_elem is not None:
            if isinstance(dop_elem, str) and len(dop_elem) in [1, 2]:
                self.dop_elem = dop_elem
            else:
                raise ValueError("Invalid dop_elem.")
        else:
            raise ValueError("dop_elem can't be None.")

        # Definindo parâmetro
        if param is not None:
            if isinstance(param, str):
                self.param=param
            else:
                raise TypeError("param must be str.")
        else:
            # Se não for passado, buscar no dicionário de parâmetros
            if self.dop_elem is not None and self.dop_elem in dops_data["param_dict"].keys():
               self.param = dops_data["param_dict"][self.dop_elem]

        # Guardar endereço principal
        self.main_dir = dirs_data["processing_output"]

        # Mode
        if mode is not None:
            if mode in ["read", "write", "look"]:
                self.mode = mode
            else:
                raise ValueError("mode must be 'read', 'look' or 'write'")
        
        # Dops_info
        if dops_info is not None:
            if dops_set.valid_dops_info(dops_info):
                self.dops_info = dops_info
            else:
                raise ValueError("Not valid dops_info")
            
        else:
            raise ValueError("dops_info cannot be None")

        ## Base struts
        for base in self.dops_info["bases"]:

            # Gerando bases
            base_struct = structure(dir = self.infer_base_dir(base),
                                    read_from_dir=True)
            
            self.bases[base] = base_struct

        # Para cada estrutura base
        for base in self.dops_info["bases"]:

            # Para cada site
            for site in self.dops_info["sites"][base]:

                kwargs = {

                    "material"       : dops_info["material"],
                    "dop_elem"       : self.dop_elem, 
                    "base"           : base,
                    "site"           : site,
                }
                
                opt_dir = self.infer_dops_dir(**kwargs)

                kwargs["dir"] = opt_dir

                replace_atom_id = self.atom_id(base, site) 

                # Gerar estrutura dopada
                doped_struct = self.bases[base].dope_out(
                                                dop_elem=self.dop_elem, 
                                                id=replace_atom_id,
                                                kwargs=kwargs)
                
                # Definindo parâmetro
                doped_struct.param = self.param

                if mode == "look":

                    # Anexar à lista
                    self.structs.append(doped_struct)

                # Se modo de escrita
                elif mode == "write":

                    # Anexar à lista
                    self.structs.append(doped_struct)

                    # Se o arquivo de backup existe
                    backups_exists = (doped_struct.dir / main_config["resume_dir_name"]).is_dir()

                    # Se o backup não existe, ou foi pedido para forçar a escrita
                    if (not backups_exists) or force_write:
                        self.gen_files(doped_struct)

                    else:
                        print(f"WARNING: not writing over {doped_struct.dir}, cause this optimization was restarted")

                # Se modo leitura
                elif mode == "read":

                    # Se não tiver convergido, pular iteração
                    if not doped_struct.report(return_conv=True):

                        print(doped_struct.dir / main_config["output_name"])
                        if not os.path.exists(doped_struct.dir / main_config["output_name"]):
                            print(f"Skipping structure with not found output: {doped_struct.dir}")
                        
                        else:
                            print(f"Not reading not converged structure: {doped_struct.dir}")
                        
                        continue
                    
                    # Sobrescrever com estrutura dopada
                    doped_struct = structure(**kwargs, read_from_dir=True)

                    # Definindo parâmetro
                    doped_struct.param = self.param

                    # Anexar à lista
                    self.structs.append(doped_struct)

    # Mapeia funções
    def map(self, method: "function", suffix: str, **kwargs) -> None:
        """
        Executa uma função para cada estrutura do conjunto de 
        estruturas.
        
        Parameters
        ---------

        method : function
            Função a ser executada em todas as estruturas.
        
        suffix : str or Path
            Extensão do arquivo de saída.


        Raises
        ------
        ValueError
            read mode is needed for proper processing results.
        """

        if self.mode != "read":
            raise ValueError("read mode is needed for proper processing"
                             + "results.")

        # Para cada structure
        for struct in self.structs:

            # Checagem de convergência com aviso no nome do arquivo?
            # Diretório de testes
            dir_out=self.infer_dir_out(struct=struct, method=method, 
                                       suffix=suffix)

            # Executando função
            method(struct, dir_out, **kwargs)

    # Mapeia funções a apresenta retornos em csv
    def map_to_csv(self, method: "function", **kwargs) -> None:
        """
        Executa uma função para cada estrutura do conjunto, recolhendo 
        os retornos para escrever em um arquivo csv final.

        Parameters
        ----------

        method : function
            Função a ser executada em todas as estruturas.

        Raises
        ------

        AssertionError
            There is no structs on this struct_set.
        """

        # Lista para receber valores
        values = []

        # Se não tiver estruturas, reportar
        if len(self.structs) == 0:
            raise AssertionError("There is no structs on this set_struct.")

        # Para cada structure
        for struct in self.structs:

            # Executar função e guardar valor
            value = method(struct, **kwargs)

            if hasattr(self, "name"):
                
                # Dict de identificação
                id = { "Name" : self.name, "Base" : struct.dir.parent.name, "Site" : struct.dir.name }

            # Se não tiver
            else:
                
                material = struct.material
                dop_elem = struct.dop_elem
                site = struct.site

                # Dict de identificação
                id = { "material" : material, "dop_elem" : dop_elem, "Base" : struct.dir.parent.name, "Site" : site }

            # Formatação para dados em csv
            if method.__name__ == "homo_lumo":

                value = {

                    "homo" : value[0],
                    "lumo" : value[1],
                    "homo_lumo" : value[2]
                }
            
            # Formatação para dados em csv
            elif method.__name__ == "formation_energy":

                value = {

                    "energy" : value[0],
                    "per_atom" : value[1],
                }
            
            # Juntando identificação e valores 
            values.append(dict(id, **value))

        # Diretório final
        dir_out=self.infer_dir_out(struct=struct, method=method, as_csv=True)
                                   
        # Criando diretório
        Path.mkdir(dir_out.parent, exist_ok=True, parents=True)

        # Escrevendo arquivo CSV
        import csv
        with open(dir_out, 'w') as file:
            # Escreve dicionário como csv, usa como nomes de coluna as chaves de um item
            writer = csv.DictWriter(file, fieldnames=values[0].keys())
            # Escreve cabeçalho, e escreve valores
            writer.writeheader()
            writer.writerows(values)

        # Relatar
        print(dir_out)

    def map_opt(self, only_report: bool=True, overwrite: bool=False, 
                verbose: bool=True, skip_hard_to_conv_SCC: bool=True, 
                resume_unfineshed: bool=False, inverse_order: bool=False
                ) -> None:

        """
        Otimiza um conjunto inteiro de dopagens em fila.

        Parameters
        ----------

        only_report : bool, default = True
            Não realiza otimização, apenas faz o report do estado dos
            dos arquivos.

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

        inverse_order : bool, default=False
            Se a ordem das otimizações devem seguir a ordem inversa da 
            lista de estruturas.
        """


        # Salvando lista de estruturas
        structs = self.structs

        # Se for pra rodar na ordem inversa
        if inverse_order:
            # Inverter lista
            structs = structs[::-1]

        # Pra cada estrutura...
        for struct in structs:

            if only_report:
                struct.report(verbose=verbose)

            else:
                struct.opt(overwrite=overwrite, verbose=verbose, 
                           skip_hard_to_conv_SCC=skip_hard_to_conv_SCC, 
                           resume_unfineshed=resume_unfineshed)
       
    ###########################################################################

    def atom_id(self, base, site: dict) -> int:
        """
        Retorna o ID do átomo correspondente ao sítio de dopagem da 
        estrutura.

        Parameters
        ----------

        base : str
            Base da estrutura.

        site : str
            Sítio de dopagem da estrutura.
        
        Returns
        -------

        int
            ID do átomo buscado.
        """

        return int(self.dops_info["sites_id"][base][site])

    ###########################################################################

    def infer_base_dir(self, base) -> Path:
        """
        Infere o diretório da estrutura base que origina a estrutura 
        dopada em questão.

        => Formato

        [bases root dir] / [material] / [base]

        Parameters
        ----------

        base : srt
            Base da qual se busca o diretório.

        Returns
        -------

        Path
            Diretório da base correspondente.
        """

        return Path(dirs_data["bases"]).expanduser() / self.dops_info["material"] / base

    def infer_dir_out(self, struct, method, suffix=None, as_csv=False) -> Path:

        """
        Faz a inferência do diretório de escrita dos arquivos, baseado
        em outros dados do conjunto. 
        
        Os diretórios seguem a seguinte estrutura:

        => Diretório CSV

        [raíz] / [method name] / [method name]-[nome do conjunto].csv  
        
        => Diretório dos arquivos individuais

        [raíz] / [method name] / [nome do conjunto] / [nome do conjunto]
        -[nome do diretório mãe]-[nome do diretório].[sufixo]
        
        Parameters
        ----------

        struct : structure
            Estrutura que terá diretório inferido.

        method : function
            Função a ser mapeada nas estruturas.

        suffix : str
            Extensão do arquivo de saída 

        as_csv : bool
            Se for True, infere diretório de um arquivo csv, que 
            representa o conjunto como um todo.
        
        Returns
        -------

        Path
            Caminho inferido para o arquivo de saída. 
        """

        material = struct.material
        base = struct.base
        dop_elem = struct.dop_elem
        site = struct.site

        if as_csv:            

            # ROOT / METHOD / FILE(method+name)
            return (self.main_dir.expanduser() / material / method.__name__  / f"{method.__name__}-{dop_elem}").with_suffix(".csv")

        else:
            
            # Se não tiver sufixo, reportar
            if suffix is None:
                raise ValueError("suffix is needed")
            
            # ROOT / METHOD / NAME / FILE(dir_name)
            return (self.main_dir.expanduser() / material / method.__name__  / dop_elem / f"{dop_elem}-{base}-{site}").with_suffix(suffix)

    @staticmethod
    def valid_dops_info(dops_info) -> bool:

        """
        Valida um dicionário dops_info, que descreve a estrutura do
        conjunto de dopagens.

        Tal dicionário deve seguir o seguinte formato:

        self.dops_info = {

            "graphs_ylims"  : { str : { str : list[float] }}
            "material"      : list[str],
            "bases"         : list[str],
            "sites"         : { str : list[str] },
            "sites_id"      : { str : { str : int }},
            "n_atoms"       : int,
            "n_ligs"        : int
        }

        Parameters
        ----------

        dops_info : dict
            Dicionário dops_info a ser validado.

        Returns
        -------

        bool
            Se é o dicionário é valido ou não.
        """

        dops_info_model_keys = [

            "graphs_config",
            "material",
            "bases",
            "sites",
            "sites_id",
            "n_atoms",
            "n_ligs"
        ]

        # Se as chaves principais forem diferentes, não é válido
        if not sorted(dops_info.keys()) == sorted(dops_info_model_keys):
            
            return False
        
        # Se as chaves de "sites" forem diferentes das "bases"
        if not sorted(dops_info["sites"].keys()) == sorted(dops_info["bases"]):
            
            return False
        
        # Se o mesmo não ocorrer com sites_id
        if not (sorted(dops_info["sites_id"].keys()) == 
                sorted(dops_info["bases"])):

            return False
        
        # Se para cada base do "sites_id" os keys de sites não correspondetem 
        # aos sites correspondente
        if not all([(sorted(dops_info["sites_id"][base].keys()) == 
              sorted(dops_info["sites"][base])) for 
              base in dops_info["bases"]]):

            return False

        # Se resistir a todos os testes, é válido
        return True

    @staticmethod
    def infer_dops_dir(material, dop_elem, site, base) -> Path:
        """
        Faz inferência do diretório de dopagem da estrutura.

        Parameters
        ----------

        struct : structure
            Estrutura sobre a qual se quer fazer a inferência.
        """

        # Retornar caminho da dopagem
        return Path(dirs_data["dopings_opt"]).expanduser() / material / dop_elem / base / site

    @staticmethod
    def gen_files(struct, read_initial_charges=False, redo=False) -> None:

        """
        Gera os arquivos de otimização de uma estrutura.

        Parameters
        ----------

        struct : structure
            Estrutura para a qual deseja-se gerar arquivos de 
            otimização.

        read_initial_charges : bool
            Se deve ser especificado no arquivo hsd que as cargas 
            iniciais devem ser lidas a partir de arquivo antes de 
            iniciar a otimização.

        redo : bool
            Se os arquivos estão sendo refeitos, e portanto o limite de
            passos na otimização deve ser maior que o normal.
        """

        ################### Funções para lidar com arquivo hsd ################

        def hsd_to_dict(dir_in):

            import hsd

            # Lendo arquivo genérico
            with open(dir_in) as file:

                # Lendo para um dicionário
                hsd_dict = hsd.load(file)

            # Retornando dicionário
            return hsd_dict

        def write_hsd_dict(hsd_dict, dir_out):

            import hsd

            # Carregando informação como string formatada
            txt = hsd.dump_string(hsd_dict)
            # Substituindo informação desejada
            txt = txt.replace('xyzFormat {}', 'xyzFormat {\n    <<< "inp.xyz"\n  }')

            # Escrevendo em arquivo
            with open(dir_out, "w") as file:
                file.write(txt)

        #################### Gerar frame ######################################

        struct.frame(struct.dir / "inp.xyz")
        
        #################### Gerar hsd ########################################

        # Target
        target = struct.dir / "dftb_in.hsd"

        # Lendo para um dicionário o modelo principal
        if struct.dop_elem != "Li":
            hsd_dict = hsd_to_dict(dirs_data["generic_hsd"][struct.param])

        else:
            hsd_dict = hsd_to_dict(dirs_data["generic_hsd"]["3ob+Li"])


        # Obtendo valores únicos de elementos na estrutura
        elements = list(set([ atom.elem for atom in struct.atoms ]))

        # Importando dados
        max_ang_momentum = atoms_data["max_ang_momentum"]

        # Definindo momentos angulares máximos para cada elemento
        for elem in elements:

            # Criando entrada, com momento angular entre aspas
            hsd_dict["Hamiltonian"]["DFTB"]["MaxAngularMomentum"][elem] = f'"{max_ang_momentum[elem]}"'
            
        # Se for 3ob, adicionar as derivadas de hubbard para cada elemento
        if struct.param == "3ob" and struct.dop_elem != "Li":
            for elem in elements:
                hsd_dict["Hamiltonian"]["DFTB"]["HubbardDerivs"][elem] = atoms_data["hubbard_derivs"][elem]
                
        # SKF DIR
        skf_dir = dirs_data["skf_dirs"][struct.param]
        hsd_dict["Hamiltonian"]["DFTB"]["SlaterKosterFiles"]["Type2FileNames"]["Prefix"] = '"' + str(skf_dir.resolve()) + '/"'

        # Temperatura como forma de melhorar a convergência
        hsd_dict["Hamiltonian"]["DFTB"]["Filling"] = {"Fermi" : {"Temperature[K]" : 2000}}

        ############################################################################

        # Ler cargas iniciais, se for passado como parâmetro
        if read_initial_charges:

            # Ler cargas iniciais
            hsd_dict["Hamiltonian"]["DFTB"]["ReadInitialCharges"] = "Yes"

            # Ler cargas como texto
            hsd_dict["Options"]["ReadChargesAsText"] = "Yes"

        ############################################################################

        # Se estiver refazendo os arquivos pra resumir a otimização
        if redo:

            # Usar como limite de passos o valor designado nas variáveis de configuração
            hsd_dict["Driver"]["GeometryOptimization"]["MaxSteps"] = main_config["redo_steps_limits"]

        ############################################################################

        # Escrevendo arquivo de propriedades
        write_hsd_dict(hsd_dict, target)

        # Relatar
        print(target)
