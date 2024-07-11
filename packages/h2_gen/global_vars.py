import json
from pathlib import Path

# Diretório do arquivo de dados dos átomos
config_dirs = {

    "atoms_data_dir"        : "../INPUT/config_files/atoms_data.json",
    "config_file_dir"       : "../INPUT/config_files/viz_config.json",
    "dops_set_data_dir"     : "../INPUT/config_files/dops_data.json",
    "graphine_data"         : "../INPUT/config_files/graphine_data.json",
    "dirs_data_dir"         : "../INPUT/config_files/dirs_data.json",
    "main_config"           : "../INPUT/config_files/main_config.json",
    "h2_gen_data"           : "../INPUT/config_files/h2_gen_data.json"
}

# Lê diretórios e converte em Path
with open(config_dirs["dirs_data_dir"]) as file:
    dirs_data = json.load(file)

    for key in dirs_data.keys():

        # Transformando primeiro nível em Path
        if isinstance(dirs_data[key], str):
            dirs_data[key] = Path(dirs_data[key]).expanduser()

        # Tranformando segundo nível em Path
        if isinstance(dirs_data[key], dict):
            for key2 in dirs_data[key].keys():
                if isinstance(dirs_data[key][key2], str):
                    dirs_data[key][key2] = Path(dirs_data[key][key2]).expanduser()

# Lê os dados de configurações gerais
with open(config_dirs["main_config"]) as file:
    main_config = json.load(file)

# Lê os dados de configurações gerais
with open(config_dirs["graphine_data"]) as file:
    graphine_data = json.load(file)

# Lê os dados do conjunto de dopagens usadas
with open(config_dirs["dops_set_data_dir"]) as file:
    dops_data = json.load(file)

# Lê dados do arquivo
with open(config_dirs["atoms_data_dir"]) as file:
    atoms_data = json.load(file)

# Lê dados do arquivo
with open(config_dirs["config_file_dir"]) as file:
    config = json.load(file)

# Lê dados do arquivo
with open(config_dirs["h2_gen_data"]) as file:
    h2_gen_data = json.load(file)