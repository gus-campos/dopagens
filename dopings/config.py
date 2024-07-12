import json
from pathlib import Path

# Para build local dos docs
#config_files_dir = Path("input/config_files")
# Para build no Read The Docs
config_files_dir = Path("../../input/config_files")

# Para uso nos scripts
#config_files_dir = Path("../input/config_files")

###############################################################################

# Lê diretórios e converte em Path
with open(config_files_dir / "dirs_data.json") as file:
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

with open(config_files_dir / "main_config.json") as file:
    main_config = json.load(file)

with open(config_files_dir / "graphine_data.json") as file:
    graphine_data = json.load(file)

with open(config_files_dir / "dops_data.json") as file:
    dops_data = json.load(file)

with open(config_files_dir / "atoms_data.json") as file:
    atoms_data = json.load(file)

with open(config_files_dir / "viz_config.json") as file:
    viz_config = json.load(file)

with open(config_files_dir / "h2_gen_data.json") as file:
    h2_gen_data = json.load(file)