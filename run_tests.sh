#!/bin/bash
# Obtiene el directorio donde está este script (la carpeta clustRX)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Exporta PYTHONPATH para que los tests encuentren el paquete "clustrx"
export PYTHONPATH="$DIR"

# Ejecuta pytest usando el binario de tu conda env
/home/mario/miniconda3/envs/amp_miner_env/bin/python3 -m pytest "$DIR/tests/" -v
