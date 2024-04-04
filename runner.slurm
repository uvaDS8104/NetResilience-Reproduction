#!/bin/bash

#SBATCH -A shakeri-lab
#SBATCH -p gpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=32000
#SBATCH -c 1
#SBATCH -t 72:00:00
#SBATCH -o ./logs/%A.out
#SBATCH -e ./logs/%A.err

# set up environment
CONDA_ENV_PATH="/scratch/jgh2xh/envs/netres12"
PIP_EX="$CONDA_ENV_PATH/bin/pip"
PYTHON_EX="$CONDA_ENV_PATH/bin/python"

# confirm python executable is from environment
echo $PYTHON_EX $@
$PYTHON_EX $@

echo
echo "All done"