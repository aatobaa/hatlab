#!/bin/bash
#SBATCH --job-name=example_sbatch_bp
#SBATCH --output=example_sbatch_bp.out
#SBATCH --error=example_sbatch_bp.err
#SBATCH --time=24:00:00
#SBATCH --partition=sandyb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load matlab
matlab -nodisplay -r 'analyze_BAO_LEARN M1TM_20111014; exit;'
