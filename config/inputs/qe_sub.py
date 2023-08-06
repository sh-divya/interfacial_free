slurm_str = """
date
source /data/apps/go.sh
module load intel/2020.2
module load intel-mpi/2020.2
module load quantum-espresso/6.6
mpirun -np $CORES$ pw.x -npool $NK$ < $INPUT$ > $OUTPUT$
date
"""
# SBATCH --job-name="$NAME$"
# SBATCH --output="$NAME$.o%j"
# SBATCH --nodes=$NODE$
# SBATCH --ntasks=$CORES$
# SBATCH --partition=defq
# SBATCH --time=$TIME$
# SBATCH --account=pclancy3
#   starting_magnetization(1) =   4.1666666667d-01

