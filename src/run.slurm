#!/bin/bash
#SBATCH -J rot-init       # name of the job
#SBATCH -N 1                    # number of nodes
#SBATCH --ntasks-per-node=1   # number of MPI tasks per node
#SBATCH --threads-per-core=1    # no hyperthreading
#SBATCH -t 48:00:00             # walltime
#SBATCH -p public               # partition
#SBATCH --array=0-0
##SBATCH --exclusive             # exclusive mode

#SBATCH -o log/%j-init-o.log # standard output
#SBATCH -e log/%j-init-e.log # standard error


TAMBARRAY=(294 283 306)

CURRENT_DIR=$(pwd)


TAMB=${TAMBARRAY[$SLURM_ARRAY_TASK_ID]}

echo "TAMB=${TAMB}"

export OMP_NUM_THREADS=1
module load hpcx

MESH_INDEX=M1
SOLVER_TYPE=simple
DISC=P1-P2P1

# ts.cfg
mpiexec -bind-to core feelpp_toolbox_heatfluid \
        --directory=${CURRENT_DIR}/feelppdb/init/${MESH_INDEX}-${DISC}/${TAMB} \
        --case.discretization=${DISC} \
        --config-files eye.cfg pc_${SOLVER_TYPE}.cfg \
        --heat-fluid.json.patch='{ "op": "replace", "path": "/Meshes/heatfluid/Import/filename", "value": "$cfgdir/mesh/Mr/'${MESH_INDEX}'/Eye_Mesh3D_p$np.json" }' \
        --heat-fluid.json.patch='{ "op": "replace", "path": "/Parameters/T_amb", "value": "'${TAMB}'" }' \
        --heat-fluid.scalability-save=1 --heat-fluid.heat.scalability-save=1 --heat-fluid.fluid.scalability-save=1
