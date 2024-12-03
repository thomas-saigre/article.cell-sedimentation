#!/bin/bash
#SBATCH -J rotation        # name of the job
#SBATCH -N 1                    # number of nodes
#SBATCH --ntasks-per-node=1   # number of MPI tasks per node
#SBATCH --threads-per-core=1    # no hyperthreading
#SBATCH -t 48:00:00             # walltime
#SBATCH -p public               # partition
#SBATCH --array=0-0
##SBATCH --exclusive             # exclusive mode

#SBATCH -o log/%j-rotation-o.log # standard output
#SBATCH -e log/%j-rotation-e.log # standard error


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
mpiexec -bind-to core /data/home/saigre/develop/build/default/research/article.cell-sedimentation.arvo2025/src/feelpp_eye_heatfluid_wss \
        --directory=${CURRENT_DIR}/feelppdb/rotation/${MESH_INDEX}-${DISC}/${TAMB} \
        --case.discretization=${DISC} \
        --config-files eye.cfg pc_${SOLVER_TYPE}.cfg ts.cfg \
        --heat-fluid.json.filename=\$cfgdir/eye.json \$cfgdir/mesh.json \$cfgdir/bc-velocity.json \$cfgdir/post_process.json \$cfgdir/initial_condition.json \
        --heat-fluid.json.patch='{ "op": "replace", "path": "/Meshes/heatfluid/Import/filename", "value": "$cfgdir/mesh/Mr/'${MESH_INDEX}'/Eye_Mesh3D_p$np.json" }' \
        --heat-fluid.json.patch='{ "op": "replace", "path": "/Parameters/T_amb", "value": "'${TAMB}'" }' \
        --heat-fluid.scalability-save=1 --heat-fluid.heat.scalability-save=1 --heat-fluid.fluid.scalability-save=1
