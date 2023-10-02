#!/bin/bash

#SBATCH --mail-user=
#SBATCH --job-name=waterdimer_training
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=2000
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodelist=gpu10

hostname

module load gcc/gcc4.8.5-openmpi1.10-cuda9.2

#input file
input=run_waterdimer.inp

#restart folder (in case a run is restarted)
restart=

#datasets folder
datasets=datasets
#neural network code
neural_network=neural_network
#training handling code
training=training
#atomlabels file
atom_labels=atom_labels.tsv
#pythonscript for training
train=train.py
#environment
envi=env-gpu #GPU
#envi=env     #CPU
 
startfolder=`pwd`
scratch=/scratch/$USER.$SLURM_JOB_ID

#create scratch folder
if [ -d "$scratch" ]; then
   echo "scratch directory exists already"
else
   echo "creating scratch directory"
   mkdir $scratch
fi

#copy existing restart folder if present
if [ -d "$restart" ]; then
   cp -r $restart $scratch
fi

#link/copy data to scratch folder and go there
cp -r $train $scratch
cp -r $input $scratch
cp -r $atom_labels $scratch
cp -r $neural_network $scratch
cp -r $training $scratch
ln -s /home/kaeser/phd_projects/transfer_learning_comparison/b3lyp/nn_training1/$datasets $scratch/$datasets 
cd $scratch

#make necessary folders and load environment
source $HOME/$envi/bin/activate

#run actual jobs
#export CUDA_VISIBLE_DEVICES=0
./$train @$input 
cp -r $scratch $startfolder

#remove scratch folder
#cd $startfolder
#rm -r $scratch


