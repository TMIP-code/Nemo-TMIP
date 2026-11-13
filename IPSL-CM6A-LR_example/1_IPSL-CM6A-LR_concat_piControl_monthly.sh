#!/bin/bash -l  
# Use the current working directory
#SBATCH -D ./
# Use the current environment for this job.
##SBATCH --export=ALL
# Define job name
#SBATCH -J cdo_timeavg
# Define a standard output file. When the job is running, %u will be replaced by user name,
# %N will be replaced by the name of the node that runs the batch script, and %j will be replaced by job id number.
#SBATCH -o ./SLURM_OUT/cdo_timeavg.%u.%N.%j.out
# Define a standard error file
#SBATCH -e ./SLURM_OUT/cdo_timeaavg.%u.%N.%j.err
# Request the partition
#SBATCH -p nodes 
# Request the number of nodes
#SBATCH -N 1 
# Request the total number of cores
#SBATCH -n 1
# This asks for 0 days, 1 hour, 0 minutes and 0 seconds of time.
#SBATCH -t 00:59:59
# Specify memory
#SBATCH --mem=30000M


module load cdo/2.5.4-openmpi5.0.8-gcc14.2.0

#Use cdo to generate monthly climatologies from outputs downloaded from esgf
#here, I'm taking the climatology over 2050-2099 in the piControl experiment for r1i2p1f1

cdo ymonmean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/mlotst_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_gn*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/mlotst.nc"
cdo ymonmean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/umo_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_gn*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/umo.nc"
cdo ymonmean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/vmo_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_gn*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/vmo.nc"
cdo ymonmean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/agessc_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/agessc.nc"
cdo ymonmean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/thetao_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_gn*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/thetao.nc"
cdo ymonmean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/so_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_gn*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/so.nc"
cdo ymonmean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/uo_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_gn*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/uo.nc"
cdo ymonmean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/vo_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_gn*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/vo.nc"

#check metric files to see if they are static (Ofx) or time evolving (Omon)
cp ../data/raw/areacello_Ofx_IPSL-CM6A-LR_piControl_r1i1p1f1_gn.nc ../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/areacello.nc

#cp ../data/raw/volcello_Ofx_IPSL-CM6A-LR_piControl_r1i2p1f1_gn.nc ../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/volcello.nc
cdo timemean -seldate,2050-01-01,2099-12-31 -mergetime "../data/raw/thkcello_Omon_IPSL-CM6A-LR_piControl_r1i2p1f1_gn*.nc" "../data/processed/cyclo_stationary/IPSL-CM6A-LR_r1i2p1f1_monthly/thkcello.nc"
