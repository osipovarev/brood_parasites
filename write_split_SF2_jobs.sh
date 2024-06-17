#!/usr/bin/env bash

#CHROMS=../../Genome/chrom.sizes
CHROMS=$1
RECMAP=$2


for C in $(awk '$2>1000000 {print $1}' $CHROMS );
do
 # end=$(tail -1 $C/header.$C.SF2.input | cut -f1);
 # last=$(tail -1 $C/norecmap.header.$C.SF2.out | cut -f1 | cut -d. -f1);

 # if [ $end != $last ];
 # then
    
    awk 'NR%250000==1 { file = FILENAME "_" sprintf("%04d", NR+249999) } { print > file }' $C/header.$C.SF2.input;
    for N in $(ls $C/header.$C.SF2.input_* | rev |  cut -d_ -f1 | rev); 
    do
      sed -i '1i position\tx\tn\tfolded' $C/header.$C.SF2.input_${N}
      
      cat template_sbatch.sh <(echo -e "fill_in_coord_gaps.py -f <(zgrep $C $RECMAP| cut -f2,4 | sort -n -k2) > $C/$C.recombination_map.txt")  > $C/sbatch.$C.recomb.sh; \
      cat template_sbatch.sh <(echo -e "SweepFinder2 -lg  1000 $C/header.$C.SF2.input_${N} $C/$C.SF2.spect $C/norecmap.header.$C.SF2.$N.out") > $C/sbatch.$C.$N.sweepfinder.norecmap.sh;
      cat template_sbatch.sh <(echo -e "SweepFinder2 -lrg 1000 $C/header.$C.SF2.input_${N} $C/$C.SF2.spect $C/$C.recombination_map.txt $C/header.$C.SF2.$N.out") > $C/sbatch.$C.$N.sweepfinder.sh;
      
      echo -e "sbatch $C/sbatch.$C.recomb.sh";
      echo -e "sbatch $C/sbatch.$C.$N.sweepfinder.norecmap.sh";
      echo -e "sbatch $C/sbatch.$C.$N.sweepfinder.sh";
    done;
#  fi;
done
