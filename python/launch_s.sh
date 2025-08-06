#!/usr/bin/env bash

for beta_r in 0.5 1.0 5.0 10.0 
do
    for s1 in 0.75 1.00 1.25 1.50 1.75 2.00 2.25 2.50 2.75 3.00
    do
        for beta_l in 1.0
        do
            
            for e1 in 0.00 
            do

                    sed -i "s/beta_r=.*/beta_r=${beta_r}/g" final_data_s.sub
                    sed -i "s/s1=.*/s1=${s1}/g" final_data_s.sub
                    sed -i "s/beta_l=.*/beta_l=${beta_l}/g" final_data_s.sub
                    sed -i "s/e1=.*/e1=${e1}/g" final_data_s.sub
                    sbatch final_data_s.sub
            done
        
        done
    done    
done
