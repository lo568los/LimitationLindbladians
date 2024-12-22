#!/usr/bin/env bash

for beta_r in 0.5 1 5 10 
do
    for g1 in 0.01 0.01873817 0.03511192 0.06579332 0.12328467 0.23101297 0.43287613 0.81113083 1.51991108 2.84803587 5.33669923 10.0
    do
        for beta_l in 1
        do
            for ham_type in 1 2  
            do
                for e1 in 0 0.01
                do

                    sed -i "s/beta_r=.*/beta_r=${beta_r}/g" final_data_g.sub
                    sed -i "s/g1=.*/g1=${g1}/g" final_data_g.sub
                    sed -i "s/ham_type=.*/ham_type=${ham_type}/g" final_data_g.sub
                    sed -i "s/beta_l=.*/beta_l=${beta_l}/g" final_data_g.sub
                    sed -i "s/e1=.*/e1=${e1}/g" final_data_g.sub
                    sbatch final_data_g.sub
                done
            done
        done
    done    
done
