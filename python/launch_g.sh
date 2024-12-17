#!/usr/bin/env bash

for beta_r in 0.5 1 5 10
do
    for beta_l in 1
    do
        for ham_type in 1 2 
        do
            for e in 0 0.01
            do

                sed -i "s/beta_r=.*/beta_r=${beta_r}/g" final_data_g.sub
                sed -i "s/ham_type=.*/ham_type=${ham_type}/g" final_data_g.sub
                sed -i "s/beta_l=.*/beta_l=${beta_l}/g" final_data_g.sub
                sed -i "s/e1=.*/e1=${e1}/g" final_data_g.sub
                sbatch final_data_g.sub
            done
        done
    done
done
