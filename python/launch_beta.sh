#!/usr/bin/env bash

for beta_r in 0.5 1.0 5.0 10.0
do
    for beta_l in 0.100 0.15199111 0.23101297 0.35111917 0.53366992 0.81113083 1.23284674 1.87381742 2.84803587 4.32876128 6.57933225 10.000
    do
        for g1 in 0.016
        do
            for ham_type in 1 2 
            do
                for e1 in 0.00 0.01
                do

                    sed -i "s/beta_r=.*/beta_r=${beta_r}/g" final_data_beta.sub
                    sed -i "s/beta_l=.*/beta_l=${beta_l}/g" final_data_beta.sub
                    sed -i "s/ham_type=.*/ham_type=${ham_type}/g" final_data_beta.sub
                    sed -i "s/g1=.*/g1=${g1}/g" final_data_beta.sub
                    sed -i "s/e1=.*/e1=${e1}/g" final_data_beta.sub
                    sbatch final_data_beta.sub
                done
            done
        done
    done
done
