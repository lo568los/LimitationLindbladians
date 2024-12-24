#!/usr/bin/env bash

for beta_r in 0.5 1.0 5.0 10.0 
do
    for g1 in 0.0010 0.00187382 0.00351119 0.00657933 0.01232847 0.0231013 0.04328761 0.08111308 0.15199111 0.28480359 0.53366992 1.0000
    do
        for beta_l in 1
        do
            for ham_type in 1 2  
            do
                for e1 in 0.00 0.01
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
