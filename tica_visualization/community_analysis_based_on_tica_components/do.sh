#step1:
###get the tica contribution for each atom pair
paste ../../important_pair_cac_star0.name ../eigvector_corr80.txt| awk '{print $1,$2,sqrt($11*$11)}' > atom1_atom2_from1_tic1_square.contribution

##step2:extract the tica importance matrix


##step3: hierarchical clustering


##step 4: optimization
