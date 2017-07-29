# in this script, *_2norm means square (probability)


cat tica_eigenvectors-4 | sort -n -k1>/home/group/leijp/work_July12_analyze_bootstra_ww/GS_NO_bootstrap_tIC1_nosort
b=`cat GS_NO_bootstrap_tIC3_nosort |  awk 'BEGIN{a=0} {a=a+$4*$4;} END{print a}' `;awk -v var=$b '{print $1,$2,$3,$4*$4/var}' GS_NO_bootstrap_tIC3_nosort>GS_NO_bootstrap_tIC3_nosort_2norm

#In /home/group/leijp/bootstrap_GS/probability
paste ../bootstrap_1/pairwise_dist5 final_stat_tIC_3_min_max_mean_median_std | cat -n | awk '{print $1, $2, $3, $7}'>/home/group/leijp/work_July12_analyze_bootstra_ww/GS_YES_bootstrap_tIC3_nosort
awk '{print $4}' GS_NO_bootstrap_tIC1_nosort_norm2 >temp0;paste ~/leijp/bootstrap_TS2/final_stat_tIC_1_min_max_mean_median_std temp0 | sort -k3 -g >temp;xmgrace -block temp xy -bxy 0:4 -block temp xy -bxy 0:6 -block temp -settype xydy -bxy 0:3:5



#In /home/group/leijp/bootstrap_TS2/eigen20
cat tica_eigenvectors-4 | sort -n -k1>/home/group/leijp/work_July12_analyze_bootstra_ww/TS2_NO_bootstrap_tIC1_nosort
b=`cat TS2_NO_bootstrap_tIC1_nosort |  awk 'BEGIN{a=0} {a=a+$4*$4;} END{print a}' `;awk -v var=$b '{print $1,$2,$3,$4*$4/var}' TS2_NO_bootstrap_tIC1_nosort>TS2_NO_bootstrap_tIC1_nosort_2norm

#In /home/group/leijp/bootstrap_TS2/min_max_mean_median_std
paste ../bootstrap_1/pairwise_dist5 final_stat_tIC_3_min_max_mean_median_std | cat -n | awk '{print $1, $2, $3, $7}'>/home/group/leijp/work_July12_analyze_bootstra_ww/GS_YES_bootstrap_tIC3_nosort



#drawing figures...........................

#set x-axis as 26450:26580
#For GS
#xmgrace to draw the figure, show the robustness of the bootstrapping results and compare with no-bootstrap results, here all the eigenvectors are a**2 to be the probability
rm temp;rm temp0;awk '{print $4}' GS_NO_bootstrap_tIC2_nosort_2norm >temp0;paste ~/leijp/bootstrap_GS/probability/final_stat_tIC_2_min_max_mean_median_std temp0 | sort -k3 -g >temp;xmgrace a.xvg -block temp xy -bxy 0:4 -block temp xy -bxy 0:6 -block temp -settype xydy -bxy 0:3:5

#For TS2
#xmgrace to draw the figure, show the robustness of the bootstrapping results and compare with no-bootstrap results, here all the eigenvectors are a**2 to be the probability
rm temp;rm temp0;awk '{print $4}' TS2_NO_bootstrap_tIC2_nosort_2norm >temp0;paste ~/leijp/bootstrap_TS2/min_max_mean_median_std/final_stat_tIC_2_min_max_mean_median_std temp0 | sort -k3 -g >temp;xmgrace a.xvg -block temp xy -bxy 0:4 -block temp xy -bxy 0:6 -block temp -settype xydy -bxy 0:3:5
