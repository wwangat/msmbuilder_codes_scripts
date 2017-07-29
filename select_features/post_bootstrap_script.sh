for tIC in `seq 1 20`;do a=1;while read -r line;do awk -v var=$tIC '{print $var}' $line>temp_${a};a=$(($a+1));done<for_statistics.txt ;paste temp_*>for_stat_tIC_${tIC};done;rm temp*
for((tIC=1;tIC<=19;tIC++));do paste bootstrap_1/pairwise_dist5  final_stat_tIC_${tIC}_min_max_mean_median_std | awk '{for(v=1;v<=NF;v++) printf "%s\t",$v;print '\n'}' | sort -g -k5 -r | awk '{print $5, $7}';done
cat final_stat_tIC_1_min_max_mean_median_std | sort -k3 -g >temp;xmgrace -block temp xy -bxy 0:4 -block temp -settype xydy -bxy 0:3:5
