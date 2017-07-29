paste kcenters_none_stationary_population kcenters_transpose_stationary_population kcenters_mle_stationary_population | awk -F '+' '{print $1,$3}' | awk '{print $1,$3,$4}' | awk -F '(' '{print $2}' >temp

#gracebat -block temp -bxy 0:1 -block temp -bxy 0:2 -block temp -bxy 0:3 -a.xvg 
sort -g -k3 temp>anquan
gracebat -block anquan -bxy 0:1 -block anquan -bxy 0:2 -block anquan -bxy 0:3 a.xvg -hdevice PNG -printfile none_transpose_mle_stat_pop.png

cat -n temp | awk '{print $1,$2-$4}' | awk '{if($2<0) print $1,-$2;else print $0}' | sort -g -k2 -r | head
