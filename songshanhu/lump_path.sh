for file in tpt_temp.log*

do
echo ${file}
cat ${file} | grep '\[' > input

cat input | awk -F '\[' '{print $2}' | awk -F '\]' '{print $1}' > temp_paths
cat input | awk '{print $2}' > temp_flux
python lump_the_tpt_paths.py | sed 's#\]##g' > new_paths
cat new_paths | sort -s | uniq > tt

while read line
do
  echo $line;
  paste temp_flux new_paths | awk -v string=${line} 'BEGIN{sum=0;}{if($2==string) sum+=$1}END{print sum}'
done < tt > mm
cat mm | awk '{if(NR%2==1) printf "path-[%s ], pop:",$1;else print $1}'

done > tpt_pathlump.log
