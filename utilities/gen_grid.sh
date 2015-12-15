k=$1; s=$2; python gen_grid.py $k $k $s > ${k}x${k}_$s.adj; python adj_to_hb.py ${k}x${k}_$s.adj > ${k}x${k}_$s.rb
