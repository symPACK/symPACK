k=$1; s=$2; ./3d_adjacency $k $k $k $s > ${k}x${k}x${k}_$s.adj; python adj_to_hb.py ${k}x${k}x${k}_$s.adj > ${k}x${k}x${k}_$s.rb
