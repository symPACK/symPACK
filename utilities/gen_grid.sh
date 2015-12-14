k=$1; s=$2; ./torus $k $k $s > ${k}x${k}_$s.adj; python adj_to_hb.py ${k}x${k}_$s.adj > ${k}x${k}_$s.rb
