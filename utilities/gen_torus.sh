k=$1; s=$2; ./torus $k $k $s > tor_${k}x${k}_$s.adj; python adj_to_hb.py tor_${k}x${k}_$s.adj > tor_${k}x${k}_$s.rb
