#phonopy -f DIR_POSCAR*/vasprun.xml
phonopy -p band.conf -s
bandplot --gnuplot band.yaml >bands.gp
phonopy -t -s -p  mesh.conf
propplot --gnuplot thermal_properties.yaml>thermal_properties.gp
phonopy -s -p pdos.conf
phonopy -s -p dos.conf


sed 's/ANIME =/ANIME = 0.0 0.0 0.0/g' mod.conf>mod_tmp.conf
phonopy -t mod_tmp.conf
cp anime.ascii anime_tot.ascii
sed 's/ANIME =/ANIME = 0.5 0.0 0.0/g' mod.conf>mod_tmp.conf
phonopy -t mod_tmp.conf
grep "#" anime.ascii>>anime_tot.ascii
sed 's/ANIME =/ANIME = 0.0 0.5 0.0/g' mod.conf>mod_tmp.conf
phonopy -t mod_tmp.conf
grep "#" anime.ascii>>anime_tot.ascii
sed 's/ANIME =/ANIME = 0.0 0.0 0.5/g' mod.conf>mod_tmp.conf
phonopy -t mod_tmp.conf
grep "#" anime.ascii>>anime_tot.ascii
sed 's/ANIME =/ANIME = 0.5 0.5 0.0/g' mod.conf>mod_tmp.conf
phonopy -t mod_tmp.conf
grep "#" anime.ascii>>anime_tot.ascii
sed 's/ANIME =/ANIME = 0.5 0.0 0.5/g' mod.conf>mod_tmp.conf
phonopy -t mod_tmp.conf
grep "#" anime.ascii>>anime_tot.ascii
sed 's/ANIME =/ANIME = 0.0 0.5 0.5/g' mod.conf>mod_tmp.conf
phonopy -t mod_tmp.conf
grep "#" anime.ascii>>anime_tot.ascii
sed 's/ANIME =/ANIME = 0.5 0.5 0.5/g' mod.conf>mod_tmp.conf
phonopy -t mod_tmp.conf
grep "#" anime.ascii>>anime_tot.ascii

