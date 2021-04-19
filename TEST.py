##Créer fasta 1000 G: 
scp /home/saraoukkal/Documents/Stage_M1/G1000.fa soukkal@clus7.igfl.ens-lyon.fr:/media/disk1/soukkal/StageM2/Stage_M1/test/


##Convertir fasta en bed: 
python3 Genome_into_bases.py -fa /media/disk1/soukkal/StageM2/Stage_M1/test/G1000.fa -out /media/disk1/soukkal/StageM2/Stage_M1/test/G1000.txt

##Dupliquer la 3eme colonne:

##Trier le fichier:
sort -k1,1 -k2,2n /media/disk1/soukkal/StageM2/Stage_M1/test/G1000.txt
 > /media/disk1/soukkal/StageM2/Stage_M1/test/G1000_sorted.txt


##Créer fichier inter barrières:

##Positionner autour des barrières:
python3 AB_bar_dist.py -bar /media/disk1/soukkal/StageM2/Stage_M1/test/test_bar.txt -AB /media/disk1/soukkal/StageM2/Stage_M1/test/G1000_sorted.txt  -out /media/disk1/soukkal/StageM2/Stage_M1/test/G1000_count_bar.txt

##Faire plots:
python3 plot_sites.py -AB /media/disk1/soukkal/StageM2/Stage_M1/test/G1000_count_bar.txt
 -out /media/disk1/soukkal/StageM2/Stage_M1/test/

