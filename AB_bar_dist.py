
#!/usr/bin/python3

#Importer packages: 
import argparse 
from collections import Counter


def load_bar(input_bar):
	"""
	Charge le fichier .bed des barrières dans un dictionnaire
	"""
	
	barriers ={} 
	f=open(input_bar,"r")
	bars=f.readlines()	
	
	for l in bars: 
		line=l.strip().split("\t")
	
		chrom=line[0] #Chromosome

		if chrom not in barriers.keys():
			barriers[chrom]=[]
		
		bar={}
		bar["st1"]=int(line[1]) #Début de la première barrière
		bar["end1"]=int(line[2]) #Fin de la première barrière 
		bar["st2"]=int(line[3]) #Début de la deuxième barrière
		bar["end2"]=int(line[4]) #Fin de la deuxième barrière
			
		barriers[chrom].append(bar) #Ajouter l'inter barrière au dictionnaire
		
	
		
	return barriers #dictionnaire de chromosomes, chaque chromosome est une liste contenant des dictionnaires pour chaque barrière, un dictionnaire de barrière contient les positions start et end
	
	
	
	
def base_bar(barriers, input_AB):
	"""
	Charge le fichier des bases non mutées et les place en fonction des barrières
	"""
	dico={}	
	AB=open(input_AB,"r")
	done_chrom=[] 
	c = 0 #compteur de bases

	for l in AB: 
		if c%100000 == 0 : 
			print(c) #affiche C toutes les 100000 nt 
		line=l.strip().split("\t")
		chrom=line[0] #chromosome du nt 
		pos=int(line[1]) #position du nt
		nuc=line[2] #Nucléotide à l'état actuel
		EA=line[3] #nucléotide à l'état ancestral
				
		if chrom not in done_chrom:
			done_chrom.append(chrom)
			print(chrom) #imprime le nouveau chromosome pris en charge 
			index=0 #réinitialise l'index 
			
		for i in range(index,len(barriers[chrom]),1):
			st1=barriers[chrom][i]["st1"]	#start de la barrière x 
			end1=barriers[chrom][i]["end1"]	#end de la barrière x 
			st2=barriers[chrom][i]["st2"]
			end2=barriers[chrom][i]["end2"]
			base=EA.upper()	#détermine le type de base ancestrale (fichier bases ancestrales)
			
			mid_bar1=end1-((end1-st1)/2)
			mid_bar2=st2+((end2-st2)/2)
			mid_inter_bar=end1+((st2-end1)/2)
			
			if mid_bar1 < pos and pos < mid_bar2:
				if pos <= end1 or pos <= mid_inter_bar:
					dist=pos-end1
				elif pos <=st2 or pos <= mid_bar2: 
					dist=st2-pos 
				
				if dist not in dico.keys(): #Si cette distance n'a pas encore été croisée on l'ajoute au dictionnaire 
					dico[dist]=Counter()	
				
				dico[dist][base]+=1	
				index=i #mets à jour l'index 
				
		#Si la base est après la deuxième barrière on continue de parcourir les barrières
				
		c += 1 #mets à jour le compteur de bases totales 
					
	return dico								


		
						
def write_out(base_count,output_file):	
	"""
	Ecrit dans le fichier output le compte de chaque type de bases par position par rapport aux barrières.
	"""
	
	out=open(output_file,"w")
	
	out.write("{}\t{}\t{}\t{}\t{}\n".format("dist","A","C","G","T"))#Ecrit un header au fichier 
	
	for dist in sorted(base_count.keys()):
		out.write("{}\t{}\t{}\t{}\t{}\n".format(str(dist),str(base_count[dist]["A"]),str(base_count[dist]["C"]),str(base_count[dist]["G"]),str(base_count[dist]["T"])))


				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	##fichier .bed des barrières chez le chimpanzé: 
	parser.add_argument('-bar', '--input_bar', type=str, help='Path to barriers intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Chimp_Step2_results/selected_inter_bar_sorted.be")
	##fichier des mutations du chimpanzé :					
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to mutation positions and nuc ancestral state', default ="/home/soukkal/Bureau/Projet/Step3_results/AB_chimp.txt")			

	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step3_results/AB_count_bar.txt")			

	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	print("loading barriers file")
	barriers=load_bar(args.input_bar)
	print("counting bases around barriers")
	base_count=base_bar(barriers, args.input_AB)
	print("writing counts in the output file")
	write_out(base_count,args.output)
	print("done")
	

if "__main__" == __name__:
	main()
