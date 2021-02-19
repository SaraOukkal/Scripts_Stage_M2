
#!/usr/bin/python3

#Importer packages: 
import argparse 


def load_bar(input_bar):
	"""
	Charge le fichier .bed des barrières dans un dictionnaire
	"""
	
	barriers ={} 
	f=open(input_bar,"r")
	bars=f.readlines()	
	
	for l in bars: 
		line=l.strip().split("\t")
	
		chrom=line[0]

		if chrom not in barriers.keys():
			barriers[chrom]=[]
		
		bar={}
		bar["st"]=int(line[1])
		bar["end"]=int(line[2])
			
		barriers[chrom].append(bar)
		
	
		
	return barriers #dictionnaire de chromosomes, chaque chromosome est une liste contenant des dictionnaires pour chaque barrière, un dictionnaire de barrière contient les positions start et end
	
	
	
	
def ancestral_base_bar(barriers, input_AB):
	"""
	Charge le fichier des bases non mutées
	"""
	dico={}	
	AB=open(input_AB,"r")
	done_chrom=[] 
	out_bar_bases=0 #initialise le compteur de bases non prises en compte 
	c = 0 #compteur de bases
	
	for l in AB: 
		if c%100000 == 0 : 
			print(c) #affiche C toutes les 100000 nt 
		line=l.strip().split("\t")
		chrom=line[0] #chromosome du nt 
		pos=int(line[1]) #position du nt
		nuc_C=line[2] #nucléotide chez le chimpanzé
		
		threshold=barriers[chrom][0]["st"]-1000 #limite basse au début du chromosome (sert à gagner du temps en début de chromosome)
		
		if chrom not in done_chrom:
			done_chrom.append(chrom)
			print(chrom) #imprime le nouveau chromosome pris en charge 
			index=0 #réinitialise l'index 
			
		for i in range(index,len(barriers[chrom])-1,1):
			st=barriers[chrom][i]["st"]	#start de la barrière x 
			end=barriers[chrom][i]["end"]	#end de la barrière x 
			st2=barriers[chrom][i+1]["st"]	#start de la barrière x+1 (suivante) 
			base=nuc_C.upper()	#détermine le type de base	
			
			if pos < threshold: #Pour le cas ou il y a des bases au début du chromosome avant la première barrière, trop loin pour qu'elles soient comptabilisées, cela permet de passer rapidement ces mutations et ne pas parcourir l'entiereté des barrières du chromosome inutilement
				out_bar_bases+=1
				break

			elif pos >=st and pos <=end : #si la base est dans la barrière
				dist1=pos-st
				dist2=end-pos		
				if dist1<=dist2: #choisis la distance la plus courte entre celle vers le start et celle vers le end 
					dist=-dist1		
				else:
					dist=-dist2
				if dist >= -50:
					if dist not in dico.keys(): #Si cette distance n'a pas encore été croisée on l'ajoute au dictionnaire 
						dico[dist]={}	
						
						dico[dist]["A"]=0	#Mets les compteurs de tous les types de nucléotides à 0 dans le cas d'une nouvelle distance 
						dico[dist]["C"]=0
						dico[dist]["G"]=0
						dico[dist]["T"]=0
						
					
					dico[dist][base]+=1 #Ajoute 1 au type de base concerné
					index=i #mets à jour l'index 
					break
					
			elif pos >= st-1000 and pos < st : #Si la base est  à - de 1000 nucléotides à gauche de la barrière
				dist=st-pos
				if dist >= -50:
					if dist not in dico.keys():
						dico[dist]={}	
						
						dico[dist]["A"]=0	#Mets les compteurs de tous les types de nucléotides à 0 dans le cas d'une nouvelle distance 
						dico[dist]["C"]=0
						dico[dist]["G"]=0
						dico[dist]["T"]=0			
						
					dico[dist][base]+=1
					index=i
					break
					
			elif pos <= end+1000 and pos > end: #Si la base est à - de 1000 nucléotides à droite de la barrière
				dist=pos-end
				if dist >= -50:
					if dist not in dico.keys():
						dico[dist]={}	
						
						dico[dist]["A"]=0	#Mets les compteurs de tous les types de nucléotides à 0 dans le cas d'une nouvelle distance 
						dico[dist]["C"]=0
						dico[dist]["G"]=0
						dico[dist]["T"]=0
						
					dico[dist][base]+=1
					index=i
					break
			
			elif pos < st2-1000 and pos > end+1000: #Si la base est entre deux barrières mais à + de 1000 nucléotides des deux, on ne la comptabilise pas mais on mets à jour l'index et on sort de la boucle
				index=i
				out_bar_bases+=1
				break
				
				
		c += 1 #mets à jour le compteur de bases totales 
					
	print("bases non concernées : ", out_bar_bases)
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
	parser.add_argument('-bar', '--input_bar', type=str, help='Path to barriers intervals', default ="/home/soukkal/Bureau/Projet/Step2_results/C_barriers_in_inter_sorted.bed")
	##fichier des mutations du chimpanzé :					
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to mutation positions and nuc ancestral state', default ="/home/soukkal/Bureau/Projet/Step3_results/AB_chimp.txt")			

	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step3_results/AB_count_bar.txt")			

	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	print("loading barriers file")
	barriers=load_bar(args.input_bar)
	print("counting bases around barriers")
	base_count=ancestral_base_bar(barriers, args.input_AB)
	print("writing counts in the output file")
	write_out(base_count,args.output)
	print("done")
	

if "__main__" == __name__:
	main()
