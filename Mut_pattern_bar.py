#!/usr/bin/python3

#Importer packages: 
import argparse 


def load_bar(input_bar):
	"""
	Charge le fichier .bed des barrières dans un dictionnaire
	"""
	
	barriers ={} 
	f=open(input_bar, "r")
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
	
	
	
	
def mut_bar(barriers, input_mut):
	dico={}	
	mut=open(input_mut,"r")
	done_chrom=[] 
	out_bar_mut=0 #initialise le compteur de mutations non prises en compte 
	c = 0 #compteur de mutations
	
	for l in mut: 
		if c%100000 == 0 : 
			print(c) #affiche C toutes les 100000 mutations 
		line=l.strip().split("\t")
		chrom=line[0] #chromosome de la mutation 
		pos=int(line[1]) #position de la mutation 
		nuc_C=line[2] #nucléotide chez le chimpanzé
		EA=line[3] #nucléotide ancestral 
		
		threshold=barriers[chrom][0]["st"]-1000 #limite basse au début du chromosome (sert à gagner du temps en début de chromosome)
		
		if chrom not in done_chrom:
			done_chrom.append(chrom)
			#print(chrom) #imprime le nouveau chromosome pris en charge 
			index=0 #réinitialise l'index 
			
		for i in range(index,len(barriers[chrom])-1,1):
			st=barriers[chrom][i]["st"]	#start de la barrière x 
			end=barriers[chrom][i]["end"]	#end de la barrière x 
			st2=barriers[chrom][i+1]["st"]	#start de la barrière x+1 (suivante) 
			mutation=nuc_C.upper()+">"+EA.upper()	#détermine le type de mutation 	
			
			if pos < threshold: #Pour le cas ou il y a des mutations au début du chromosome avant la première barrière, trop loin pour qu'elles soient comptabilisées, cela permet de passer rapidement ces mutations et ne pas parcourir l'entiereté des barrières du chromosome inutilement
				out_bar_mut+=1
				break

			elif pos >=st and pos <=end : #si la mutation est dans la barrière
				dist1=pos-st
				dist2=end-pos		
				if dist1<=dist2: #choisis la distance la plus courte entre celle vers le start et celle vers le end 
					dist=-dist1		
				else:
					dist=-dist2
				
				if dist >= -50: #Prend uniquement jusqu'a -50nt dans la barrière 
					if dist not in dico.keys(): #Si cette distance n'a pas encore été croisée on l'ajoute au dictionnaire 
						dico[dist]={}	
						
						dico[dist]["A>T"]=0	#Mets les compteurs de tous les types de mutations à 0 dans le cas d'une nouvelle distance 
						dico[dist]["A>C"]=0
						dico[dist]["A>G"]=0
						dico[dist]["C>T"]=0
						dico[dist]["C>A"]=0
						dico[dist]["C>G"]=0
						dico[dist]["G>T"]=0
						dico[dist]["G>A"]=0
						dico[dist]["G>C"]=0
						dico[dist]["T>A"]=0
						dico[dist]["T>C"]=0
						dico[dist]["T>G"]=0
					
					dico[dist][mutation]+=1 #Ajoute 1 au type de mutation concerné 
					index=i #mets à jour l'index 
					break
					
			elif pos >= st-1000 and pos < st : #Si la mutation est  à - de 1000 nucléotides à gauche de la barrière
				dist=st-pos
				if dist >= -50:
					if dist not in dico.keys():
						dico[dist]={}	
						
						dico[dist]["A>T"]=0
						dico[dist]["A>C"]=0
						dico[dist]["A>G"]=0
						dico[dist]["C>T"]=0
						dico[dist]["C>A"]=0
						dico[dist]["C>G"]=0
						dico[dist]["G>T"]=0
						dico[dist]["G>A"]=0
						dico[dist]["G>C"]=0
						dico[dist]["T>A"]=0
						dico[dist]["T>C"]=0
						dico[dist]["T>G"]=0
						
					dico[dist][mutation]+=1
					index=i
					break
					
			elif pos <= end+1000 and pos > end: #Si la mutation est à - de 1000 nucléotides à droite de la barrière
				dist=pos-end			
				if dist >= -50:
					if dist not in dico.keys():
						dico[dist]={}	
						
						dico[dist]["A>T"]=0
						dico[dist]["A>C"]=0
						dico[dist]["A>G"]=0
						dico[dist]["C>T"]=0
						dico[dist]["C>A"]=0
						dico[dist]["C>G"]=0
						dico[dist]["G>T"]=0
						dico[dist]["G>A"]=0
						dico[dist]["G>C"]=0
						dico[dist]["T>A"]=0
						dico[dist]["T>C"]=0
						dico[dist]["T>G"]=0
						
					dico[dist][mutation]+=1
					index=i
					break
			
			elif pos < st2-1000 and pos > end+1000: #Si la mutation est entre deux barrières mais à + de 1000 nucléotides des deux, on ne la comptabilise pas mais on mets à jour l'index et on sort de la boucle
				index=i
				out_bar_mut+=1
				break
				
			
				
				
		c += 1 #mets à jour le compteur de mutations totales 
					
	print("Mutations non concernées : ", out_bar_mut)
	return dico								


		
						
def write_out(mut_count,output_file):	
	"""
	Ecrit dans le fichier output le compte de chaque type de mutations par position par rapport aux barrières.
	"""
	
	out=open(output_file,"w")
	
	out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("dist","A>T","A>C","A>G","C>T","C>A","C>G","G>T","G>A","G>C","T>A","T>C","T>G"))#Ecrit un header au fichier 
	
	for dist in sorted(mut_count.keys()):
		out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(str(dist),str(mut_count[dist]["A>T"]),str(mut_count[dist]["A>C"]),str(mut_count[dist]["A>G"]),str(mut_count[dist]["C>T"]),str(mut_count[dist]["C>A"]),str(mut_count[dist]["C>G"]),str(mut_count[dist]["G>T"]),str(mut_count[dist]["G>A"]),str(mut_count[dist]["G>C"]),str(mut_count[dist]["T>A"]),str(mut_count[dist]["T>C"]),str(mut_count[dist]["T>G"])))


				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	##fichier .bed des barrières chez le chimpanzé: 
	parser.add_argument('-bar', '--input_bar', type=str, help='Path to barriers intervals', default ="/home/soukkal/Bureau/Projet/Step2_results/C_barriers_in_inter_sorted.bed")
	##fichier des mutations du chimpanzé :					
	parser.add_argument('-mut', '--input_mut', type=str, help='Path to mutation positions and nuc ancestral state', default ="/home/soukkal/Bureau/Projet/Step3_results/mut_EA_sorted.txt")			

	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step3_results/mut_count_bar.txt")			

	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	print("loading barriers file")
	barriers=load_bar(args.input_bar)
	print("counting mutations around barriers")
	mut_count=mut_bar(barriers, args.input_mut)
	print("writing counts in the output file")
	write_out(mut_count,args.output)
	print("done")
	

if "__main__" == __name__:
	main()
