
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
	
		chrom=line[0]

		if chrom not in barriers.keys():
			barriers[chrom]=[]
		
		bar={}
		bar["st1"]=int(line[1])
		bar["end1"]=int(line[2])
		bar["st2"]=int(line[3])
		bar["end2"]=int(line[4])

		barriers[chrom].append(bar)
		
	
		
	return barriers #dictionnaire de chromosomes, chaque chromosome est une liste contenant des dictionnaires pour chaque barrière, un dictionnaire de barrière contient les positions start et end
	
	
	
	
def mut_bar(barriers, input_mut):
	"""
	Charge le fichier des mutations ponctuelles et les place en fonction des barrières
	"""
	dico={}	
	mut=open(input_mut,"r")
	done_chrom=[] 
	c = 0 #compteur de mutations
	
	for l in mut: 
		if c%100000 == 0 : 
			print(c) #affiche C toutes les 100000 mutations 
		line=l.strip().split("\t")
		chrom=line[0] #chromosome de la mutation 
		pos=int(line[1]) #position de la mutation 
		nuc_C=line[2] #nucléotide chez le chimpanzé
		EA=line[3] #nucléotide ancestral 
		
		threshold=barriers[chrom][0]["st"]-1000 #limite basse au début du chromosome donc ce qu'il y a 1000nt avant la première barrière (sert à gagner du temps en début de chromosome)
		
		if chrom not in done_chrom:
			done_chrom.append(chrom)
			#print(chrom) #imprime le nouveau chromosome pris en charge 
			index=0 #réinitialise l'index 
			
		for i in range(index,len(barriers[chrom])-1,1):
			st1=barriers[chrom][i]["st1"]	#start de la barrière 1
			end1=barriers[chrom][i]["end1"]	#end de la barrière 1 
			st2=barriers[chrom][i]["st2"]
			end2=barriers[chrom][i]["end2"]
			mutation=nuc_C.upper()+">"+EA.upper()	#détermine le type de mutation 	
			
			if pos < threshold: #Pour le cas ou il y a des mutations au début du chromosome avant la première barrière, trop loin pour qu'elles soient comptabilisées, cela permet de passer rapidement ces mutations et ne pas parcourir l'entiereté des barrières du chromosome 
				break

			elif pos >=st1 and pos <=end1 : #si la mutation est dans la première barrière
				dist1=st1-pos
				dist2=pos-end1		
				dist=max(dist1,dist2)
				if dist >= -50: #Prend uniquement jusqu'a -50nt dans la barrière 
					if dist not in dico.keys(): #Si cette distance n'a pas encore été croisée on l'ajoute au dictionnaire 
						dico[dist]=Counter()	
					dico[dist][mutation]+=1 #Ajoute 1 au type de mutation concerné 
					index=i #mets à jour l'index 
					break		
			
			elif pos >=st2 and pos <=end2: #Si la mutation est dans la deuxième barrière
				dist1=st2-pos
				dist2=pos-end2		
				dist=max(dist1,dist2)
				if dist >= -50:
					if dist not in dico.keys(): #Si cette distance n'a pas encore été croisée on l'ajoute au dictionnaire 
						dico[dist]=Counter()	
					dico[dist][base]+=1 #Ajoute 1 au type de base concerné
					index=i #mets à jour l'index 
					break			
					
			elif pos <= end1+1000 and pos > end1: #Si la mutation est à - de 1000 nucléotides après la première barrière
				dist=pos-end			
				if dist not in dico.keys():
					dico[dist]=Counter()	
				dico[dist][mutation]+=1
				index=i
				break

			elif pos >= st2-1000 and pos < st2 : #Si la mutation est  à - de 1000 nucléotides avant la deuxième barrière
				dist=st2-pos
				if dist not in dico.keys():
					dico[dist]=Counter()	
				dico[dist][mutation]+=1
				index=i
				break

			elif pos < st2-1000 and pos > end+1000: #Si la mutation est entre les deux barrières mais à + de 1000 nucléotides des deux, on ne la comptabilise pas mais on mets à jour l'index et on sort de la boucle
				index=i
				break
				
		c += 1 #mets à jour le compteur de mutations totales 

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
	parser.add_argument('-bar', '--input_bar', type=str, help='Path to barriers intervals', default ="/media/disk1/soukkal/StageM2/Stage_M1/Chimp_Step2_results/selected_inter_bar_sorted.be")
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
