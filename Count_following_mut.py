#!/usr/bin/python3

#Importer packages: 
import argparse 

def mut_count(input_mut):
	"""
	Charge le fichier des mutations ponctuelles et compte les mutations successives
	"""

	mutations=open(input_mut,"r")
	done_chrom=[] 
	mut=mutations.readlines()
	fol_count=0
	
	for i in range(len(mut)-1): 
		line=mut[i].strip().split("\t")
		chrom=line[0] #chromosome de la mutation 
		pos=int(line[1]) #position de la mutation 
		nuc_C=line[2] #nucléotide chez l'espèce d'interêt
		EA=line[3] #nucléotide ancestral 
		line2=mut[i-1].strip().split("\t")
		previous_pos=int(line2[1])
		line3=mut[i+1].strip().split("\t")
		next_pos=int(line3[1])
		
		if chrom not in done_chrom:
			done_chrom.append(chrom)
			print(chrom) #imprime le nouveau chromosome pris en charge 
		
		if previous_pos == pos-1:
			fol_count+=1
		elif next_pos == pos+1:
			fol_count+=1
	
	print("Mutations following each other", fol_count)

										
	
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	parser.add_argument('-mut', '--input_mut', type=str, help='Path to mutation positions and nuc ancestral state', default ="/home/soukkal/Bureau/Projet/Step3_results/mut_EA_sorted.txt")			
	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	mut_count(args.input_mut)

if "__main__" == __name__:
	main()
