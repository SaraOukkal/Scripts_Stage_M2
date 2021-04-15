#!/usr/bin/python3

#Importer packages: 
import argparse 
	
def load_AB(input_AB, nonCpG, CpG):
	"""
	Charge les bases dans un dictionnaire par chromosome et lance remove_CpG
	"""
	
	data=open(input_AB,"r") #Fichier des bases ancestrales + bases actuelles
	nonCpG=open(nonCpG,"w") #Fichier de sortie mutations non CpG
	CpG=open(CpG,"w") #Fichier de sortie mutations CpG 
	done_chrom=[]
	
	for l in data:
		line=l.strip().split("\t")
		chrom=line[0] #chromosome du nt 
		pos=int(line[1]) #position du nt
		nuc=line[2] #Nucléotide actuel
		EA=line[3] #nucléotide à l'état ancestral
		
		base=[]	
		base.append(chrom)
		base.append(pos)
		base.append(nuc)
		base.append(EA)
		
		if chrom not in done_chrom:
			if chrom != "chr1" :
				print("start calculating mutations")
				remove_CpG(chromosome, nonCpG, CpG)
				print("done calculating mutations")
			chromosome=[]
			print(chrom)
			done_chrom.append(chrom)
		
		chromosome.append(base)
		
	print("start calculating last mutations")
	calcul_mut(chromosome, nonCpG, CpG)
	print("done calculating last mutations")		
		
			
	
def remove_CpG(chromosome, nonCpG, CpG):
	"""
	Trie les bases ancestrales en CpG et nonCpG et les écrit dans deux fichiers
	"""
	
	for i in range(len(chromosome)):
		chrom=chromosome[i][0]
		pos=int(chromosome[i][1])
		nuc=chromosome[i][2]
		EA=chromosome[i][3]
			
		if i == len(chromosome)-1: #Si c'est la fin du chromosome
			break
			
		pos2=int(chromosome[i+1][1]) #Base qui suit dans le fichier 
		EA2=chromosome[i+1][3]
		
		if EA == "C": #Si c'est un C
			if pos2 == pos+1: #Si la base suivante du fichier la suit dans la séquence 
				if EA2 == "G": #Si c'est un CpG 
					#print("CpG")
					CpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, nuc, EA)) 
			
				else: #Si ce n'est pas un CpG
					nonCpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, nuc, EA)) 
		
			else: #Si la base suivante de la séquence n'est pas une base ancestrale
				continue #Passe directement à l'itération suivante
					
		else: #Si ce n'est pas un C 
			nonCpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, nuc, EA)) 
		
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to ancestral bases', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")						
	
	#fichier de sortie: 
	parser.add_argument('-nonCpG', '--nonCpG', type=str, help='Path to non CpG output file', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")			
	parser.add_argument('-CpG', '--CpG', type=str, help='Path to CpG output file', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")		
	
	args = parser.parse_args()
	print("start loading AB")
	load_AB(args.input_AB, args.nonCpG, args.CpG)

	
	
if "__main__" == __name__:
	main()
