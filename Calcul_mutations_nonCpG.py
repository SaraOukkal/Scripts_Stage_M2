#!/usr/bin/python3

#Importer packages: 
import argparse 


def calcul_mut(input_AB, nonCpG, CpG): 
	data=open(input_AB,"r") #Fichier des bases ancestrales + bases actuelles
	nonCpG=open(nonCpG,"w") #Fichier de sortie mutations non CpG
	CpG=open(CpG,"w") #Fichier de sortie mutations CpG 
	done_chrom=[]
	c=0
	
	AB=data.readline()
	
	for l in AB:
		c+=1
	
	for i in range(c): 
		line=AB[i].strip().split("\t")
		chrom=line[0] #chromosome du nt 
		pos=int(line[1]) #position du nt
		nuc=line[2] #Nucléotide de l'espèce d'interêt 
		nuc_EA=line[3] #nucléotide à l'état ancestral
		
		line2=AB[i+1].strip().split("\t") #Nucléotide suivant pour repérer les CpG 
		pos2=int(line2[1])
		nuc2=line2[2]

		if chrom not in done_chrom:
			done_chrom.append(chrom)
			print(chrom)
		
		if nuc != nuc_EA: #S'il y a une mutation
			if nuc == "C": #Si c'est un C
				if pos2 == pos+1: #Si la base suivante du fichier la suit dans la séquence 
					if nuc2 == "G": #Si c'est un CpG 
						CpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, nuc, nuc_EA)) 
					else: #Si ce n'est pas un CpG
						nonCpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, nuc, nuc_EA)) 
				
				else: #Si la base suivante du fichier ne la suit pas dans la séquence 
					continue #Passe directement à l'itération suivante
			else: #Si ce n'est pas un C
				nonCpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, nuc, nuc_EA)) 
			

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to ancestral bases', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")						
	
	#fichier de sortie: 
	parser.add_argument('-nonCpG', '--nonCpG', type=str, help='Path to non CpG output file', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")			
	parser.add_argument('-CpG', '--CpG', type=str, help='Path to CpG output file', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")		
	
	args = parser.parse_args()
	
	calcul_mut(args.input_AB, args.nonCpG, args.CpG)
	
	
if "__main__" == __name__:
	main()
