#!/usr/bin/python3

#Importer packages: 
import argparse 
	
def load_AB(input_AB):
	data=open(input_AB,"r") #Fichier des bases ancestrales + bases actuelles
	
	done_chrom=[]
	AB={}
	
	for l in data:
		base=[]		
		
		line=l.strip().split("\t")
		chrom=line[0] #chromosome du nt 
		pos=int(line[1]) #position du nt
		nuc=line[2] #Nucléotide actuel
		EA=line[3] #nucléotide à l'état ancestral

		base.append(pos)
		base.append(nuc)
		base.append(EA)
		
		if chrom not in done_chrom:
			print(chrom)
			done_chrom.append(chrom)
			AB[chrom]=[]
		
		AB[chrom].append(base)
		
	return AB
			
	
def calcul_mut(AB, nonCpG, CpG):
	nonCpG=open(nonCpG,"w") #Fichier de sortie mutations non CpG
	CpG=open(CpG,"w") #Fichier de sortie mutations CpG 
		
	for chrom in AB.keys():
		for i in AB[chrom]:
			pos=int(AB[chrom][i][0])
			nuc=AB[chrom][i][1]
			EA=AB[chrom][i][2]
			pos2=int(AB[chrom][i+1][0])
			nuc2=AB[chrom][i+1][1]
			
			if nuc != EA: #si il y a une mutation
				if nuc == "C": #Si c'est un C
					if pos2 == pos+1: #Si la base suivante du fichier la suit dans la séquence 
						if nuc2 == "G": #Si c'est un CpG 
							print("CpG")
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
	AB=load_AB(args.input_AB)
	print("load AB done")
	print("start calculating mutations")
	calcul_mut(AB, args.nonCpG, args.CpG)
	print("done calculating mutations")
	
	
if "__main__" == __name__:
	main()
