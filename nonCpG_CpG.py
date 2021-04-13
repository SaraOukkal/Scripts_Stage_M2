#!/usr/bin/python3

#Importer packages: 
import argparse 
	
def load_AB(input_AB):
	data=open(input_AB,"r") #Fichier des bases ancestrales + bases actuelles
	
	done_chrom=[]
	AB={}
	
	for l in data:
		line=l.strip().split("\t")
		chrom=line[0] #chromosome du nt 
		pos=int(line[1]) #position du nt
		nuc=line[2] #Nucléotide actuel
		nuc_EA=line[3] #nucléotide à l'état ancestral
		
		if chrom not in done_chrom:
			print(chrom)
			done_chrom.append(chrom)
			AB[chrom]={}
			
		AB[chrom][pos]=[]
		AB[chrom][pos].append(nuc)
		AB[chrom][pos].append(nuc_EA)
		
	return AB
			
	
def calcul_mut(AB, nonCpG, CpG):
	nonCpG=open(nonCpG,"w") #Fichier de sortie mutations non CpG
	CpG=open(CpG,"w") #Fichier de sortie mutations CpG 
		
	for chrom in AB.keys():
		for pos in AB[chrom].keys():
			if AB[chrom][pos][0] != AB[chrom][pos][1]: #si nuc différent de nuc_EA y a une mutation
				if AB[chrom][pos][0] == "C": #Si c'est un C
					if pos+1 in AB[chrom].keys(): #Si la base suivante du fichier la suit dans la séquence 
						print(pos+1)
						if AB[chrom][pos+1][0] == "G": #Si c'est un CpG 
							print("CpG")
							CpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, AB[chrom][pos][0], AB[chrom][pos][1])) 
				
						else: #Si ce n'est pas un CpG
							print("nonCpG")
							nonCpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, AB[chrom][pos][0], AB[chrom][pos][1])) 
			
					else: #Si la base suivante de la séquence n'est pas une base ancestrale
						continue #Passe directement à l'itération suivante
						
				else: #Si ce n'est pas un C
					nonCpG.write("{}\t{}\t{}\t{}\n".format(chrom, pos, AB[chrom][pos][0], AB[chrom][pos][1])) 
			else:
				print("pas de mutation")
			
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
