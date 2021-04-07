#!/usr/bin/python3

#Importer packages: 
import argparse 


def calcul_mut(input_AB, output): 
	f=open(input_AB,"r")
	out=open(output,"w")
	AB=f.readlines()
	done_chrom=[]
	
	for l in AB: 
		line=l.strip().split("\t")
		chrom=line[0] #chromosome du nt 
		pos=int(line[1]) #position du nt
		nuc=line[2] #Nucléotide de l'espèce d'interêt 
		nuc_EA=line[3] #nucléotide à l'état ancestral
		mut_count=0

		if chrom not in done_chrom:
			done_chrom.append(chrom)
			print(chrom)
		
		if nuc != nuc_EA: 
			print('yay')
			mut_count+=1
			out.write("{}\n".format(l)) #Ecrit la ligne entière 
	
	print(mut_count)

def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-AB', '--input_AB', type=str, help='Path to ancestral bases', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")						
	
	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")			

	
	args = parser.parse_args()
	
	calcul_mut(args.input_AB, args.output)
	
	
if "__main__" == __name__:
	main()
