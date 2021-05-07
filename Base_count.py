#!/usr/bin/python3

#Importer packages: 
import argparse 

def base_count(input_base):
	"""
	Charge le fichier des bases et les compte en fonction du type 
	"""

	bases=open(input_base,"r")
	done_chrom=[] 
	A_count=0
	T_count=0
	C_count=0
	G_count=0
	
	for l in bases:
		line=l.strip().split("\t")
		chrom=line[0] #chromosome de la mutation 
		pos=int(line[1]) #position de la mutation 
		nuc=line[2] #nucléotide 

		if chrom not in done_chrom:
			done_chrom.append(chrom)
			print(chrom) #imprime le nouveau chromosome pris en charge 
		
		if nuc == "A":
			A_count+=1
		elif nuc == "T":
			T_count+=1
		elif nuc == "C":
			C_count+=1
		elif nuc == "G": 
			G_count+=1
	
	print("A", A_count)
	print("T", T_count)
	print("C", C_count)
	print("G", G_count) 

										
	
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	parser.add_argument('-in', '--input', type=str, help='Path to homo poly nucleotides', default ="/home/soukkal/Bureau/Projet/Step3_results/mut_EA_sorted.txt")			
	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	base_count(args.input)

if "__main__" == __name__:
	main()
