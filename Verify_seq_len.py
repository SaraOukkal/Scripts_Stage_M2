#!/usr/bin/python3

#Importer packages: 
import argparse 


def verify_seq_length(input1,input2,input3,input4):
	"""
	Vérifie que les séquences font la même taille dans les fichiers des 4 espèces
	"""
	
	H=open(input1,"r")
	C=open(input2,"r")
	G=open(input3,"r")
	O=open(input4,"r")
	
	H_len=calculate_seq_length(H)
	print("finished calculating Human seq lengths")
	H.close()
	
	C_len=calculate_seq_length(C)
	print("finished calculating Chimpanzee seq lengths")
	C.close()
	
	G_len=calculate_seq_length(G)
	print("finished calculating Gorilla seq lengths")
	G.close()
	
	O_len=calculate_seq_length(O)
	print("finished calculating Orangutan seq lengths")
	O.close()
	
	if H_len == C_len and H_len == G_len and H_len == O_len:
		print("Sequences of 4 species are the same length")		
			
def calculate_seq_length(input_seq):
	"""
	Calcule la longueur de chaque séquence et l'ajoute à une liste 
	"""
	
	seq_len=[]
	for line in input_seq:
		if not line.startswith(">"):
			seq_len.append(len(line))
			
	return seq_len
							
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-H', '--input_H', type=str, help='Path to Human sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/H_seq.fa")					
	parser.add_argument('-C', '--input_C', type=str, help='Path to Chimp sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/C_seq.fa")			
	parser.add_argument('-G', '--input_G', type=str, help='Path to Gorilla sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/G_seq.fa")			
	parser.add_argument('-O', '--input_O', type=str, help='Path to Orangutan sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/O_seq.fa")			
			
	args = parser.parse_args()
	
	#Vérifie que les séquences font la même taille dans les 4 fichiers:
	verify_seq_length(args.input_H, args.input_C, args.input_G, args.input_O)


if "__main__" == __name__:
	main()
