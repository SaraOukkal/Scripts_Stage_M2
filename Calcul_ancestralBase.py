#!/usr/bin/python3

#Importer packages: 
import argparse 



def same_base(fa_file1, fa_file2, fa_file3, fa_file4, output_file):
	"""
	Rapporte les nucléotides qui sont conservés entres les 4 espèces et leurs positions
	"""
		
	fini=False
	
	out=open(output_file,"w")
	
	H_file=open(fa_file1,"r")
	C_file=open(fa_file2,"r")
	G_file=open(fa_file3,"r")
	O_file=open(fa_file4,"r")
	
	done_chrom=[]
	nucleotides=["A","C","T","G"]
	
	while not fini: 
		H=H_file.readline()
		C=C_file.readline()
		G=G_file.readline()
		O=O_file.readline()
		
		if H=="":
			fini=True
		
		if H.startswith(">"): #Lignes >chr:st-end
			line_H=H.strip().split(":")
			chrom_H=line_H[0].replace(">","")  #Chromosome chez l'Humain 
			mut=[]
			
			if chrom_H not in done_chrom:
				done_chrom.append(chrom_H)
				print(chrom_H) #Sert à savoir ou on en est dans le fichier (au changement de chromosome)
		
			line_C=C.strip().split(":")
			chrom_C=line_C[0].replace(">","") #Chromosome chez le chimpanzé
			interval=line_C[1].split("-")	#Intervalle st-end chimpanzé
			
			if "+" in C: 
				 start=int(interval[0])	#Start chimpanzé pour brin +
				 strand="+" #Brin chimpanzé
			else: 
				start =int(interval[1].replace("(",""))	#Start chimpanzé pour le brin -
				strand ="-" 
				 
			#~ print(chrom)
		
		else:	#Lignes de séquence
			
			for n in range(len(H)):#Pour chaque position de nucléotide dans la séquence
				mut=[] #Réinitialise la liste mut à chaque nucléotide 
				if H[n].upper() in nucleotides: 
					mut.append(chrom_C)
					if strand == "+":		
						pos=start+n #Calcule la position du nucléotide 
								
					else:
						pos=start-(n+1) #n+1 car la position end dans les intervalles n'est pas inclue dans la séquence
					
					#.upper() car les nucléotides des éléments répétés sont en minuscule
					if H[n].upper()==C[n].upper(): #Si même nucléotide entre Homme et chimpanzé
						if  C[n].upper()==G[n].upper() and C[n].upper()==O[n].upper(): #Même nucléotide partout
						
							mut.append(pos)
							mut.append(C[n].upper())
							
							out.write("{}\t{}\t{}\n".format(mut[0],mut[1],mut[2])) 
							mut=[] #Réinitialise la liste mut à chaque nucléotide 
						
			start=pos #Mets à jour la position start pour les cas ou on a un retour à la ligne pour une même séquence dans le fichier fasta


def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-H', '--input_H', type=str, help='Path to Human sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/H_seq.fa")					
	parser.add_argument('-C', '--input_C', type=str, help='Path to Chimp sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/C_seq.fa")			
	parser.add_argument('-G', '--input_G', type=str, help='Path to Gorilla sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/G_seq.fa")			
	parser.add_argument('-O', '--input_O', type=str, help='Path to Orangutan sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/O_seq.fa")			
	
	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step4_results/AB_chimp.txt")			

	
	args = parser.parse_args()
	
	print("start searching for ancestral bases")
	same_base(args.input_H, args.input_C, args.input_G, args.input_O,args.output)
	
	
if "__main__" == __name__:
	main()
