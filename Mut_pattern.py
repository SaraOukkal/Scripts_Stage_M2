#!/usr/bin/python3

#Importer packages: 
import argparse 



def mut_pattern(fa_file1, fa_file2, fa_file3, fa_file4, output_file):
	
	fini=False
	
	out=open(output_file,"w")
	
	H_file=open(fa_file1,"r")
	C_file=open(fa_file2,"r")
	G_file=open(fa_file3,"r")
	O_file=open(fa_file4,"r")
	
	done_chrom=[]
	
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
				mut=[] #Réinitialise la liste mut à chaque nouvelle mutation 
				mut.append(chrom_C)
				if strand == "+":		
					pos=start+n #Calcule la position du nucléotide 
							
				else:
					pos=start-(n+1) #n+1 car la position end dans les intervalles n'est pas inclue dans la séquence
				
				#.upper() car les nucléotides des éléments répétés sont en minuscule
				if H[n].upper()!=C[n].upper(): #Si nucléotide différent entre Homme et chimpanzé
					if  C[n].upper() != G[n].upper() or  C[n].upper() != O[n].upper() : #S'il y a peut être eu mutation dans la lignée du chimpanzé
	
						mut.append(pos)
						mut.append(C[n])
						
						if H[n].upper()==G[n].upper() and H[n].upper()==O[n].upper(): #Si H = G = O alors on détermine l'état ancestral
							mut.append(H[n])
						else:
							mut.append("N") #Sinon on ne peut pas conclure sur l'état ancestral, on écrit N
						
						mut.append(H[n])
						mut.append(G[n])
						mut.append(O[n])
						
			
						out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(mut[0],mut[1],mut[2],mut[3],mut[4],mut[5],mut[6])) 
			
			start=pos #Mets à jour la position start pour les cas ou on a un retour à la ligne pour une même séquence dans le fichier fasta


def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers fasta des séquences qui s'alignent chez les 4 espèces:
	parser.add_argument('-H', '--input_H', type=str, help='Path to Human sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/H_seq.fa")					
	parser.add_argument('-C', '--input_C', type=str, help='Path to Chimp sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/C_seq.fa")			
	parser.add_argument('-G', '--input_G', type=str, help='Path to Gorilla sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/G_seq.fa")			
	parser.add_argument('-O', '--input_O', type=str, help='Path to Orangutan sequences', default ="/home/soukkal/Bureau/Projet/Step3_results/O_seq.fa")			
	
	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step3_results/mut_pattern.txt")			

	
	args = parser.parse_args()
	
	print("start searching for mutations")
	mut_pattern(args.input_H, args.input_C, args.input_G, args.input_O,args.output)
	


	
if "__main__" == __name__:
	main()
