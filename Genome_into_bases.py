
#!/usr/bin/python3

#Importer packages: 
import argparse 



def same_base(fa_file, output_file):
	"""
	Rapporte les positions des bases alignées
	"""
	
	out=open(output_file,"w")
	data=open(fa_file,"r")
	
	done_chrom=[]
	nucleotides=["A","C","T","G"]
	
	start=0
	
	for l in data: 
		
		if l.startswith(">"): #Lignes >chr:st-end
			chrom=l.strip().replace(">","")  #Chromosome
			base=[]
			start=0
			
			if chrom not in done_chrom:
				done_chrom.append(chrom)
				print(chrom) #Sert à savoir ou on en est dans le fichier (au changement de chromosome)
			
		else:	#Lignes de séquence
			for n in range(len(l)):#Pour chaque position de nucléotide dans la séquence
				base=[] #Réinitialise la liste base à chaque nucléotide 
				if l[n].upper() in nucleotides: 
					pos=start+n #Calcule la position du nucléotide 
					base.append(chrom) #Chromosome chez l'espèce d'interet 	
					base.append(pos)
					base.append(l[n].upper())

					out.write("{}\t{}\t{}\n".format(base[0],base[1],base[2])) 
					base=[] #Réinitialise la liste à chaque nucléotide 
				elif l[n].upper() =='N' :
					pos=start+n
				else:
					break
		
		start=pos #Mets à jour la position start pour les cas ou on a un retour à la ligne pour une même séquence dans le fichier fasta


def main(): 
	parser = argparse.ArgumentParser()
	
	#fichier fasta des séquences:
	parser.add_argument('-fa', '--input_fa', type=str, help='Path to fasta file', default ="/home/soukkal/Bureau/Projet/Step3_results/H_seq.fa")					
		
	
	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step4_results/Ali_bases.txt")			

	
	args = parser.parse_args()
	
	print("start processing")
	same_base(args.input_fa,args.output)
	print("finished")
	
	
if "__main__" == __name__:
	main()
