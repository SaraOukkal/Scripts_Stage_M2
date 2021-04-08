
#!/usr/bin/python3

#Importer packages: 
import argparse 



def load_fasta(fa_file):
	"""
	Lis le fichier fasta et stock dans un dictionnaire les séquences des chromosomes
	"""
	print("start loading fasta file")
	data=open(fa_file,"r")
	genome={}
	
	for l in data: 
		if l.startswith(">"): #Lignes >chr:st-end
			chrom=l.strip().replace(">","")  #Chromosome
			print(chrom) 

		else:	#Lignes de séquence
			print(l)
			
			if chrom in genome.keys():
				genome[chrom]+=l.strip('\n')
				print(genome[chrom])
				break
				
			else:
				genome[chrom]=l.strip('\n')
				
	return genome
			
			
def read_fasta(genome, output):
	print("start reading fasta file")	
	nucleotides=["A","C","T","G"]
	chromosomes={}
	base=[]
	
	for chrom in genome.keys(): 
		for i in range(len(genome[chrom])):
			if genome[chrom][i].upper() in nucleotides: 
				pos=i #L'index = position du nucléotide
				nuc=genome[chrom][i].upper()
				
				base=[]
				base.append(chrom) 
				base.append(pos)
				base.append(nuc)
	
				out.write("{}\t{}\t{}\n".format(base[0],base[1],base[2])) 



def main(): 
	parser = argparse.ArgumentParser()
	
	#fichier fasta des séquences:
	parser.add_argument('-fa', '--input_fa', type=str, help='Path to fasta file', default ="/home/soukkal/Bureau/Projet/Step3_results/H_seq.fa")					
		
	
	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/home/soukkal/Bureau/Projet/Step4_results/Ali_bases.txt")			

	
	args = parser.parse_args()
	
	
	Genome=load_fasta(args.input_fa)
	read_fasta(Genome, args.output)

	print("finished writing output")
	
	
if "__main__" == __name__:
	main()
