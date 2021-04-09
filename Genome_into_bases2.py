
#!/usr/bin/python3

#Importer packages: 
import argparse 



def load_fasta(fa_file):
	"""
	Charge les sequences fasta dans un dictionnaire par chromosome
	"""
	print("start loading fasta file")
	data=open(fa_file,"r")
	genome={} #Initialise le dico qui contiendra les chromosomes en clé
	seq=[] #Initialise la liste qui contiendra les séquences d'un chromosome
	a=0
	
	for l in data: 
		line=l.strip() #Retire les \n en fin de ligne
		if line.startswith(">"): #Lignes >chr:st-end
			if a == 1:
				genome[chrom]=''.join(seq) #Relie la liste des séquences d'un chromosome en une string
				seq=[]
				
			chrom=line.replace(">","")  #Chromosome
			a=1
			print(chrom) 

		else:	#Lignes de séquence
			seq.append(line) #Ajoute la séquence à la liste du chromosome
				
	return genome
			
			
def read_fasta(genome, output):
	"""
	Lis les séquences et écrit les bases dans le fichier de sortie
	"""
	print("start reading fasta file")	
	nucleotides=["A","C","T","G"]
	chromosomes=["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
	base=[] 
	
	for chrom in genome.keys(): 
		if chrom in chromosomes:
		for i in range(len(genome[chrom])):
			if genome[chrom][i].upper() in nucleotides: 
				pos=i #L'index = position du nucléotide
				nuc=genome[chrom][i].upper()
				
				base=[]
				base.append(chrom) #Chromosome
				base.append(pos) #Position de la base
				base.append(nuc) #Type de nucléotide
	
				out.write("{}\t{}\t{}\n".format(base[0],base[1],base[2])) #Ecrit en format Chr \t position \t nucléotide \n



def main(): 
	parser = argparse.ArgumentParser()
	
	#fichier fasta des séquences:
	parser.add_argument('-fa', '--input_fa', type=str, help='Path to fasta file', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step3_results/hg38.fa")					
		
	
	#fichier de sortie: 
	parser.add_argument('-out', '--output', type=str, help='Path to output file', default ="/media/disk1/soukkal/StageM2/Stage_M1/Human_Step5_results/genome_bases_2.txt")			

	
	args = parser.parse_args()
	
	
	Genome=load_fasta(args.input_fa)
	read_fasta(Genome, args.output)

	print("finished writing output")
	
	
if "__main__" == __name__:
	main()
