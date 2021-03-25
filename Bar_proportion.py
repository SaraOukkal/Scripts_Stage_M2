#!/usr/bin/python3

#Importer packages: 
import argparse 


def bar_genome(input_bar):
	data=open(input_bar,"r")
	tot_size=0
	
	for l in data: 
		line=l.strip().split('\t')
		chrom_br=line[0] #Chromosome
		st_br=int(line[1]) #Start barrière
		end_br=int(line[2]) #End barrière	
		bar_size=end_br-st_br
		tot_size+=bar_size
		
	print("Barriers total size:", tot_size)
	return tot_size
	

def bar_intervals(input_bar_inter):
	data=open(input_bar_inter,"r")
	tot_ali_size=0
	
	for l in data: 
		line=l.strip().split('\t')
		chrom_br=line[0] #Chromosome
		st_br=int(line[1]) #Start barrière
		end_br=int(line[2]) #End barrière	
		bar_size=end_br-st_br
		tot_ali_size+=bar_size
		
	print("Barriers in intervals total size:", tot_ali_size)
	return tot_ali_size
	
def proportion(bar_genome, bar_intervals):
	prop=bar_intervals/bar_genome 
	print("Proportion of barriers included in aligned intervals:", prop)
	percentage=(bar_intervals/bar_genome)*100
	print("Percentage of barriers included in aligned intervals:", percentage)


			
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	parser.add_argument('-bar', '--input_bar', type=str, help="")					
	parser.add_argument('-bar_inter', '--input_bar_inter', type=str, help="")	

	args = parser.parse_args()
	
	bar=bar_genome(args.input_bar)
	bar_inter=bar_intervals(args.input_bar_inter)
	proportion(bar,bar_inter)

if "__main__" == __name__:
	main()
