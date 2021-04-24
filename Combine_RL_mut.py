
#!/usr/bin/python3

#Importer packages: 
import argparse 

def combine(Left, Right):
	data1=open(Left,"r")
	data2=open(Right,"r")
	
	L=data1.readlines()
	R=data2.readlines()
	dico={}
	index=0
	
	for a in L: 
		if not a.startswith("d"):
			lineL=a.strip().split("\t")
			distL=int(lineL[0])
			AT_L=int(lineL[1])
			AC_L=int(lineL[2])
			AG_L=int(lineL[3])
			CT_L=int(lineL[4])
			CA_L=int(lineL[5])
			CG_L=int(lineL[6])
			GT_L=int(lineL[7])
			GA_L=int(lineL[8])
			GC_L=int(lineL[9])
			TA_L=int(lineL[10])
			TC_L=int(lineL[11])
			TG_L=int(lineL[12])
			
			for i in range(index,len(R)):
				if not R[i].startswith("d"):
					lineR=R[i].strip().split("\t")
					distR=int(lineR[0])
					if distL == distR:
						index=i 
						AT_R=int(lineL[1])
						AC_R=int(lineL[2])
						AG_R=int(lineL[3])
						CT_R=int(lineL[4])
						CA_R=int(lineL[5])
						CG_R=int(lineL[6])
						GT_R=int(lineL[7])
						GA_R=int(lineL[8])
						GC_R=int(lineL[9])
						TA_R=int(lineL[10])
						TC_R=int(lineL[11])
						TG_R=int(lineL[12])
						
						AT=AT_L + TA_R
						AC=AC_L + TG_R
						AG=AG_L + TC_R
						CT=CT_L + GA_R
						CA=CA_L + GT_R
						CG=CG_L + GC_R
						GT=GT_L + CA_R
						GA=GA_L + CT_R
						GC=GC_L + CG_R
						TA=TA_L + AT_R
						TC=TC_L + AG_R
						TG=TG_L + AC_R
						
						dico[distR]={}
					
						dico[distR]["AT"]=AT
						dico[distR]["AC"]=AC
						dico[distR]["AG"]=AG
						dico[distR]["CT"]=CT
						dico[distR]["CA"]=CA
						dico[distR]["CG"]=CG
						dico[distR]["GT"]=GT
						dico[distR]["GA"]=GA
						dico[distR]["GC"]=GC
						dico[distR]["TA"]=TA
						dico[distR]["TC"]=TC
						dico[distR]["TG"]=TG
						break 

						
	return dico
	
def write_out(dico, output):
	
	out=open(output,"w")		
			
	out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("dist","A>T","A>C","A>G","C>T","C>A","C>G","G>T","G>A","G>C","T>A","T>C","T>G"))#Ecrit un header au fichier 
	
	for dist in sorted(dico.keys()):
		out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(str(dist),str(dico[dist]["AT"]),str(dico[dist]["AC"]),str(dico[dist]["AG"]),str(dico[dist]["CT"]),str(dico[dist]["CA"]),str(dico[dist]["CG"]),str(dico[dist]["GT"]),str(dico[dist]["GA"]),str(dico[dist]["GC"]),str(dico[dist]["TA"]),str(dico[dist]["TC"]),str(dico[dist]["TG"])))		
						
			
				
def main(): 
	parser = argparse.ArgumentParser()
	
	#fichiers input:
	##fichiers compte d'AB R et L:				
	parser.add_argument('-L', '--Left', type=str, help='Path to bases count around left borders')			
	parser.add_argument('-R', '--Right', type=str, help='Path to bases count around left borders')		
	#fichier de sortie: 
	parser.add_argument('-RL', '--RightLeft', type=str, help='Path to combined output file')	
	
	args = parser.parse_args()
	
	#Lancer fonctions: 
	both_borders=combine(args.Left, args.Right)
	write_out(both_borders, args.RightLeft)
	

if "__main__" == __name__:
	main()
