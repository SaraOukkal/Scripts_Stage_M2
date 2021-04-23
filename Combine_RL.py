
#!/usr/bin/python3

#Importer packages: 
import argparse 

def combine(Left, Right):
	data1=open(Left,"r")
	data2=open(Right,"r")
	
	L=data1.readlines()
	R=data2.readlines()
	dico={}
	
	for a in L: 
		if not a.startswith("d"):
			lineL=a.strip().split("\t")
			distL=int(lineL[0])
			A_L=int(lineL[1])
			C_L=int(lineL[2])
			G_L=int(lineL[3])
			T_L=int(lineL[4])
			
			for i in range(index,len(R)):
				if not R[i].startswith("d"):
					lineR=R[i].strip().split("\t")
					distR=int(lineR[0])
					if distL == distR:
						index=i 
						T_R=int(lineR[1])
						G_R=int(lineR[2])
						C_R=int(lineR[3])
						A_R=int(lineR[4])
						
						A=A_L + A_R
						C=C_L + C_R
						G=G_L + G_R
						T=T_L + T_R
						
						dico[distR]={}
						dico[distR]["A"]=A
						dico[distR]["C"]=C
						dico[distR]["G"]=G
						dico[distR]["T"]=T
						
	return dico
	
def write_out(dico, output):
	
	out=open(output,"w")		
			
	out.write("{}\t{}\t{}\t{}\t{}\n".format("dist","A","C","G","T"))#Ecrit un header au fichier 
	
	for dist in sorted(dico.keys()):
		out.write("{}\t{}\t{}\t{}\t{}\n".format(str(dist),str(dico[dist]["A"]),str(dico[dist]["C"]),str(dico[dist]["G"]),str(dico[dist]["T"])))		
						
			
				
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

	

if "__main__" == __name__:
	main()
