import matplotlib.pyplot as plt 

def file_(file_name):
	file0=open(file_name,"r")
	x={}
	i=0
	for line in file0:
		L=line.split()
		if len(L)==2:
			time=i*0.1
			qhat=float(L[1])
			if time not in x:
				x[time]=0
			x[time]+=qhat
			i+=1
		else:
			i=0
	x1=[]
	y=[]
	for j in x:
		x1.append(j)
		y.append(x[j]/1e4)
	return x1,y

def MAIN():
	x0,y0=file_("/Users/ritobandatta/OneDrive - Wayne State University/v2-related-/doc_file/qhat_0-5_saved_profile_10_hydro.txt")
	x1,y1=file_("/Users/ritobandatta/OneDrive - Wayne State University/v2-related-/doc_file/qhat_0-5_concurrent_profile_10_hydro.txt")

	plt.plot(x0,y0,label="Saved profiles",linestyle="--",marker="o",color="red")
	plt.plot(x1,y1,label="Generated Hydro",linestyle="--",marker="o",color="blue")
	plt.legend()
	plt.xlabel("time(fm)")
	plt.ylabel("$\\hat{q}/\\sqrt{2}$")
	plt.xlim(0.45,9)
	plt.title("0-5% $\\hat{q}$ comparison; 100GeV Gluon fired using PGun at 1000 different $\\phi$s; 10Hydro Profiles;")
	plt.show()

if __name__=="__main__":
	MAIN()


