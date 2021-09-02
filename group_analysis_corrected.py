import os,sys
import numpy as np
import nibabel as nib
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

directory=sys.argv[1]
outfile=sys.argv[2]

def partial_group_analysis(directory,images,start_vol, end_vol):
	arr=[]
	
	for img in images:
		#print(img.get_filename())
		arr.append(img.slicer[:,:,:,start_vol:end_vol].dataobj)
		
	#print(len(arr))
	#print(arr[0].shape)
	#print("append complete")
	
	t,p=stats.ttest_1samp(arr,0,axis=0)

	del arr

	print("Tests complete")
	z=stats.norm.ppf(p/2)
	z=-1*z*np.sign(t)
	z[np.isnan(z)]=0
	#print(z.shape)
	p[np.isnan(p)]=1
	t[np.isnan(t)]=0
	fin_q=[]
	fin_p=-1*np.log10(p)*np.sign(t)
	#print("fin_p shape ", fin_p.shape)
	shape = fin_p.shape
	sign = np.sign(t)
	p=p.reshape(-1,shape[3])
	#print("p shape ", p.shape)
	p = np.transpose(p)
	for pval in p:
		#print(pval.shape)
		r,q=fdrcorrection(pval,alpha=0.05,method="indep",is_sorted=False)
		fin_q.append(q)
	fin_q = np.transpose(fin_q)
	fin_q=np.reshape(fin_q,(shape[0],shape[1],shape[2],shape[3]))
	fin_q = -1*np.log10(fin_q)*sign
	fin_q[np.isnan(fin_q)]=0
	# corr_z=z*np.transpose(rr)
	#print(fin_q[fin_q>1], fin_q[fin_q>1].shape)
	#print(fin_q.shape)

	return fin_p,fin_q,z,t



def group_analysis(directory):
	# 1000 * 100 * 100 * 100 * 12 = 96 gb
	images = []
	min_4th_dim = 10000
	#count = 0 
	for file in os.listdir(directory):
		#count+=1
		#if(count>30):
		#	break
		img = nib.load(directory+"/"+file)
		img.set_filename(file)
		images.append(img) 
		if min_4th_dim>img.shape[-1]:
			min_4th_dim = img.shape[-1]
		print(img.get_filename())
		del img
	print(len(images))
	filename=directory+"/"+images[0].get_filename()
	print(filename)
	print("images loaded")
	
	ps=[]
	qs=[]
	zs=[]
	ts=[]
	
	start_vol = 0
	end_vol = 3 
	while(start_vol<min_4th_dim):
		pvals1,qvals1,zvals1,tvals1=partial_group_analysis(directory,images,start_vol,min(end_vol,min_4th_dim))
		ps.append(pvals1)
		qs.append(qvals1)
		zs.append(zvals1)
		ts.append(tvals1)
		print(start_vol,end_vol,pvals1.shape)
		start_vol += 3
		end_vol += 3	

	pvals = np.array(ps[0])
	qvals = np.array(qs[0])
	zvals = np.array(zs[0])
	tvals = np.array(ts[0])
	
	for i in range(1, len(ps)) :
		pvals = np.concatenate((pvals,ps[i]),axis=3)
		qvals = np.concatenate((qvals,qs[i]),axis=3)
		zvals = np.concatenate((zvals,zs[i]),axis=3)
		tvals = np.concatenate((tvals,ts[i]),axis=3)
		print("number : ",i,ps[i].shape,pvals.shape)
		
	print(pvals.shape,qvals.shape,zvals.shape,tvals.shape)
	del ps,qs,zs,ts
	return pvals,qvals,zvals,tvals,filename


	
def save_files(p4d,q4d,z4d,t4d,outfile,affine,header,shape,dtype):
	print("started saving files")
	p4d=p4d.astype(dtype)
	p4dimg=nib.Nifti1Image(p4d, affine,header)
	del p4d
	nib.save(p4dimg,outfile+'/group_p_values.nii.gz')
	del p4dimg
	print("saved group_p_values")
	q4d=q4d.astype(dtype)
	q4dimg=nib.Nifti1Image(q4d, affine,header)
	del q4d
	nib.save(q4dimg,outfile+'/group_q_values.nii.gz')
	del q4dimg
	print("saved group_q_values")
	z4d=z4d.astype(dtype)
	z4dimg=nib.Nifti1Image(z4d, affine,header)
	del z4d
	nib.save(z4dimg,outfile+'/group_z_values.nii.gz')
	del z4dimg
	print("saved group_z_values")
	t4d=t4d.astype(dtype)
	t4dimg=nib.Nifti1Image(t4d, affine,header)
	del t4d
	nib.save(t4dimg,outfile+'/group_t_values.nii.gz')
	del t4dimg
	print("saved group_t_values")
	return

pvals,qvals,zvals,tvals,filename=group_analysis(directory)

print("group analysis done, now loading file data")

fmri_4D=nib.load(filename)
shape=np.array(fmri_4D.get_fdata()).shape
dtype=fmri_4D.get_data_dtype()
header=fmri_4D.header
affine=fmri_4D.affine
del fmri_4D

print("saving files")

if not os.path.exists(outfile):
	os.makedirs(outfile)

save_files(pvals,qvals,zvals,tvals,outfile,affine,header,shape,dtype)

print("done")


