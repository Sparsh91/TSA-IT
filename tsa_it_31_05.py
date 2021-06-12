import numpy as np
np.seterr(divide='ignore', invalid='ignore')
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import nilearn
from nilearn.masking import apply_mask
from nilearn.input_data import NiftiMasker
import nibabel as nib
import sys
import math
import time
import os

help_msg="usage : tsa-it [-h] [-f] -i Inputfile -o Outputfile [-m mask] [-w window] [-n] [-b] [-t trials_file] [-x v1,v2,v3..] [--version] \n Required Arguments: \n -i <INPUT> (4D input files in (NIfTI/Analyze/DICOM/) format)\n -o <OUTPUT> (Path of the output directory in which the output files (p-values, corrected p-values, z-values, corrected z-values, t-values, mean, standard deviation of percentage signal change from baseline) will be saved; the directory is created if does not already exist) \n  -t or -x option (see below for their functionality) \n Optional arguments:\n-f, --file  If this option is given, then the -i <INPUT> is a file containing the names of input files and -t <trials.txt> is a file containing the list of EV files corresponding to the input files (-x option cannot be used in this case). These files contain multiple runs for the same subject and the trials in these need to be combined to produce the final inter-trial TSA output\n -m <MASK>, --mask <MASK> 3D mask file for masking input data. If <MASK> is a number then input files are thresholded by that number to produce the mask\n -w start end, --window start end The time window (in seconds) from start to end, relative to the trials (start can also be negative, which means that the window starts before the trial; end can also be negative or positive) for which to compute the ITS (defaults to 25 seconds after the end of the task condition, which is the maximum duration for which hemodynamic response function can influence the BOLD signal)\n -n, --normal       Consider a null hypothesis of normal distribution model on subjects to calculate ITS probabilities (default)\n -b, --binomial   Consider a null hypothesis binomial model on subjects to calculate ITS probabilities \n -t <trials.txt>  EV file in single-column or three-column format. \n -x v1,v2,v3,  ...    A list of volume numbers indicating the beginning of trials (the volume indexing starts from zero). Cannot be used with -t. \n --version Show program's version number and exit"

def time_series(filename,mask):
	nifti_masker=None
	if mask==None:
		# nifti_masker=NiftiMasker(mask_strategy='epi')
		# nifti_masker.fit(filename)
		# mask=nifti_masker.mask_img_
		# masked_data=apply_mask(filename,mask)
		img = nib.load(filename)
		masked_data=((img.get_fdata()).moveaxis(-1,0)).reshape(img.shape[3],-1)
		print("Warning: No mask is provided, performing analysis on the entire brain \n")


	elif type(mask)==str:
		nifti_masker=NiftiMasker(mask_img=mask)
		nifti_masker.fit(filename)
		mask=nifti_masker.mask_img_
		masked_data=apply_mask(filename,mask)

	elif type(mask)==float:
		fmri_4D  = nib.load(filename)
		data=fmri_4D.get_fdata()
		brain=data[:,:,:,0]
		brain=np.sort(brain[brain>mask])
		x,y=brain[0],brain[-1]
		del brain,data,fmri_4D
		threshold=mask/y
		nifti_masker=NiftiMasker(mask_strategy='epi',mask_args=dict(lower_cutoff=threshold))
		nifti_masker.fit(filename)
		mask=nifti_masker.mask_img_
		masked_data=apply_mask(filename,mask)

	return masked_data,nifti_masker

def trials_single(trial_file,n):
	ev1=open(trial_file,"r")
	initial_trials=[]
	line=ev1.readline()
	min_duration=[]
	bool_var=True
	while line:
		trial=line.split()
		if len(trial) !=1 and len(trial)!=3:	
			print("Error: ",trial_file," is neither in 1-column format nor 3-column format")	
			print(help_msg)	
			sys.exit()
		initial_trials.append(trial[0])
		line=ev1.readline()
		if len(trial)==3:
			min_duration.append(float(trial[1]))
			if min(min_duration) < float(trial[1]):
				bool_var=False
	if bool_var==False:
		print("Warning: the trial durations are different")
	if len(min_duration)==0:
		min_duration.append(0)
	#print(min(min_duration))
	return initial_trials,min(min_duration),bool_var

def trials_multiple(trial_file):
	files=open(trial_file,"r")
	trials=[]
	file=files.readline()
	min_durations=[]
	n=1
	while file:
		initial_trials,min_duration,bool_var=trials_single(file.split()[0],n)
		n=n+1
		min_durations.append(min_duration)
		trials.append(initial_trials)
		file=files.readline()
		if bool_var==False:
			fin_bool=False
	

	return trials,min(min_durations)

def correct_start(trial_starts,trial_ends):
	diff=trial_ends-trial_starts
	diff=diff-min(diff)
	return trial_starts+diff

def signal_change_util(ts,trials,start,end,tr):
	ts=ts.transpose()
	fin_arr=[]
	trial_starts=np.floor((trials+start)/float(tr)).astype(np.int32)
	# trial_ends=np.floor((trials+end)/float(tr)).astype(np.int32)
	trial_ends=trial_starts + math.floor((end-start)/float(tr))
	print(trials)
	
	print(trial_ends)
	trial_starts=correct_start(trial_starts,trial_ends)
	print(trial_starts)

	x1,x2,end=trials+start,trials+end,trial_ends-1
	# print(x1[x1>0])
	# print(x2[x2<tr*ts.shape[1]])
	# print(trial_starts[trial_starts>0])
	# print(end[end<tr*ts.shape[1]])
	for voxel in ts:
		arr=[]
		for  (trial_start,trial_end) in zip(trial_starts,trial_ends):
			if trial_end < ts.shape[1] and trial_start > 0:
				initial=voxel[(trial_start)-1]
				window=np.array(voxel[(trial_start):(trial_end)])
				percent_change=(window-initial)/initial*100
				arr.append(percent_change)	

		if len(arr)>0:
			arr=np.array(arr).astype(np.float32).transpose()
			fin_arr.append(arr)	

	return fin_arr

def signal_change(ts,trials,start,end,tr):
	changes=[]
	for t_s,trial in zip(ts,trials):
		change=signal_change_util(t_s,np.array(trial).astype(np.float32),start,end,tr)
		print(np.array(change).shape)
		changes.append(change)

	changes=np.concatenate(changes,axis=2)
	return changes

def fdr(logp):
	sign_mask=np.sign(logp)
	logp=np.abs(logp)
	sorted_logp=np.sort(logp)[::-1]
	sorted_logp_indices=np.argsort(-1*logp)
	del logp
	q_vals=np.zeros((sorted_logp.shape[0]))
	length=sorted_logp.shape[0]
	for i in range(sorted_logp.shape[0]):
		q_vals[i]=sorted_logp[i]+ math.log10(i+1) - math.log10(length)

	del sorted_logp
	q_min=-math.inf
	for i in range(length-1,0):
		if q_vals[i]< q_min:
			q_vals[i]=q_min
		else:
			q_min=q_vals[i]

	for i in range(len(q_vals)):
		if q_vals[i] < -1*math.log10(0.05) :
			q_vals[i] = 0		
	q_star=np.zeros((length))
	for val,index in zip(q_vals,sorted_logp_indices):
		q_star[index]=val

	q_star=q_star*sign_mask

	return q_star


def hypothesis_testing_normal(signal):
	t,p=stats.ttest_1samp(signal,0,axis=2)
	z=stats.norm.ppf(p/2)
	z=-1*z*np.sign(t)
	z[np.isnan(z)]=0
	t[np.isnan(t)]=0
	p[np.isnan(p)]=1
	means=np.mean(signal,axis=2)
	std_dev=np.std(signal,axis=2)
	fin_p=-1*np.log10(p)*np.sign(t)
	#shape = fin_p.shape
	#print("aa",shape)
	fin_q=[]
	#p=p.reshape(-1,shape[3])
	p=p.transpose()
	print(p[p>0], p[p>0].shape)
	print(p[p<0], p[p<0].shape)

	# rr=[]
	for pval in p:
		r,q=fdrcorrection(pval,alpha=0.05,method="indep",is_sorted=False)
		fin_q.append(q)	
	# 	rr.append(r)
	#for logp in p:
	#	q=fdr(logp)
	#	fin_q.append(q)
	fin_q = np.transpose(fin_q)
	fin_q = -1*np.log10(fin_q)*np.sign(t)
	#fin_q=np.reshape(fin_q,(shape[0],shape[1],shape[2],shape[3]))
	fin_q[np.isnan(fin_q)]=0
	# corr_z=z*np.transpose(rr)
	print(fin_q[fin_q>1], fin_q[fin_q>1].shape)
	return fin_p,t,z,fin_q,means,std_dev

def hypothesis_testing_binomial(signal):
	fin_p=[]
	for voxel in signal:
		p=[]
		for TR in voxel:
			tr=TR[TR>0]
			pval=stats.binom_test(len(tr),len(TR),p=0.5,alternative='two-sided')
			p.append(pval)
		fin_p.append(p)
		
	t=stats.t.ppf(np.array(fin_p)/2,len(fin_p)-2)
	z=stats.norm.ppf(np.array(fin_p)/2)
	z[np.isnan(z)]=0
	t[np.isnan(t)]=0
	fin_p=np.array(fin_p)
	fin_p[np.isnan(fin_p)]=1
	fin_p=fin_p.transpose()
	fin_q=[]
	rr=[]
	for pval in fin_p:
		r,q=fdrcorrection(pval,alpha=0.05,method="indep",is_sorted=False)
		fin_q.append(q)
		rr.append(r)
	fin_q=np.transpose(fin_q)	
	fin_q=-1*np.log10(fin_q)*np.sign(t)
	fin_q[np.isnan(fin_q)]=1
	fin_p=fin_p.transpose()
	fin_p=-1*np.log10(fin_p)*np.sign(t)
	means=np.mean(signal,axis=2)
	std_dev=np.std(signal,axis=2)
	corr_z=z*np.transpose(rr)
	return fin_p,t,z,fin_q,means,corr_z,std_dev


def save_files_1(p,q,z,t,mean,std_dev,pos,outfile,ext):
	t=pos.inverse_transform(t.transpose())
	p=pos.inverse_transform(p.transpose())
	q=pos.inverse_transform(q.transpose())
	z=pos.inverse_transform(z.transpose())
	mean=pos.inverse_transform(mean.transpose())
	std_dev=pos.inverse_transform(std_dev.transpose())

	t.to_filename(outfile+'/t_statistic'+ext)
	del t
	p.to_filename(outfile+'/p_statistic'+ext)
	del p
	q.to_filename(outfile+'/q_statistic'+ext)
	del q
	z.to_filename(outfile+'/z_statistic'+ext)
	del z
	mean.to_filename(outfile+'/mean_signal_change'+ext)
	del mean
	std_dev.to_filename(outfile+'/standard_deviation'+ext)

	return

def save_files_2(p,q,z,t,mean,std_dev,outfile,ext,shape,header,affine):
	t=t.reshape(shape[0],shape[1],shape[2],shape[3])
	p=p.reshape(shape[0],shape[1],shape[2],shape[3])
	q=q.reshape(shape[0],shape[1],shape[2],shape[3])
	z=z.reshape(shape[0],shape[1],shape[2],shape[3])
	mean=mean.reshape(shape[0],shape[1],shape[2],shape[3])
	std_dev=std_dev.reshape(shape[0],shape[1],shape[2],shape[3])

	t_img=nib.Nifti1Image(t, affine,header)
	del t
	p_img=nib.Nifti1Image(p, affine,header)
	del p
	z_img=nib.Nifti1Image(z, affine,header)
	del z
	q_img=nib.Nifti1Image(q, affine,header)
	del q
	mean_img=nib.Nifti1Image(mean, affine,header)
	del mean
	std_img=nib.Nifti1Image(std_dev, affine,header)
	del std_dev

	nib.save(t_img,outfile+'/t_statistic'+ext)
	nib.save(z_img,outfile+'/z_statistic'+ext)
	nib.save(q_img,outfile+'/q_statistic'+ext)
	nib.save(p_img,outfile+'/p_statistic'+ext)
	nib.save(mean_img,outfile+'/mean_signal_change'+ext)
	nib.save(std_img,outfile+'/standard_deviation'+ext)

	return


if __name__ == '__main__':
	
	time0=time.time()
	args=sys.argv
	del args[0]
	i=0
	filename=None
	outfile=None
	start,end=None,None
	mask=None
	distribution="normal"
	file=""
	trial_file=None
	volumes=[]
	version="1.0"
	a,b,c,d,e,f=False,False,False,False,False,False
	while(i<len(args)):
		if args[i] in ("--input","-i"):
			a=True
			i=i+1
			filename=args[i]
			if i+1<len(args) and args[i+1][0]!="-":	
				print("Error: single file expected \n")	
				print(help_msg)	
				sys.exit()
		if args[i] in ("--output","-o"):
			b=True
			i=i+1
			outfile=args[i]
		if args[i] in ("--window","-w"):
			c=True
			if i+3<len(args):
				if args[i+1][0]!="-" and args[i+2][0]!="-" and args[i+3][0]!="-":
					print("Error: only two arguments required to specify the window size with -w \n")
					print(help_msg)
					sys.exit()
			try:
				i=i+1
				start=float(args[i])
				i=i+1
				end=float(args[i])
			except:
				print("Error: provide start-time and end-time to specify the window size with -w \n")
				print(help_msg)
				sys.exit()
		if args[i] in ("-m","--mask"):
			d=True
			i=i+1
			mask=args[i]
		if args[i] in ("--binomial","-b"):
			distribution="binomial"
		if args[i] in ("--normal","-n"):
			distribution="normal"
		if args[i] in ("--help","-h"):
			 print(help_msg)
			 sys.exit()
		if args[i] in ("--version"):
			print(version)
			sys.exit()
		if args[i] in ("--file","-f"):
			file="file"
		if args[i] in ("-t","--trials"):
			e=True
			i=i+1
			trial_file=args[i]
		if args[i]=="-x":
			f=True
			i=i+1
			volumes=np.array(args[i].split(sep=",")).astype(np.int32)

		i=i+1

	if a==False:	
		print("Error: functional file is a mandatory input argument \n")	
		print(help_msg)	
		sys.exit()	

	if b==False:	
		print("Error: Output directory path is a mandatory input argument \n")	
		print(help_msg)	
		sys.exit()	

	if e==False and f==False:	
		print("Error: one of -t and -x must be specified \n")	
		print(help_msg)	
		sys.exit()

	if os.path.isfile(filename)==False:
		print(filename," Error: ",filename," does not exist \n")
		print(help_msg)
		sys.exit()
	if trial_file!=None and len(volumes)!=0:
		print("Error: -t and -x options cannot be used together \n")
		sys.exit()

	if file=="file" and len(volumes)!=0:
		print("Error: -x option cannot be used along with -f option \n")
		sys.exit()

	if len(volumes)==0 and os.path.isfile(trial_file)==False:	
		print(" Error: ",trial_file," does not exist \n")	
		print(help_msg)	
		sys.exit()

	if mask!=None:
		try:
			mask=float(mask)
		except:
			mask=str(mask)

	
	# if file=="file" and type(mask)!=str:	
	# 	print(" Error: Please provide a mask if using the -f option \n")	
	# 	sys.exit()

	if(type(mask)==str):
		if os.path.isfile(mask)==False:	
			print("The mask file does not exist")	
			#print(help_msg)	
			sys.exit()
		else:
			try:
				x=nib.load(mask).get_fdata().shape
				if len(x)!=3:
					sys.exit()

			except:
				print("Error : The mask has to be a 3D file \n")
				print(help_msg)
				sys.exit()
		

	if file!="file":	
		try:	
			x=nib.load(filename).get_fdata().shape
			if len(x)!=4:	
				sys.exit()
		except:
			print("Error: ",filename,"  is not a functional file \n")	
			print(help_msg)	
			sys.exit()	
			
	if file=="file":
		flag=0
		try:
			f=open(filename,"r")
			x = f.readline()
			while x:
				flag=0
				if not(os.path.isfile(x.split()[0])):
					flag=1
					print("Error: ",x," file does not exist \n")
					sys.exit()
				else:
					y=nib.load(x.split()[0]).get_fdata().shape
					#print(y)
					flag=1
					if len(y)!=4:
						print("Error: ",x," is not a functional file \n")	
						print(help_msg)	
						sys.exit()
		
				x = f.readline()
		except:
			if flag==0:
				print("Error: at least one of the input files is not a functional file \n")
			sys.exit()


	if file=="file" and e==True:
		f=open(trial_file,"r")
		x=f.readline()
		while x:
			if os.path.isfile(x.split()[0])==False:
				print("Error: ",x.split()[0]," does not exist \n")
				sys.exit()
			x=f.readline()
		


	if start!=None and end!=None and start>=end:
		print("Error: The start time must be smaller than the end time \n")
		sys.exit()
	
	if not os.path.exists(outfile):
		os.makedirs(outfile)


	time1=time.time()
	print("initial",time1-time0,filename)

	inputfile=""
	t_series=[]
	trials=[]
	pos=""
	min_dur=0
	# try:
	if file=="":
		inputfile=filename
		ts,pos=time_series(filename,mask)
		t_series.append(ts)
		if trial_file!=None:
			trial,min_dur,_=trials_single(trial_file,"")
			#print(min_dur)
			trials.append(trial)
		elif len(volumes)!=0:
			trials.append(volumes)
		trials=np.array(trials).astype(np.float32)
	# except:
	# 	print(help_msg)
	# 	sys.exit()

	# try:
	if file=="file":
		files=open(filename,"r")
		file=files.readline()
		affine=nib.load(file.split()[0]).affine
		while file:
			affn=nib.load(file.split()[0]).affine
			if affn.all()!=affine.all():
				print("Error: Affine of images don't match \n")
				sys.exit()
			file=files.readline()
			
		if type(mask)!=str:
			mask=0
			print("Warning: Threshold option for mask cannot be used along with -f option. Using the entire brain data for testing.")
		files=open(filename,"r")
		file=files.readline()
		n=0
		while file:
			inputfile=file.split()[0]
			ts,pos=time_series(file.split()[0],mask)
			t_series.append(ts)
			n=n+1
			file=files.readline()
		trials,min_dur=trials_multiple(trial_file)
		print(n,len(trials))
		if n!=len(trials):
			print("Error : The number of 4D files is not equal to the number of EV files \n")
			sys.exit()
	# except:
	# 	print(help_msg)
	# 	sys.exit()


	# try:
	fmri_4D=nib.load(inputfile)
	# print(fmri_4D.header)
	header=fmri_4D.header
	affine=fmri_4D.affine
	shape_data=fmri_4D.get_fdata().shape
	tr=float(header.get_zooms()[-1])
	del fmri_4D
	# except:
	# 	print(help_msg)
	# 	sys.exit()
	
	if start==None and end==None:
		start=0
		end=25+min_dur

	# try:
	if len(volumes)==0:
		fin=[]
		for trial in trials:
			trial=np.array(trial).astype(np.float32)
			rem=trial%tr
			for i in range(len(rem)):
				if rem[i]==0:
					trial[i]=trial[i]-tr
			fin.append(trial.tolist())
		trials=fin

	elif len(volumes)!=0:
		trials=(trials+1)*tr

	# except:
	# 	print(help_msg)
	# 	sys.exit()

	# try:
	time2=time.time()
	print("ts",time2-time1)

	print(start,end)
	signal=signal_change(t_series,trials,start,end,tr)
	

	time3=time.time()
	del t_series
	print("signal",time3-time2,filename)
	fin_p,fin_t,fin_z,fin_q,fin_mean,corr_z=[],[],[],[],[],[]
	if distribution=="normal":
		fin_p,fin_t,fin_z,fin_q,fin_mean,std_dev=hypothesis_testing_normal(signal)
	elif distribution=="binomial":
		fin_p,fin_t,fin_z,fin_q,fin_mean,corr_z,std_dev=hypothesis_testing_binomial(signal)	
	
	

	del signal
	time4=time.time()
	print("hypothesis",time4-time3)

	ext=None
	if ".nii" in inputfile:
		ext=".nii"
	if ".nii.gz" in inputfile:
		ext=".nii.gz"
	if ".img" in inputfile:
		ext=".img"
	if ".mnc" in inputfile:
		ext=".mnc"
	if ".dcm" in inputfile:
		ext=".dcm"

	ext=".nii.gz"

	if d==True:
		save_files_1(fin_p,fin_q,fin_z,fin_t,fin_mean,std_dev,pos,outfile,ext)
	else:
		shape_save = (shape_data[0],shape_data[1],shape_data[2],fin_p.shape[1])
		save_files_2(fin_p,fin_q,fin_z,fin_t,fin_mean,std_dev,outfile,ext,shape_save,header,affine)

	time5=time.time()
	print("files",time5-time4)
	print("total",time5-time0,filename)

	print("tr",tr)


	# except:
	# 	print(help_msg)
	# 	sys.exit()