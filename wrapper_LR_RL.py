#!/usr/bin/python

import os,sys
import math
basedir="/mnt/scratch1/rahul/HCP_DATA"

outputdir=sys.argv[1]
output_z_dir=sys.argv[2]
task = sys.argv[3]
ev = sys.argv[4]
stuff = sys.argv[5]
mask="/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask.nii.gz"

include_list = os.listdir(basedir)
files4d_LR,evfiles_LR,files4d_RL,evfiles_RL,subjects=[],[],[],[],[]

f=open("HCP_subj_new2.txt","r")
for patient in f:
	item = patient[0:6]
	path_LR=basedir + "/"+item + "/MNINonLinear/Results/"+task+"_LR"
	path_RL=basedir + "/"+item + "/MNINonLinear/Results/"+task+"_RL"
	if os.path.exists(path_LR) and os.path.exists(path_RL):
		print(path_LR)
		print(path_RL)	
		files4d_LR.append(path_LR+ "/"+task+"_LR.nii.gz")
		files4d_RL.append(path_RL+ "/"+task+"_RL.nii.gz")
		evfiles_LR.append(path_LR +"/EVs/"+ev+".txt")
		evfiles_RL.append(path_RL +"/EVs/"+ev+".txt")
		subjects.append(item)

nc=20
k=math.floor(len(subjects)/nc)
k=1

if not os.path.exists(outputdir):
	os.makedirs(outputdir)

if not os.path.exists(output_z_dir):
	os.makedirs(output_z_dir)

print(len(files4d_LR))
print(len(files4d_RL))
print(len(evfiles_LR))
print(len(evfiles_LR))
print(len(subjects))


if not os.path.exists(stuff):
	os.makedirs(stuff)

evfiles,files4d=[],[]

for l,r,ld,rd,sb in zip(evfiles_LR,evfiles_RL,files4d_LR,files4d_RL,subjects):
	if not os.path.exists(stuff+"/"+sb):
		os.makedirs(stuff+"/"+sb)
	evfile=open(stuff+ "/" + sb + "/" + task + "_" + ev + "_ev.txt","w")
	fmrifile=open(stuff + "/" + sb + "/" + task  + "_" + ev + "_files.txt","w")
	evfiles.append(stuff+ "/" + sb + "/" + task  + "_" + ev + "_ev.txt")
	files4d.append(stuff + "/" + sb + "/" + task + "_" + ev  + "_files.txt")
	evfile.write(l)
	evfile.write('\n')
	evfile.write(r)
	evfile.write('\n')
	fmrifile.write(ld)
	fmrifile.write('\n')
	fmrifile.write(rd)
	fmrifile.write('\n')
	evfile.close()
	fmrifile.close()


nc=20
k=math.floor(len(subjects)/nc)
#k = 1

if not os.path.exists(outputdir):
	os.makedirs(outputdir)

if not os.path.exists(output_z_dir):
	os.makedirs(output_z_dir)


for i in range(k):	
	print(i)
	print(i*nc+19)
	os.system('python tsa_it_31_05.py -f -i {0} -o {1} -t {2} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {4} -o {5} -t {6} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {7} -o {8} -t {9} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {10} -o {11} -t {12} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {13} -o {14} -t {15} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {16} -o {17} -t {18} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {19} -o {20} -t {21} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {22} -o {23} -t {24} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {25} -o {26} -t {27} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {28} -o {29} -t {30} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {31} -o {32} -t {33} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {34} -o {35} -t {36} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {37} -o {38} -t {39} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {40} -o {41} -t {42} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {43} -o {44} -t {45} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {46} -o {47} -t {48} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {49} -o {50} -t {51} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {52} -o {53} -t {54} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {55} -o {56} -t {57} -m {3} -w 0 25 &'
			  'python tsa_it_31_05.py -f -i {58} -o {59} -t {60} -m {3} -w 0 25 & wait'.format(files4d[i*nc], outputdir+"/"+subjects[i*nc],evfiles[i*nc],mask,
			  	                                                               files4d[i*nc+1], outputdir+"/"+subjects[i*nc+1],evfiles[i*nc+1],
			  	                                                               files4d[i*nc+2], outputdir+"/"+subjects[i*nc+2],evfiles[i*nc+2],
			  	                                                               files4d[i*nc+3], outputdir+"/"+subjects[i*nc+3],evfiles[i*nc+3],
			  	                                                               files4d[i*nc+4], outputdir+"/"+subjects[i*nc+4],evfiles[i*nc+4],
			  	                                                               files4d[i*nc+5], outputdir+"/"+subjects[i*nc+5],evfiles[i*nc+5],
			  	                                                               files4d[i*nc+6], outputdir+"/"+subjects[i*nc+6],evfiles[i*nc+6],
                                                                         files4d[i*nc+7], outputdir+"/"+subjects[i*nc+7],evfiles[i*nc+7],
			  	                                                               files4d[i*nc+8], outputdir+"/"+subjects[i*nc+8],evfiles[i*nc+8],
			  	                                                               files4d[i*nc+9], outputdir+"/"+subjects[i*nc+9],evfiles[i*nc+9],
			  	                                                               files4d[i*nc+10], outputdir+"/"+subjects[i*nc+10],evfiles[i*nc+10],
			  	                                                               files4d[i*nc+11], outputdir+"/"+subjects[i*nc+11],evfiles[i*nc+11],
			  	                                                               files4d[i*nc+12], outputdir+"/"+subjects[i*nc+12],evfiles[i*nc+12],
			  	                                                               files4d[i*nc+13], outputdir+"/"+subjects[i*nc+13],evfiles[i*nc+13],
			  	                                                               files4d[i*nc+14], outputdir+"/"+subjects[i*nc+14],evfiles[i*nc+14],
			  	                                                               files4d[i*nc+15], outputdir+"/"+subjects[i*nc+15],evfiles[i*nc+15],
			  	                                                               files4d[i*nc+16], outputdir+"/"+subjects[i*nc+16],evfiles[i*nc+16],
			  	                                                               files4d[i*nc+17], outputdir+"/"+subjects[i*nc+17],evfiles[i*nc+17],
			  	                                                               files4d[i*nc+18], outputdir+"/"+subjects[i*nc+18],evfiles[i*nc+18],
			  	                                                               files4d[i*nc+19], outputdir+"/"+subjects[i*nc+19],evfiles[i*nc+19]))


for j in range(k*nc,len(subjects)):
	print(j)
	os.system('python tsa_it_31_05.py -f -i {0} -o {1} -t {2} -m {3} -w 0 25 & wait'.format(files4d[j], outputdir+"/"+subjects[j],evfiles[j],mask))
	
for patient in os.listdir(outputdir):
	print(patient)
	zfile=outputdir+"/"+patient + "/z_statistic.nii.gz"
	os.system('cp {0} {1} & wait'.format(zfile,output_z_dir+"/"+patient+"_z_statistic.nii.gz"))
	






