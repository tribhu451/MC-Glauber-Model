


import matplotlib.pyplot as plt
import sys
#plt.style.use('seaborn-whitegrid')


 


# (1) 
file1 = open(sys.argv[1], 'r')
Lines = file1.readlines()

count = 0
for line in Lines:
	count += 1
	if count > 1 :
		values = [float(s) for s in line.split()]
		if values[0] == 40 :
			ref_mult = values[2]

	else :
		continue
                    
#ASE1 = [LE1, UE1]


file1 = open(sys.argv[1], 'r')
Lines = file1.readlines()

count = 0
X1, Y1, LE1, UE1 , UB1, LB1 = [],[],[],[],[],[]
for line in Lines:
	count += 1
	if count > 1 :
		values = [float(s) for s in line.split()]
		X1.append( (values[0]+values[1] ) / 2)                       
		Y1.append(values[2]/ref_mult) 
		#LE1.append(values[2]) 
		#UE1.append(values[2]) 
		#UB1.append(values[1]+values[2]) 
		#LB1.append(values[1]-values[2]) 
 
	else :
		continue
                    
#ASE1 = [LE1, UE1]


# .............................................................

# (2) 
file1 = open(sys.argv[2], 'r')
Lines = file1.readlines()

count = 0
for line in Lines:
	count += 1
	if count > 1 :
		values = [float(s) for s in line.split()]
		if values[0] == 45 :
			ref_mult = values[5]

	else :
		continue
                    
#ASE1 = [LE1, UE1]


file1 = open(sys.argv[2], 'r')
Lines = file1.readlines()

count = 0
X2, Y2, LE2, UE2 , UB2, LB2 = [],[],[],[],[],[]
for line in Lines:
	count += 1
	if count > 1 and count < 11  :
		values = [float(s) for s in line.split()]
		X2.append( values[0] )                       
		Y2.append(values[5]/ref_mult) 
		#LE2.append(values[2]) 
		#UE2.append(values[2]) 
		#UB2.append(values[1]+values[2]) 
		#LB2.append(values[1]-values[2]) 
 
	else :
		continue
                    
#ASE2 = [LE2, UE2]



# .............................................................

# (2) 
file1 = open(sys.argv[3], 'r')
Lines = file1.readlines()

count = 0
for line in Lines:
	count += 1
	if count > 1 :
		values = [float(s) for s in line.split()]
		if values[0] == 45 :
			ref_mult = values[5]

	else :
		continue
                    
#ASE1 = [LE1, UE1]


file1 = open(sys.argv[3], 'r')
Lines = file1.readlines()

count = 0
X3, Y3, LE3, UE3 , UB3, LB3 = [],[],[],[],[],[]
for line in Lines:
	count += 1
	if count > 1 and count < 11  :
		values = [float(s) for s in line.split()]
		X3.append( values[0] )                       
		Y3.append(values[5]/ref_mult) 
		#LE2.append(values[2]) 
		#UE2.append(values[2]) 
		#UB2.append(values[1]+values[2]) 
		#LB2.append(values[1]-values[2]) 
 
	else :
		continue
                    
#ASE2 = [LE2, UE2]



X4,Y4 = [],[]
X4.append(-10)
Y4.append(1)
X4.append(100)
Y4.append(1)






#         PLOTS
fig = plt.figure(figsize=(10, 8))
ax = fig.gca()

ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(axis='both', which='minor', labelsize=15)


#plt.fill_between(X1,LB1, UB1 ,facecolor='r',alpha=0.5, color = 'black')
#plt.fill_between(X2,LB2, UB2 ,facecolor='r',alpha=0.5, color = 'blue')
#plt.fill_between(X3,LB3, UB3 ,facecolor='r',alpha=0.5, color = 'green')
#plt.fill_between(X4,LB4, UB4 ,facecolor='r',alpha=0.3, color = 'gray')
#plt.fill_between(X5,LB5, UB5 ,facecolor='r',alpha=0.3, color = 'gray')

plt.scatter(X1, Y1,marker = 'o', s =100, color = 'black',linewidth=2.0,facecolors='none',label = 'ALICE')
plt.scatter(X2, Y2,marker = '*', s =60, color = 'red',linewidth=1.0,facecolors='none', label = 'mc glauber model')
plt.scatter(X3, Y3,marker = 's', s =60, color = 'green',linewidth=2.0,facecolors='none', label = 'optical glauber model')
#plt.scatter(X3, Y3,marker = '*', s =30, color = 'green',linewidth=1.0,facecolors='none', label = 'double FO (160MeV + 150MeV )')
#plt.scatter(X4, Y4,marker = 's', s =30, color = 'navy',linewidth=1.0,facecolors='none', label = 'double FO (160MeV + 140MeV )')
#plt.scatter(X5, Y5,marker = 'D', s =30, color = 'magenta',linewidth=1.0,facecolors='none', label = 'double FO (155MeV + 150MeV )')

#plt.errorbar(X1, Y1, yerr = ASE1,fmt = 'none',ms = 5, color = 'black',elinewidth=1.5)
#plt.errorbar(X2, Y2, yerr = ASE1,fmt = 'none',ms = 5, color = 'blue',elinewidth=1.0)
#plt.errorbar(X3, Y3, yerr = ASE1,fmt = 'none',ms = 5, color = 'green',elinewidth=0.5)
plt.plot(X4, Y4, ls = '--', color = 'grey')

plt.xlabel(r'$centrality (\%)$', fontsize = 20)
plt.ylabel(r'$\frac{dN_{ch}}{d\eta}$ ratio',fontsize = 20)
plt.yscale('linear')
plt.legend(fontsize = 28)

plt.xlim(0.,80.)
#plt.ylim(0.1,8)


plt.text(25, 3 , 'Pb+Pb 2760 GeV \n xhard = 0.14', fontsize=40, color = 'blue')


plt.savefig('mult_ratio.pdf', format = 'pdf')
#plt.savefig('eta_vs_ch_x.eps', format = 'eps')
#plt.savefig('eta_vs_ch_X.png', format = 'png')





