#!/usr/bin/env python

#!/usr/bin/env python
import numpy
import math
import sys
import string
from pulp import *

# version 1.0 coalition calculator
# it's not pretty, but it works -- cleaner version 1.1 in library 

#coal=list(input("Please enter the coalition in the form [x,y,z...]: "))
coal = [141,113,49,25,21]
print("Working on this coalition: ",coal)
print()

miw=list(coal)

# T is total size of parties
T=len(coal)

# generate simple majority u
u=0
for i in coal:
    u=u+i
E=u   
if (u/2) == 0: u=u/2+1
else: u=math.floor(u/2)+1

# generate variable names for parties
varnames=[]
for i in range(0,T):
    tempname=i
    varnames.append(tempname)
    tempname=""

# establish dictionary w/ list of all coalitions and weights
coal_list={}

# establish a list that will include weights for each T variables for pulp
pulp_list=[]

# establish a list of ties for pulp constraints
tie_list=[]

#establish losing list
loss_list=[]

#define recursive coalition combinatorics generator
#note this finds ALL coalitions, not just MWC's or min. weight coal.'s
def combo(tcoal, tlim, n, L1, tl):
    # tcoal is coalition list, tlim is # of parties in coal, n is starting party
    # and L1 is list of coalitions used to fill dictionary coal_list
    for j in range(n,T-tlim+1):
        tl.append(j)
        L1.append(tcoal[j])
        if tlim-1>0:
            combo(tcoal, tlim-1, j+1, L1, tl)
        if tlim==1:
            tsum=0
            for z in range(0,len(L1)): tsum=tsum+L1[z]
            coal_list[tuple(L1)]=tsum
            if E%2==0 and tsum==(u-1):
                tie_list.append(tuple(tl))
            elif tsum>=u:
                pulp_list.append(tuple(tl))
                for k in L1:
                    if tsum - k >= u:
                        pulp_list.pop()
                        break
            else:
                loss_list.append(tuple(tl))
        L1.pop()
        tl.pop() 
    return 0


#print "# Parties / Total Votes / Total votes needed for coalition: ",T," / ",E," / ",u

#L is temp list for coalition calculator combo
L=[]
#temp_list is to record variable places for repeat coalitions for pulp
temp_list=[]

for i in range(1,T):
    combo(coal,i,0,L, temp_list)
        

print("Pulp List")
print(pulp_list)
print()
print("Tie List")
print(tie_list)
print()
print("Loss List")
print(loss_list)
print()

#calculate largest possible coalition for use in rankings

max_size=0
for i in range (0, len(pulp_list)):
    temp1=len(pulp_list[i])
    if temp1>max_size: max_size=temp1

#now that we have all coalitions, check to see those that are MWC's
#create new dictionary MWC
MWC = {}
MWC=coal_list.copy()

# delete non-winning coalitions
for i in list(MWC.keys()):
    if MWC[(i)]<u: del MWC[(i)]


# delete too large coalitions
for i in list(MWC.keys()):
    temp_coal=i
    for j in temp_coal:
        if MWC[(i)] - j >= u:
            del MWC[(i)]
            break
    

print("MWC: ")
print(MWC)
print() 

#this is a kludge; there must be a better way
# assign ordinal ranks to the different parties for use later by pulp
rank=[]
for i in range(0, T):
    rank.append(0)
    
t_mwc=len(pulp_list)
k=2
t_dummy=0
while k<max_size+1:
    for i in range(0, len(pulp_list)):
        for j in pulp_list[i]:
            if len(pulp_list[i])==k:
                rank[j]=rank[j]+((max_size-k)*(max_size-k)+1)*t_mwc
                t_dummy=t_dummy+1
    ranked=1
    t_mwc=t_mwc-t_dummy
    t_dummy=0
    for l in range(0,T-1):
        if rank[l]==rank[l+1] and coal[l]!=coal[l+1]:
            k=k+1
            ranked=0
            break
    if ranked==1: break   

print("rank: ", rank)
print()



tkeys=MWC.keys()
#organize MWC's by party
#for i in range(0,T):
#    print "Player ", i+1, "'s coalitions:"
#    for j in range(0, len(tkeys)):
#        if coal[i] in tkeys[j]: print "Coalition: ",tkeys[j]," / Size: ",MWC[tkeys[j]]
#    print


num_coals=len(MWC)
#print "Number of distinct coalitions: ", num_coals

# pulp code

#create filters for weights

party_filter=[]
tie_filter=[]
loss_filter=[]

row = []
for i in range(0, len(pulp_list)):
    for j in range(0,len(varnames)):
        row.append(-1)
    party_filter.append(row)
    row=[]
    
for i in range (0, len(pulp_list)):
    for j in pulp_list[i]:
        party_filter[i][j]=1
    
row = []
for i in range(0, len(tie_list)):
    for j in range(0,len(varnames)):
        row.append(-1)
    tie_filter.append(row)
    row=[]
    
for i in range (0, len(tie_list)):
    for j in tie_list[i]:
        tie_filter[i][j]=1

row = []
for i in range(0, len(loss_list)):
    for j in range(0,len(varnames)):
        row.append(-1)
    loss_filter.append(row)
    row=[]
    
for i in range (0, len(loss_list)):
    for j in loss_list[i]:
        loss_filter[i][j]=1

# Create the 'prob' variable to contain the problem data
prob = LpProblem("Minimum Integer Weights", LpMinimize)

coal_vars = LpVariable.dicts("weight",varnames, 0, None, LpInteger)

# The objective function is added to 'prob' first
prob += lpSum([coal_vars[i] for i in varnames]), "Minimize Weights"

# The constraints are added to 'prob'
for i in range(0, len(pulp_list)):
    prob += lpSum([party_filter[i][j] * coal_vars[j] for j in varnames]) >= 1, "coalition constraint "+str(i)

for i in range(0, len(tie_list)):
    prob += lpSum([tie_filter[i][j] * coal_vars[j] for j in varnames]) == 0, "tie constraint "+str(i)

for i in range(0, len(loss_list)):
    prob += lpSum([loss_filter[i][j] * coal_vars[j] for j in varnames]) <= -1, "loss constraint "+str(i)


# do rank constraints for pulp
for k in range(0, T-1):
    if rank[k]>rank[k+1]:
        prob += coal_vars[k] - coal_vars[k+1] >= 1, "rank constraint "+str(k)
    elif rank[k]==rank[k+1]:
        prob += coal_vars[k] - coal_vars[k+1] == 0, "rank constraint "+str(k)    

if rank[T-1]>0: prob += coal_vars[T-1] >=1, "final rank constraint"

for k in varnames:
    if rank[k]>0:
        prob += coal_vars[k]>=1, "final constraint "+str(k)

#The constraints are entered
#prob += w1 + w2 - w3 - w4 >=1, "coalition 1"
#prob += w1 - w2 + w3 - w4 >=1, "coalition 2"
#prob += w1 - w2 - w3 + w4 >=1, "coalition 3"
#prob += -w1 + w2 + w3 + w4 >=1, "coalition 4"
#prob += w1 - w2 >=0, "coalition 5"
#prob += w2 - w3 >=0, "coalition 6"
#prob += w3 - w4 >=0, "coalition 7"
#prob += w4 >=0, "coalition 8"


# The problem data is written to an .lp file
prob.writeLP("bargain1.lp")

# The problem is solved using PuLP's choice of Solver
prob.solve()

# The status of the solution is printed to the screen
# print "Status:", LpStatus[prob.status]

# Each of the variables is printed with it's resolved optimum value
for v in prob.variables():
#    print v.name, "=", v.varValue
    tindex = int(v.name.lstrip("weight_"))
    miw[tindex]=v.varValue

    
# The optimised objective function value is printed to the screen
# print "Total = ", value(prob.objective)

print("base: ", coal)
print("miw: ", miw)
