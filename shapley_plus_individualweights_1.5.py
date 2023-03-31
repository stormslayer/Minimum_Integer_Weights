#!/usr/bin/env python

# cite: "Pivotality Measures", 2023, scott de marchi, scott.demarchi@gmail.com

import numpy
import math
import sys
import re
import itertools
from pulp import *
import csv 
import sys
import glob
import pandas as p

# GLOBALS
# main parameters for Shapley+
# used in shapley_gauss -- this is the Shapley+ measured used in our empirical work
R = 10000
Enoise = .165 # this is derived from all cases post-WW2 for use in Shapley_gauss 
 
 
# ------------------------------------------------------------------------------------------------------
# FUNCTIONS

#define recursive coalition combinatorics generator
#note this finds ALL coalitions, not just MWC's or min. weight coal.'s
def combo(tcoal, tlim, n, L1, tl):
    # tcoal is coalition list, tlim is # of parties in coal, n is starting party
    # and L1 is list of coalitions used to fill dictionary coal_list
    global all_list, pulp_list
    for j in range(n,T-tlim+1):
        tl.append(j)
        L1.append(tcoal[j])
        if tlim-1>0:
            combo(tcoal, tlim-1, j+1, L1, tl)
        if tlim==1:
            tsum=0
            for z in range(0,len(L1)): tsum=tsum+L1[z]
            coal_list[tuple(L1)]=tsum
            all_list.append(tuple(tl))
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


def scombo(tcoal, tlim, n, L1, tl):
    # tcoal is coalition list, tlim is # of parties in coal, n is starting party, lenT is length of coalition passed
    # and L1 is list of coalitions used to fill dictionary coal_list
    global s_list
    for j in range(n,T-tlim+1):
        tl.append(j)
        L1.append(tcoal[j])
        if tlim-1>0:
            scombo(tcoal, tlim-1, j+1, L1, tl)
        if tlim==1:
            tsum=0
            for z in range(0,len(L1)): tsum=tsum+L1[z]
            s_list.append(tuple(tl))
        L1.pop()
        tl.pop() 
    return 0

# find shapley_plus neighborhood method

def pshapley(c):
    global s_list

    # reset s_list -- change later to make local !!!
    s_list = []
    unity = []
    for i in range(0,T):
        unity.append(i)
    s_list.append(tuple(unity))   
    #L is temp list for coalition calculator combo
    L=[]
    #temp_list is to record variable places for repeat coalitions for pulp
    temp_list=[]
    for i in range(1,T):
        scombo(c,i,0,L, temp_list)
    
    # find new u given s_list
    # generate simple majority u
    su=0
    for i in c:
        su=su+i
 
    if (su/2) == 0: su=su/2+1
    else: su=math.floor(su/2)+1

    # init and set new shapley to 0
    s = []
    for i in range(0, T):
        s.append(0.0)
        
    for i in s_list:
        temp_sum = 0
        for z in i:
            temp_sum = temp_sum + c[z]
        if temp_sum >=su:
            for j in i:
                if (temp_sum - c[j]) < su:
                    s[j] = s[j] + float(math.factorial(len(i)-1)*math.factorial(T-len(i)))/math.factorial(T)
    return s


def find_neighbors(c):
    # find 1 step neighborhood -- call multiple times to find n steps
    
    # p retains all coalitions 1 step away from c
    p = []
    
    # first to +1 step
    for i in range(0,T):
        # perturb up 1 step
        c[i] = c[i] + 1
        for j in range(0,T):
            # now, adjust another party down 1 step to balance weights
            if j != i and c[j]>0:
                c[j] = c[j] - 1
                p.append(list(c))
                c[j] = c[j] + 1
        # return to baseline
        c[i] =  c[i] - 1

    # next do -1 step
        # first to +1 step
    for i in range(0,T):
        # perturb up 1 step
        if c[i] > 1:
            c[i] = c[i] - 1
            for j in range(0,T):
                # now, adjust another party up 1 step to balance weights
                if j != i and c[j]>0:
                    c[j] = c[j] + 1
                    p.append(list(c))
                    c[j] = c[j] - 1
            # return to baseline
            c[i] =  c[i] + 1
        
    return p


def shapley_neighbor(step):
    
    # keep track of number of coalitions new_coal from perturbing by step -- to normalize shapley_neighbor at the end
    count = 0
    splus = []  # define vector of shapley_plus weights
    global all_list, temp_list, L, shapley_step
    # temp list and total list of all neighbors within step distance of original vector of weights
    temp_list = []
    total_list = []
    
    
    # expand the coalitions all_list given different step sizes by calling find_neighbors
    # start with 1 step
    temp_list = find_neighbors(coal)
    total_list = list(temp_list)
    
    # then proceed if step > 1 by iteratively calling find_neighbors()
    for z in range(1, step):
        temp_list = list(total_list)
        for i in temp_list:
            t2 = find_neighbors(list(i))
            for j in t2:
                total_list.append(list(j))
    
        # get rid of duplicates    
        for i in range(0, len(total_list)):
            total_list[i] = tuple(total_list[i])
        total_list = list(set(total_list))
    
    #print "total_list ", total_list
    
    # now iterate through total_list and add to shapley weights
    for i in range(0, len(total_list)):
            # call shapley w/ new vector of weights
            new_coal = list(total_list[i])
            temp_s = pshapley(new_coal)
            #print new_coal, temp_s
            for k in range(0, T):
                shapley_step[k] = shapley_step[k] + temp_s[k]
            count = count + 1

    # last, add in unperturbed vector of weights
    new_coal = list(coal)
    temp_s = pshapley(new_coal)
    for k in range(0, T):
        shapley_step[k] = shapley_step[k] + temp_s[k]
    count = count + 1
    #print new_coal, temp_s
    # now normalize shapley_neighbor by count
    # maybe normalize further so vector sums to 1 !!!
    for i in range(0, T):
        shapley_step[i] = float(shapley_step[i]) / count
    
    return 0


def shapley_gauss_individual():
    # this does Shapley+ which we use for applied stats models -- more efficient than neighborhood model
    # keep track of number of coalitions new_coal from perturbing by step -- to normalize shapley_neighbor at the end
    # find shapley_plus Gaussian method w/ individual noise
    # data from Germany / Finland shows st dev is Enoise*seat_share w/ an r^2 of .87 in-sample

    global all_list, temp_list, L
    count = 0
    splus = []  # define vector of shapley_plus weights
    sumw = 0
    isigma = [] # define vector of indivdiual st. dev.'s for sigma
    
    # first calculate total sum of weights for use in normalizing after perturbing
    for i in range(0, T):
        sumw = sumw + coal[i]
    
    # next, determine column vector of sigmas for each party
    for i in range(0,T): 
        isigma.append(float(Enoise*coal[i]))
    
    # for R times, perturb the weights -- I don't like the use of R...    
    for i in range(0, R):
        new_coal = list(coal)
        count = count + 1
        for j in range(0, T):
            # draw weight for each party -- note that numpy.random.normal uses sd rather than var
            new_coal[j] = int(round(numpy.random.normal(new_coal[j], math.sqrt(isigma[j])),0))
            # boundary check
            if new_coal[j] < 0:
                new_coal[j] = 0

        # now normalize to keep same number of seats as originally using largest remainder method to assign left over seats
        tsum = 0
        p = []
        remainder = []
        for j in range(0, T):
            tsum = tsum + new_coal[j]
        for j in range(0, T):
            p.append(float(new_coal[j]) / tsum)
        for j in range(0, T):
            new_coal[j] = int(round(p[j]*sumw,0))
            remainder.append((p[j]*sumw)-new_coal[j])
        tsum = 0
        for j in range(0, T):
            tsum = tsum + new_coal[j]
        tdelta = 0
        tdelta = sumw - tsum
        # totals of new seats are off of original total
        if tdelta < 0:
            for k in range(0, abs(tdelta)):
                new_coal[remainder.index(max(remainder))] = new_coal[remainder.index(max(remainder))] + 1
                remainder[remainder.index(max(remainder))] = 0.0

        # call shapley w/ new vector of weights
        temp_s = pshapley(new_coal)
        for k in range(0, T):
            shapley_sigma[k] = shapley_sigma[k] + temp_s[k]
        #print new_coal, shapley_sigma
    # last, add in unperturbed vector of weights
    new_coal = list(coal)
    temp_s = pshapley(new_coal)
    for k in range(0, T):
        shapley_sigma[k] = shapley_sigma[k] + temp_s[k]
    count = count + 1
    #print new_coal, shapley_sigma
    # now normalize shapley_neighbor by count
    # maybe normalize further so vector sums to 1 !!!
    for i in range(0, T):
        shapley_sigma[i] = float(shapley_sigma[i]) / count
    return 0


def miw_find():
    # find minimum integer weights for a given vector of raw weights
    # this is a bit of a kludge; there must be a better way
    global pulp_list, tie_list, loss_list, miw

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

    tkeys=MWC.keys()

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
        tindex = int(re.sub('weight_','',v.name))
        miw[tindex]=v.varValue
    
    return 0


def banzhaf_find():
    # find banzhaf values
    # depends on globals all_list, coal

    tpivotal = []
    tbanzhaf = []

    for i in range(0, T):
        tpivotal.append(0)
        tbanzhaf.append(0.0)

    for i in all_list:
        temp_count = 0
        # find sum of coalition values
        #print i
        for j in i:
            temp_count = temp_count + coal[j]
        # now check to see who is pivotal 
        for j in i:
            if temp_count >= u and temp_count - coal[j] < u:
                #print "peep", j
                tpivotal[j]=tpivotal[j]+1

    # create total number of pivots from above
    psum = 0

    for i in tpivotal:
        psum = psum + i

    #create banzhaf index
    for i in range(0, len(tbanzhaf)):
        tbanzhaf[i] = float(tpivotal[i]) / psum

    return tbanzhaf


def shapley_find():
    # find shapley values
    # use all coalitions including unity

    # depends on globals all_list, coal

    tshapley = []

    for i in range(0, T):
        tshapley.append(0.0)

    for i in all_list:
        temp_sum = 0
        # add weights for coalition i
        for z in i:
            temp_sum = temp_sum + coal[z]
        # if i is winning then proceed
        if temp_sum >=u:
            for j in i:
                # check for pivotal party and then apply formula w/ permutations
                if (temp_sum - coal[j]) < u:
                    tshapley[j] = tshapley[j] + float(math.factorial(len(i)-1)*math.factorial(T-len(i)))/math.factorial(T)
    return tshapley


def slack_find(tslack):
    # this finds the value of slack for each coalition
    # there are, however, many options for slack -- using raw seats, salience weighted, using a pivotality measure, etc.
    # the unidimensional policy position here is a mean of the vector of policies from CHES 

    # depends on global temp_data

    for i in tslack.keys():
        # get coalition mean for mwc i
        tcoal_mean = 0.0
        # and get slack for each coalition -- slack starts off at the whole dollar retained before paying out policy
        tslack[i] = 1.0
        for j in i:
            tcoal_mean = tcoal_mean + temp_data[j][7]   # add mean policy position for each party
        tcoal_mean = tcoal_mean / len(i)
        for j in i:
            tslack[i] = tslack[i] - abs(tcoal_mean - temp_data[j][7])

    return tslack


def leverage_find(tleverage):    
    # now calculate leverage
    # depends on global variable slack and that slack only has MWCs

    for i in tleverage.keys():
        for j in i:
            # remove one party from MWC
            tcoal = list(i)
            tcoal.remove(j)

            # now that one party is removed, see if the remaining members of the coal can form another MWC
            # keep best alternative slack if such a MWC exists
            best_alt_slack = 0.0
            for k in tleverage.keys():
                if k != i:
                    if set(tcoal) < set(k): 
                        if slack[k] > best_alt_slack:
                            best_alt_slack = slack[k]
            tleverage[i][j] = slack[i] - best_alt_slack 
    
    print("hi")

    return tleverage


# ---------------------------------------------------------------------------------------------------------
# main code

# reset output files
sys.stdout = open('data_gmd1.csv', 'w') 
print("party, seats, miw_old, sq_cabinet, sq_pm, country, election, mean_policy, base, miw_new, banzhaf, shapley, shapley+, max_slack, max_slack_without")
# reset stdout to screen
sys.stdout = sys.__stdout__

output = open("slack_coal.csv", 'w')
output.write("case;year;mwc;slack_coal;leverage_list")
output.write('\n')
output.close()

# input nations file
case_data_input = []
model_params = p.read_excel('model_params_measures.xlsx')
for line in range(model_params.shape[0]):
    case_data_input.append(model_params.iloc[line,0])

for z in case_data_input:
    f1 = open("c:\\users\\karuthers\\dropbox\\bargaining_abm_ai\\data_input\\all_cases\\"+str(z), 'r')
    csvobj = csv.reader(f1)

    # put data into list of lists
    dataset = []
    header = 1              # use to delete first row of input file which is a header
    for row in csvobj:
        if header == 0: dataset.append(row)
        header = 0
        
    # each row has party,seats,miw,position,salience,alpha,sq_cabinet,sq_pm,country,election, govt_nmbr
    # get raw and miw weights from each country file
    temp_data = []

    # set T global for the number of parties
    T = len(dataset)

    for i in range(0, T):
        # get party, seats, miw, sq_cabinet, sq_pm, country, election year, CHES policies
        temp_data.append(['.','.','.','.','.','.','.','.'])
        temp_data[i][0]=dataset[i][0]   # party name                                                      
        temp_data[i][1]=dataset[i][1]   # seats                           
        temp_data[i][2]=dataset[i][2]   # miw                                                  
        temp_data[i][3]=dataset[i][6]   # sq_cab
        temp_data[i][4]=dataset[i][7]   # sq_pm                                                  
        temp_data[i][5]=dataset[i][8]   # country
        temp_data[i][6]=dataset[i][9]   # election year

        # check for '.' (missing policies / saliences) in input data
        temp_pol = dataset[i][3].split(",")
        temp_pol2 = []
        for j in range (0, len(temp_pol)):
            # assuming policy and salience vectors are symmetric wrt to missingness, but this doesn't change raw data in dataset[] for either
            try: 
                temp_pol2.append(float(temp_pol[j]))
            except:
                pass
        mean_policy = sum(map(float, temp_pol2))/len(temp_pol2)
        temp_data[i][7] = mean_policy   # mean policy from CHES
    
    f1.close()
        
    # now generate raw seat shares for the rest of the code to use -- this is in place of user input
    coal = []
    for i in range(0, T):
        coal.append(int(temp_data[i][1]))

    # version 1.0 coalition calculator
    # it's not pretty, but it works -- cleaner version 1.1 in library 


    # not used b/c of file input above
    #coal=input("Please enter the coalition in the form [x,y,z...]: ")
    #step = input("Please enter the maximum step size for Shapley neighborhood: ")
    #sigma = input("Please enter variance for Shapley Gaussian: ")
    #print "Working on this coalition: ",coal, " with step size ", step, " with variance ", sigma, " with R ", R
    print("Working on this coalition: ",coal)
    print("/n")

    miw=list(coal)

    # T is total size of parties -- done above for file input
    #T=len(coal)

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
    # also has MWCs with indices of parties; mwc dictionary (below) has MWCs by raw seats
    pulp_list=[]

    # establish a list of ties for pulp constraints
    tie_list=[]

    #establish losing list
    loss_list=[]

    # s_list for use w/ shapley neighborhood and gaussian
    s_list = []

    # list of all coals, not just MWC's.  Also does NOT remove duplicates like coal_list
    all_list = []
    unity = []
    for i in range(0,T):
        unity.append(i)
    all_list.append(tuple(unity))    

    #print "# Parties / Total Votes / Total votes needed for coalition: ",T," / ",E," / ",u

    #L is temp list for coalition calculator combo
    L=[]
    #temp_list is to record variable places for repeat coalitions for pulp
    temp_list=[]

    for i in range(1,T):
        combo(coal,i,0,L, temp_list)

    #calculate largest possible coalition for use in rankings

    max_size=0
    for i in range (0, len(pulp_list)):
        temp1=len(pulp_list[i])
        if temp1>max_size: max_size=temp1

    # now that we have all coalitions, check to see those that are MWC's
    # create new dictionary MWC
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

    # number of distinct minimum winning coalitions
    num_coals=len(MWC)

    # create global slack dictionary
    slack = {}
    for i in pulp_list:
        slack[i] = -1.0

    # get slack for each MWC
    slack = slack_find(slack)

    # create leverage dictionary
    leverage = {}
    temp_leverage = []
    for i in range(0, T):
        temp_leverage.append(0.0)
    # for each viable non-0 slack MWC, initialize player leverages -- all non-members of any given coalition remain at 0
    for i in slack.keys():
        leverage[i] = temp_leverage.copy()
    # get leverage now that we have slack
    leverage = leverage_find(leverage)


    # -----------------------------------------------------------------------------------------------

    # START shapely_plus work
    # u is threshold for winning
    # all_list is a list of all coalitions
    # T is the number of parties
    # coal is list of coalition weights
    # create list by party of number of pivotal outcomes
    pivotal = []
    banzhaf = []
    shapley = []
    shapley_step = []
    shapley_sigma = []

    for i in range(0, T):
        pivotal.append(0)
        banzhaf.append(0.0)
        shapley.append(0.0)
        shapley_step.append(0.0)
        shapley_sigma.append(0.0)

    # find MIWs
    miw_find()
    print("base: ", coal)
    print("miw: ", miw)

    # find banzhaf
    banzhaf = banzhaf_find()
    print("Banzhaf: ", banzhaf)

    # find shapley
    shapley = shapley_find()
    print("Shapley: ", shapley)

    # this does neighborhoods version of shapley
    #for x in range(1,step):
    #    shapley_neighbor(x)
    #    print "Shapley neighbor ",x, " steps: ", shapley_step
        
    #shapley_neighbor(10)
    #print shapley_step

    # call gaussian for INDIVIDUAL sigma for all parties
    shapley_gauss_individual()
    print("Shapley+: ", shapley_sigma)

    # compute max_slack for each party and max_slack_without -- for printing to datafile below
    max_slack = []
    max_slack_without = []
    for j in range(0, T):
        max_slack.append(0.0)
        max_slack_without.append(0.0)
    
    for j in slack.keys():
        for i in range(0, T):
            if (i in j) and slack[j] > max_slack[i]: max_slack[i] = slack[j]
            if not(i in j) and slack[j] > max_slack_without[i]: max_slack_without[i] = slack[j]


    # print in correct format -- all data
    # header first and change stdout to a file instead of the screen
    sys.stdout = open('data_gmd1.csv', 'a')
    for i in range(0, T):
        temp_str = []
        for j in temp_data[i]:
            temp_str.append(str(j))
        print (", ".join(temp_str), ",", coal[i], ",", miw[i], ",", banzhaf[i], ",", shapley[i], ",", shapley_sigma[i], ",", max_slack[i], ",", max_slack_without[i])
    # reset stdout to screen
    sys.stdout = sys.__stdout__

    # now append slack at the coalition level to a separate file
    output = open("slack_coal.csv", 'a')
    for i in slack.keys():
        output.write(str(temp_data[0][5]))  # country -- using index 0 b/c same for all parties in a given case
        output.write(";")
        output.write(str(temp_data[0][6]))  # year
        output.write(";")
        output.write(str(i))
        output.write(";")
        output.write(str(slack[i]))
        output.write(";")
        for j in range(0, T):
            output.write(str(leverage[i][j]))
            output.write(";")
        output.write('\n')
    output.close

