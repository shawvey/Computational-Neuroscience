#for the submission uncomment the submission statements
#so submission.README

from math import *
import numpy as np
from submission import *

def seven_segment(pattern):

    def to_bool(a):
        if a==1:
            return True
        return False
    

    def hor(d):
        if d:
            print(" _ ")
        else:
            print("   ")
    
    def vert(d1,d2,d3):
        word=""

        if d1:
            word="|"
        else:
            word=" "
        
        if d3:
            word+="_"
        else:
            word+=" "
        
        if d2:
            word+="|"
        else:
            word+=" "
        
        print(word)

    pattern_b=list(map(to_bool,pattern))

    hor(pattern_b[0])
    vert(pattern_b[1],pattern_b[2],pattern_b[3])
    vert(pattern_b[4],pattern_b[5],pattern_b[6])

    number=0
    for i in range(0,4):
        if pattern_b[7+i]:
            number+=pow(2,i)
    print(int(number))
    
def weight_matrix(listT):
    length = len(listT)
    weight = np.zeros([length,length])
    for i in range(length):
        for j in range(i,length):
            if i==j:
                weight[i,j] = 0
            else:
                weight[i,j] = listT[i]*listT[j]
                weight[j,i] = weight[i,j]
    return weight


def update(pattern,weight):
    nextPattern = [0,0,0,0,0,0,0,0,0,0,0]
    sums = 0
    for i in range(11):
        sums = 0
        for j in range(11):
            if i!=j:
                sums+= (pattern[j] * weight[i,j])
        if sums > 0:
            nextPattern[i]= 1
        else:
            nextPattern[i]= -1
    return nextPattern

def Cal_energy(pattern, weight):
    sums = 0
    for i in range(11):
        for j in range(11):
            sums += pattern[i] * weight[i,j] * pattern[j]
    energy = -1/2 * sums
    energy = round(energy,3)
    return energy

def hopfield_network(pattern,weight):
    seven_segment(pattern)
    PrePattern=[0 for i in range(0,11)]
    PrePatternCopy=[0 for i in range(0,11)]
    submission.seven_segment(pattern)
    submission.print_number(Cal_energy(pattern,weight))
    #here the network should run printing at each step
    #for the final submission it should also output to submission on each step
    while(pattern != PrePattern):
        PrePattern = PrePatternCopy
        pattern = update(pattern,weight)
        PrePatternCopy = pattern
        seven_segment(pattern)
        submission.seven_segment(pattern)
        ##where energy is the energy of test
        submission.print_number(Cal_energy(pattern,weight))
        

submission=Submission("Xiaoyue Xiao")
submission.header("Xiaoyue Xiao")

six=[1,1,-1,1,1,1,1,-1,1,1,-1]
three=[1,-1,1,1,-1,1,1,1,1,-1,-1]
one=[-1,-1,1,-1,-1,1,-1,1,-1,-1,-1]

seven_segment(three)
seven_segment(six)
seven_segment(one)

##this assumes you have called your weight matrix "weight_matrix"
submission.section("Weight matrix")
weight = (weight_matrix(one) + weight_matrix(three) + weight_matrix(six)) / 3.0
submission.matrix_print("W",weight)


## Test1

submission.section("Test 1")
test=[1,-1,1,1,-1,1,1,-1,-1,-1,-1]
hopfield_network(test,weight)
##for COMSM0027
##this prints a space
submission.qquad()



## Test2

submission.section("Test 2")
test=[1,1,1,1,1,1,1,-1,-1,-1,-1]
hopfield_network(test,weight)
##for COMSM0027
##this prints a space
submission.qquad()



submission.bottomer()



