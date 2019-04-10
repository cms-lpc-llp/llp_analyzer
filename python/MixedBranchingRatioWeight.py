
def getBRWeight(x, y, n_t, n_chipm):
    # x = BR(g->bbchi0)
    # y = BR(g->ttchi0)
    # 1-x-y = BR(g->tbchi+-)
    assert x+y <= 1
    assert n_t <= 4
    assert n_chipm <= 2
    
    weight = 0
    weight += (x*x)*            (n_t==0 and n_chipm==0) # g->bbchi0  g->bbchi0
    weight += (y*y)*            (n_t==4 and n_chipm==0) # g->ttchi0  g->ttchi0 
    weight += ((1-x-y)*(1-x-y))*(n_t==2 and n_chipm==2) # g->tbchi+- g->tbchi+-
    weight += (2*x*y)*          (n_t==2 and n_chipm==0) # g->ttchi0  g->bbchi0
    weight += (2*x*(1-x-y))*    (n_t==1 and n_chipm==1) # g->tbchi+- g->bbchi0
    weight += (2*y*(1-x-y))*    (n_t==3 and n_chipm==1) # g->ttchi0  g->tbchi+-
    return weight

def getBRWeightString(x,y):
    # x = BR(g->bbchi0)
    # y = BR(g->ttchi0)
    
    assert x+y <= 1
    
    weightStrings = []
    
    if x*x>0:
        weightStrings.append('(%f*%f)*(ntFromGluino==0 && nCharginoFromGluino==0)'%(x,x))            # g->bbchi0  g->bbchi0
    if y*y>0:
        weightStrings.append('(%f*%f)*(ntFromGluino==4 && nCharginoFromGluino==0)'%(y,y))            # g->ttchi0  g->ttchi0
    if (1-x-y)*(1-x-y)>0:
        weightStrings.append('(%f*%f)*(ntFromGluino==2 && nCharginoFromGluino==2)'%(1-x-y,1-x-y))    # g->tbchi+- g->tbchi+-
    if x*y>0:
        weightStrings.append('(2.*%f*%f)*(ntFromGluino==2 && nCharginoFromGluino==0)'%(x,y))         # g->ttchi0  g->bbchi0
    if x*(1-x-y)>0:
        weightStrings.append('(2.*%f*%f)*(ntFromGluino==1 && nCharginoFromGluino==1)'%(x,1-x-y))     # g->tbchi+- g->bbchi0
    if y*(1-x-y)>0:
        weightStrings.append('(2.*%f*%f)*(ntFromGluino==3 && nCharginoFromGluino==1)'%(y,1-x-y))     # g->ttchi0  g->tbchi+-

    weightString = '+'.join(weightStrings)
    weightString = '(' + weightString + ')'
    
    return weightString

if __name__ == '__main__':

    for x in [0,0.25,0.50,0.75,1]:
        for y in [0,0.25,0.50,0.75,1]:
            if x+y>1: continue
            print '(%s/%s)'%(getBRWeightString(x,y),getBRWeightString(0.25,0.25))
