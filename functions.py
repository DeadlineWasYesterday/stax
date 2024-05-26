import pandas as pd
import numpy as np

def export_block(mt,i,i2):
    t1,t2 = mt.loc[i:,['gi','gi+1','Sequence','Sq+1']],mt.loc[i2:,['gi','gi+1','Sequence','Sq+1']] 
    t1 = t1.iloc[:t2.shape[0],:] #truncate t1
    t1,t2 = t1.reset_index(drop = 1), t2.reset_index(drop = 1)
    b = []
    for t1g1,t1g1p1,t1s1,t1s1p1, t2g1,t2g1p1,t2s1,t2s1p1 in pd.concat([t1,t2], axis = 1).values:
        if (t1g1 == t2g1) and (t1g1p1 == t2g1p1) and (t1s1 == t1s1p1) and (t2s1 == t2s1p1):
            b.append(t1g1)
        else:
            b.append(t1g1)
            break
    return b

def eb_rev(mt,mu,i,i2):
    t1,t2 = mt.loc[i:,['gi','gi+1','Sequence','Sq+1']],mu.loc[:i2,['gi','gi-1','Sequence','Sq-1']][::-1]
    m = min(t1.shape[0],t2.shape[0])
    t1,t2 = t1.iloc[:m,:],t2.iloc[:m,:] #truncate
    t1,t2 = t1.reset_index(drop = 1), t2.reset_index(drop = 1)
    b = []
    for t1g1,t1g1p1,t1s1,t1s1p1, t2g1,t2g1m1,t2s1,t2s1m1 in pd.concat([t1,t2], axis = 1).values:
        if (t1g1 == t2g1) and (t1g1p1 == t2g1m1) and (t1s1 == t1s1p1) and (t2s1 == t2s1m1):
            b.append(t1g1)
        else:
            b.append(t1g1)
            break
    return b

def prep(d1):
    d1 = d1[(d1['Status'] == 'Single') | (d1['Status'] == 'Duplicated')]
    gi = d1.sort_values(['Sequence', 'Gene Start']).drop_duplicates('Gene')['Gene'].reset_index(drop=1).reset_index()
    d1 = d1.sort_values(['Sequence', 'Gene Start'])
    d1 = pd.merge(d1,gi, how = 'left')
    gii = pd.merge(d1,gi, how = 'left')['index']
    d1['gi'] = gii
    d1["Sq+1"] = d1['Sequence'].shift(-1)
    d1["gi+1"] = d1['gi'].shift(-1)
    d1["Sq-1"] = d1['Sequence'].shift(1)
    d1["gi-1"] = d1['gi'].shift(1)
    d1['index'] = d1.index
    return d1

def itr(d1):

    mt = d1[['index', 'Sequence', 'gi', 'Sq+1', 'gi+1']]
    mu = d1[['index', 'Sequence', 'gi', 'Sq-1', 'gi-1']]

    hs = [False] * mt.shape[0]
    t = []
    for i,sq,gi,sq_p1,gi_p1 in mt.values:
        if hs[i]:
            continue
        for i2,sq2,gi2,sq_p12,gi_p12 in mt.loc[i+1:,:].values: #the second hash makes i+1 redundant
            if hs[i2]:
                continue
            if (gi == gi2) and (gi_p1 == gi_p12) and (sq == sq_p1) and (sq2 == sq_p12): #match pair
                #export block
                b = export_block(mt,i,i2)
                #hs[i:i+len(b)] = [True]*len(b) #this would prevent protrusions.
                hs[i2:i2+len(b)] = [True]*len(b)
                hs[i] = [True] #block is defined by first gene. this allows protrusions, i.e., a longer match block on iter 2 after a shorter one.
                #hs[i2] = [True] #this isn't necessary, because the match block can't protrude.
                t.append(('+',i,i2,b,len(b)))

    #reverse
    hs = [False] * mt.shape[0]
    hsr = [False] * mu.shape[0] #might not be possible without two hashes
    for i,sq,gi,sq_p1,gi_p1 in mt.values:
        if hs[i] or hsr[i]: #not absolutely sure.
            continue
        for i2,sq2,gi2,sq_m12,gi_m12 in mu.values:
            if hsr[i2] or hs[i2]:
                continue
            if (gi == gi2) and (gi_p1 == gi_m12) and (sq == sq_p1) and (sq2 == sq_m12): #reverse match pair
                b = eb_rev(mt,mu,i,i2)
                if abs(i2-i) == len(b)-1: #palindromes #placed before hash to allow internal blocks
                    continue
                #hs[i:i+len(b)] = [True]*len(b)
                #hs[i2-len(b):(i2-len(b))+len(b)] = [True]*len(b)
                hs[i] = [True]
                hs[i2-len(b) : i2] = [True]*len(b) 
                hsr[i:i+len(b)] = [True]*len(b)
                hsr[i2] = [True]
                t.append(('-',i,i2,b,len(b)))
                
    return pd.DataFrame(t)


#remove confined subblocks in either orientation
def clear_confined(df):
    t = []
    for i,r in df[df[4] > 2].sort_values(4, ascending = 0).iterrows():
        if r[0] == '-':
            l1,l2 = list(range(r[1],r[1]+r[4])), list(range(r[2]-r[4]+1, r[2]+1))
        else:
            l1,l2 = list(range(r[1],r[1]+r[4])), list(range(r[2], r[2]+r[4]))

        for i2,r2 in df[df[4] <= r[4]].iterrows():

            if r2[0] == '-':
                a1 = min(r2[1],r2[1]+r2[4]-1,r2[2],r2[2]-r2[4]+1)
                a2 = max(r2[1],r2[1]+r2[4]-1,r2[2],r2[2]-r2[4]+1)
            else:
                a1 = min(r2[1],r2[1]+r2[4]-1,r2[2],r2[2]+r2[4]-1)
                a2 = max(r2[1],r2[1]+r2[4]-1,r2[2],r2[2]+r2[4]-1)

            if ((a1 in l1) and (a2 in l1)) or ((a1 in l2) and (a2 in l2)):
                t.append(i2)
    return df.drop(t)

def reformat_reverse(tf1):
    tf1.loc[tf1[0] == '-', 2] = tf1.loc[tf1[0] == '-'][2] - tf1.loc[tf1[0] == '-'][4] + 1
    tf1.loc[:,7] = tf1.loc[:,1]
    tf1.loc[tf1[1] > tf1[2], 3] = tf1.loc[tf1[1] > tf1[2], 3].apply(lambda x: x[::-1])
    tf1.loc[tf1[7] > tf1[2], 1] = tf1.loc[tf1[7] > tf1[2], 2]
    tf1.loc[tf1[7] > tf1[2], 2] = tf1.loc[tf1[7] > tf1[2], 7]
    tf1[5] = tf1[1]+tf1[4]-1
    tf1[6] = tf1[2]+tf1[4]-1
    return tf1

def self_overlaps(tf1):
    #return tf1[(tf1[6]-tf1[5]) < tf1[4]] #method 1 
    return tf1[(tf1[2] >= tf1[1]) & (tf1[2] <= tf1[5])] #method 2

def simple_stack(d1,tf1):
    t1 = pd.DataFrame(range(len(d1)))
    w = [1]*len(d1)
    for i,r in tf1.sort_values(4,ascending=0).iterrows():

        x = max(w[r[1]:r[1]+r[4]])   
        if r[0] == '+':
            t1.loc[range(r[1],r[1]+r[4]),x] = range(r[2],r[2]+r[4])
            for i2 in range(r[1],r[1]+r[4]):
                w[i2]+=1
        else:
            t1.loc[range(r[1],r[1]+r[4]),x] = range(r[2],r[2]+r[4])[::-1]
            for i2 in range(r[1],r[1]+r[4]):
                w[i2]+=1
                
    r1 = set(t1.loc[:,1:].values.flatten()[t1.loc[:,1:].values.flatten() > 0]) #entries that appear as stacks
    r2 = set(t1[t1.loc[:,1:].sum(axis = 1) == 0].index) #entries not having any stacks
    t2 = t1.loc[~t1.index.isin(r1.intersection(r2)),:]
    return [t1, t2]

def map_genes(d1, t2):
    t = {-1:-1}
    for i,r in prep(d1).iterrows():
        t[i] = r['gi']
    return t2.fillna(-1).applymap(lambda x: t[x]).replace(-1,np.nan)

