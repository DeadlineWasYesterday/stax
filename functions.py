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
def clear_confined(df): #n2 iter. repeated indexing. slow.
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
    tf1 = tf1.loc[tf1.astype(str).drop_duplicates([2,3]).index,] #[a,a,a,..] structures creating overlaps
    return tf1

def clear_localized_subblocks(df): #n2 iter. repeated indexing. slow.
    t = []
    for i,r in df.iterrows():
        for i2,r2 in df[df[4] < r[4]].iterrows():
            if r[1] <= r2[1] and r[2] <= r2[2] and r[5] >= r2[5] and r[6] >= r2[6]:
                t.append(i2)    
    return df.drop(t)


def self_overlaps(tf1):
    #return tf1[(tf1[6]-tf1[5]) < tf1[4]] #method 1 
    return tf1[(tf1[2] >= tf1[1]) & (tf1[2] <= tf1[5])] #method 2

def simple_stack(d1,tf1):
    t1 = pd.DataFrame({0: prep(d1)['Sequence'], 1 : range(len(d1))})
    w = [2]*len(d1)
    for i,r in tf1.sort_values([4,2],ascending=[0,1]).iterrows():

        x = max(w[r[1]:r[1]+r[4]])   
        if r[0] == '+':
            t1.loc[range(r[1],r[1]+r[4]),x] = range(r[2],r[2]+r[4])
            for i2 in range(r[1],r[1]+r[4]):
                w[i2] = x+1
        else:
            t1.loc[range(r[1],r[1]+r[4]),x] = range(r[2],r[2]+r[4])[::-1]
            for i2 in range(r[1],r[1]+r[4]):
                w[i2] = x+1
                
    r1 = set(t1.loc[:,2:].values.flatten()[t1.loc[:,2:].values.flatten() > 0]) #entries that appear as stacks
    r2 = set(t1[t1.loc[:,2:].sum(axis = 1) == 0].index) #entries not having any stacks
    t2 = t1.loc[~t1.index.isin(r1.intersection(r2)),:]
    return [t1, t2]

def map_genes(d1, t2):
    t = {-1:-1}
    for i,r in prep(d1).iterrows():
        t[i] = r['gi']
    return t2.fillna(-1).applymap(lambda x: t[x]).replace(-1,np.nan)



#demarcate protrusion grids universally. 
t1bm1 = t1b.shift(-1)
t,t1 = [],[]
for i in range(len(t1b)):
    if sum(abs(t1bm1.iloc[i,2:] - t1b.iloc[i,2:]) == 1) != 0:
        #start
        t1.append(i)
    else:
        if len(t1)!= 0:
            t1.append(i)
            t.append(t1)
            t1 = []
t2 = []
t1b[t1b.shape[1]] = np.nan
for i in t:
    t1b.iloc[i,-1] = t.index(i)


def extract_blocks_from_grid(grid):
    #[(column, start_index, end_index, start_entry, end_entry, elements, length)]
    t = []
    for ci in grid.columns:
        starts = grid.loc[:,ci].dropna().loc[abs(grid.loc[:,ci].dropna() - grid.loc[:,ci].dropna().shift(1)) != 1]
        ends = grid.loc[:,ci].dropna().loc[abs(grid.loc[:,ci].dropna() - grid.loc[:,ci].dropna().shift(-1)) != 1]
        for si,ei,se,ee in zip(starts.index, ends.index, starts.values, ends.values):
            t.append((ci, si, ei, se, ee, list(range(int(se),int(ee)+1)), ei-si+1))
    return t


def compare_grids(source_grid, destination_grid):
    ap = (0,0,0,0)
    for c1 in source_grid.dropna(axis = 1, thresh = 1).columns:
        for c2 in destination_grid.dropna(axis = 1, thresh = 1).columns:
            si,ml = compare_columns(source_grid[c1], destination_grid[c2])
            if ml > ap[3]:
                ap = (c1,c2,si,ml)
        #t.append(ap)
    return ap #[(col1, col2, source_index, best_match_length)]
def compare_columns(column1, column2):
    sw = False
    if len(column1) < len(column2):
        sw = True
        t = column1#.reset_index(drop = 1)
        column1 = column2#.reset_index(drop = 1)
        column2 = t  
    mx = (0,0)
    for os in range(-len(column2),len(column1)):
        if sum(column1.shift(-os).iloc[:len(column2)].reset_index(drop = 1) == column2.reset_index(drop = 1)) >= mx[1]:
            mx = (os, sum(column1.shift(-os).iloc[:len(column2)].reset_index(drop = 1) == column2.reset_index(drop = 1)))
    if sw:
        return (-mx[0],mx[1]) #negative offset
    else:
        return mx
    

gl1,gl2 = [],[]
for gr in t1b[14].dropna().unique()[::-1]:
    m = set(t1b[t1b[14] == gr].loc[:,1:13].values.flatten()).intersection(set(vc2)) #matches #this test is not exclusive to source grids
    if len(m) != 0:
        gl1.append(gr)
    else:
        gl2.append(gr)




t1c = t1b.drop(14,axis=1)

def recur(sg, t1c):
    if sg in gl2:
        return 999,t1c
    
    gs = t1c[t1c[100] == sg].iloc[:,1:t1c.shape[1]-1]
    ap,d = (0,0,0,0),0
    for dg in t1c[100].dropna().unique()[::-1]:
        if sg == dg:
            continue
        gd = t1c[t1c[100] == dg].iloc[:,1:t1c.shape[1]-1]
        p = compare_grids(gs,gd)
        if p[3] > ap[3]:
            ap = p
            d = dg
    dg = d
    gd = t1c[t1c[100] == dg].iloc[:,1:t1c.shape[1]-1]
    c1,c2,os,ml = ap
    if ap[3] == 0:
        return 998,t1c
    else:
        #shorter grid always stacks on top of the longer grid
        if len(gs) > len(gd):
            temp1,temp2 = gs,sg
            gs,sg = gd,dg
            gd,dg = temp1,temp2
            c1,c2,os = c2,c1,-os
            #print(sg,dg,c1,c2,os)
        
        #update t1c 
        #x = gd.iloc[:,:99].dropna(axis = 1, thresh = 1).shape[1]#.columns.max()+1
        #sl = gs.iloc[:,:99].dropna(axis = 1, thresh = 1).shape[1]-2#.columns.max()-1
        x = gd.loc[:,:99].dropna(axis = 1, thresh = 1).columns.max()+1
        sl = gs.loc[:,:99].dropna(axis = 1, thresh = 1).columns.max()-1
        c = 0
        
        s1 = t1c.loc[:gd.index[0],:].shape[0]-1+c#-os+c
        
        if os < 0: #gd starts before gs
            s1 = s1-os #slice lower
            fgd = gd.iloc[abs(os):]
            fsm = min(len(gs), len(fgd))
            for fi in range(fsm):
                if gs[c1].iloc[fi] == fgd[c2].iloc[fi]:
                    s1=s1+fi #slice at first match location
                    break
        elif os > 0: #gd starts after gs
            fgs = gs.iloc[abs(os):]
            fsm = len(fgs) #min(len(gd), len(fgs))
            for fi in range(fsm):
                if fgs[c1].iloc[fi] == gd[c2].iloc[fi]:
                    s1=s1+fi #slice at first match location
                    break
            
        
        for i,r in gs.iterrows():
            print(sg,dg,i,s1,os)
            print(c1,c2,r[c1], t1c.iloc[s1,c2])
            if r[c1] != t1c.iloc[s1,c2]:
                top = t1c.iloc[:s1,:]
                insert = t1c.loc[i,:]
                insert.loc[100] = dg
                insert.index = [0]+list(insert.index+x-1)[1:-1]+[100]
                top.loc[1000000+insert.name,:] = insert
                bottom = t1c.iloc[s1:]
                t1c = pd.concat([top,bottom])
                t1c = t1c.drop(i,axis = 0)
            else:
                #print(x,sl,t1c.loc[t1c.iloc[s1,:].name,x:x+sl+1].shape,len(r[:sl+1]))
                ski = t1c.iloc[s1,:].name
                print(x,sl,ski)
                #print(gd)
                #print(gs)
                for soi in range(sl+1):
                    t1c.loc[ski,x+soi] = r.iloc[soi]
                #t1c.loc[t1c.iloc[s1,:].name,x:x+sl] = r[:sl+1].values
                t1c.loc[ski,100] = dg
                t1c = t1c.drop(i,axis = 0)
            s1+=1
            if sg < dg:
                s1-=1 #how to code like an idiot
            c+=1
            print(t1c.shape)
        t1c = t1c[list(range(t1c.shape[1]-1))+[100]]
        
        return recur(dg,t1c)

for sg in gl1:
    _,t1c = recur(sg,t1c)
    print(sg)