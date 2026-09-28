import itertools
from fractions import Fraction

def closure(n, rel):
    less=[[False]*n for _ in range(n)]
    for (i,j) in rel: less[i][j]=True
    for k in range(n+2):
        for i in range(n):
            for j in range(n):
                if less[i][j]: continue
                if any(less[i][k2] and less[k2][j] for k2 in range(n)): less[i][j]=True
    return less

def is_poset_ok(less,n):
    for i in range(n):
        if less[i][i]: return False
        for j in range(n):
            if i!=j and less[i][j] and less[j][i]: return False
    return True

def trans_close_ok(less,n):
    for i in range(n):
        for j in range(n):
            if less[i][j] and any(less[j][k] and not less[i][k] for k in range(n)): return False
    return True

def exts(less,n):
    res=[]
    def rec(cur,avail):
        if not avail: res.append(tuple(cur)); return
        for x in avail:
            if not any(less[y][x] for y in avail): rec(cur+[x],[a for a in avail if a!=x])
    rec([],list(range(n)))
    return res

def e_(less,n): return len(exts(less,n))

def add_rel(less,n,x,y):
    def lt2(i,j): return less[i][j] or (i==x and j==y)
    # emulate Lean: closureStep on non-transitive base with n+2 rounds
    rel2=[(i,j) for i in range(n) for j in range(n) if lt2(i,j)]
    return closure(n,rel2)

def check_all_on(n, mask):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    rel=[pairs[k] for k in range(len(pairs)) if (mask>>k)&1]
    less=closure(n,rel)
    if not is_poset_ok(less,n) or not trans_close_ok(less,n): return None
    fails=[]
    # checkA
    maxl=[z for z in range(n) if not any(less[z][y] for y in range(n))]
    if e_(less,n)!=sum(e_([[less[i][j] for j in range(n) if j!=z] for i in range(n) if i!=z],n-1) for z in maxl):
        fails.append('A')
    # checkB: ORDERED pairs like Lean
    for x in range(n):
        for y in range(n):
            if x==y or less[x][y] or less[y][x]: continue
            crit=all(w==x or w==y or ((not less[w][x]) or less[w][y]) and ((not less[x][w]) or less[y][w]) for w in range(n))
            if crit:
                extsP=exts(less,n)
                if not any(any(M[k]==x and M[k+1]==y for k in range(n-1)) for M in extsP):
                    fails.append(('B',x,y))
    # checkConj
    found=False
    for x in range(n):
        for y in range(n):
            if x==y or less[x][y] or less[y][x]: continue
            if e_(less,n)<=3*e_(add_rel(less,n,x,y),n) and e_(less,n)<=3*e_(add_rel(less,n,y,x),n):
                found=True; break
        if found: break
    if not found: fails.append('CONJ')
    # cGuard/cCheck
    def down(z,i): return i if i<z else i-1
    for z in range(n):
        for x in range(n):
            for y in range(n):
                if x==z or y==z or x==y or less[x][y] or less[y][x]: continue
                guard=is_poset_ok(add_rel(less,n,x,y),n)
                if not guard: continue
                Q=[[less[i][j] for j in range(n) if j!=z] for i in range(n) if i!=z]
                below=[down(z,v) for v in range(n) if less[v][z]]
                above=[down(z,v) for v in range(n) if less[z][v]]
                ms=exts(Q,n-1); sx,sy=down(z,x),down(z,y)
                def slot_count(M):
                    c=0
                    for k in range(len(M)+1):
                        if all(v in M[:k] for v in below) and all(v in M[k:] for v in above): c+=1
                    return c
                den=sum(slot_count(M) for M in ms)
                num=sum(slot_count(M) for M in ms if M.index(sx)<M.index(sy))
                if den*e_(add_rel(less,n,x,y),n)!=num*e_(less,n):
                    fails.append(('C',z,x,y))
    # pGuard/pCheck
    for z in range(n):
        comp=all(w==z or less[w][z] or less[z][w] for w in range(n))
        if not comp: continue
        Q=[[less[i][j] for j in range(n) if j!=z] for i in range(n) if i!=z]
        if not is_poset_ok(Q,n-1): continue
        if e_(less,n)!=e_(Q,n-1): fails.append(('P1a',z))
        for x in range(n):
            for y in range(n):
                if x==z or y==z or x==y or less[x][y] or less[y][x]: continue
                if e_(add_rel(less,n,x,y),n)!=e_(add_rel(Q,n-1,down(z,x),down(z,y)),n-1):
                    fails.append(('P1b',z,x,y))
    return fails

total_posets=0
for mask in range(64):
    f=check_all_on(4,mask)
    if f is None: continue
    total_posets+=1
    if f:
        print("FAIL mask=%d rel=%s fails=%s"%(mask,[(i,j) for i in range(4) for j in range(i+1,4) if (mask>>[(i,j) for i in range(4) for j in range(i+1,4)].index((i,j)))&1],f[:6]))
print("posets checked:",total_posets)
