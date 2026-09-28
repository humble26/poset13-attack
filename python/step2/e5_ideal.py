import itertools, random, sys
from fractions import Fraction

def closure(n, rel):
    less=[[False]*n for _ in range(n)]
    for (i,j) in rel: less[i][j]=True
    for k in range(n):
        for i in range(n):
            if less[i][k]:
                ri,rk=less[i],less[k]
                for j in range(n):
                    if rk[j]: ri[j]=True
    return less

def is_poset(less,n):
    for i in range(n):
        if less[i][i]: return False
        for j in range(i+1,n):
            if less[i][j] and less[j][i]: return False
    return True

def count_ext(less,n):
    dp=[0]*(1<<n); dp[0]=1
    for S in range(1<<n):
        v=dp[S]
        if not v: continue
        for x in range(n):
            if (S>>x)&1: continue
            ok=True
            for y in range(n):
                if less[y][x] and not (S>>y)&1: ok=False; break
            if ok: dp[S|1<<x]+=v
    return dp[-1]

def antichain_width(less,n):
    best=0
    for S in range(1,1<<n):
        c=bin(S).count('1')
        if c<=best: continue
        els=[x for x in range(n) if (S>>x)&1]
        if all(not(less[a][b] or less[b][a]) for a,b in itertools.combinations(els,2)):
            best=c
    return best

def induced(less, n, subset):
    """induced poset on subset (list), reindexed"""
    idx={v:k for k,v in enumerate(subset)}
    m=len(subset)
    rel=[(idx[i],idx[j]) for i in subset for j in subset if less[i][j]]
    return closure(m,rel), idx

def ideals(less,n):
    """all down-closed subsets (as bitmasks over 0..n-1)"""
    res=[]
    for S in range(1<<n):
        ok=True
        for i in range(n):
            if (S>>i)&1:
                for j in range(n):
                    if less[j][i] and not (S>>j)&1: ok=False; break
            if not ok: break
        if ok: res.append(S)
    return res

def pval(less,n,x,y):
    rel=[(i,j) for i in range(n) for j in range(n) if less[i][j]]
    return Fraction(count_ext(closure(n,rel+[(x,y)]),n), count_ext(less,n))

THIRD=Fraction(1,3)

def analyze(less,n,z):
    """returns dict with classification for P, fixed z"""
    Q=[v for v in range(n) if v!=z]
    Ql,idx=induced(less,n,Q)
    m=n-1
    D=[v for v in Q if less[v][z]]
    U=[v for v in Q if less[z][v]]
    base=count_ext(less,n)
    ps={}
    for x,y in itertools.combinations(range(n),2):
        if less[x][y] or less[y][x]: continue
        ps[(x,y)]=pval(less,n,x,y)
    # z-pairs
    zp=[]
    for v in Q:
        if (z,v) in ps: zp.append(ps[(z,v)])
        elif (v,z) in ps: zp.append(ps[(v,z)])
    qp=[]
    for (x,y),p in ps.items():
        if x!=z and y!=z: qp.append(p)
    bal=lambda p: THIRD<=p<=1-THIRD
    zpair_ok=any(bal(p) for p in zp)
    qpair_ok=any(bal(p) for p in qp)
    d=max([min(p,1-p) for p in ps.values()] or [Fraction(0)])
    best_z=max([min(p,1-p) for p in zp] or [Fraction(0)])
    return dict(target=(d>=THIRD), zpair=zpair_ok, qpair=qpair_ok,
                delta=d, best_z=best_z, zmax=(len(U)==0), nD=len(D), nU=len(U))

def lemmaD_check(less,n,z):
    """e(P) = sum over ideals I of Q with D<=I, IcapU=empty of e(I)e(QI)"""
    Q=[v for v in range(n) if v!=z]
    Ql,idx=induced(less,n,Q)
    m=n-1
    Dmask=0
    for v in Q:
        if less[v][z]: Dmask|=1<<idx[v]
    Umask=0
    for v in Q:
        if less[z][v]: Umask|=1<<idx[v]
    tot=0
    for S in ideals(Ql,m):
        if (S & Dmask)==Dmask and (S & Umask)==0:
            Ilist=[Q[k] for k in range(m) if (S>>k)&1]
            Clist=[Q[k] for k in range(m) if not (S>>k)&1]
            Il,_=induced(less,n,Ilist); Cl,_=induced(less,n,Clist)
            tot+=count_ext(Il,len(Ilist))*count_ext(Cl,len(Clist))
    return tot==count_ext(less,n)

random.seed(5)
print("=== Lemma D (ideal decomposition) verification ===")
ok=bad=0
for t in range(300):
    n=random.choice([4,5,6])
    rel=[(i,j) for i in range(n) for j in range(i+1,n) if random.random()<0.4]
    L=closure(n,rel)
    if not is_poset(L,n): continue
    z=random.randrange(n)
    if lemmaD_check(L,n,z): ok+=1
    else: bad+=1; print("  LEMMA D FAIL", n, z, rel)
print("Lemma D: %d ok, %d bad"%(ok,bad))

print("=== strip formula p_P(z,u) = sum_{I in box, u notin I} e(I)e(QI)/e(P) ===")
ok=bad=0
for t in range(300):
    n=random.choice([4,5,6])
    rel=[(i,j) for i in range(n) for j in range(i+1,n) if random.random()<0.4]
    L=closure(n,rel)
    if not is_poset(L,n): continue
    z=random.randrange(n)
    Q=[v for v in range(n) if v!=z]
    Ql,idx=induced(L,n,Q); m=n-1
    Dmask=sum(1<<idx[v] for v in Q if L[v][z])
    Umask=sum(1<<idx[v] for v in Q if L[z][v])
    inc=[v for v in Q if not(L[v][z] or L[z][v])]
    if not inc: continue
    u=random.choice(inc)
    tot=0
    for S in ideals(Ql,m):
        if (S & Dmask)==Dmask and (S & Umask)==0 and not (S>>idx[u])&1:
            Ilist=[Q[k] for k in range(m) if (S>>k)&1]
            Clist=[Q[k] for k in range(m) if not (S>>k)&1]
            Il,_=induced(L,n,Ilist); Cl,_=induced(L,n,Clist)
            tot+=count_ext(Il,len(Ilist))*count_ext(Cl,len(Clist))
    lhs=pval(L,n,z,u)  # P(z before u)
    rhs=Fraction(tot,count_ext(L,n))
    if lhs==rhs: ok+=1
    else: bad+=1; print("  STRIP FAIL", n, z, u, lhs, rhs)
print("strip formula: %d ok, %d bad"%(ok,bad))
