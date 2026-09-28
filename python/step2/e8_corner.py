import itertools
from fractions import Fraction
from pathlib import Path
exec(open(Path(__file__).with_name("e5_ideal.py"), encoding="utf-8").read().split("random.seed(5)")[0])
third=Fraction(1,3); twothird=Fraction(2,3)

def honest_partition(Ql, m):
    """find a partition of width-2 Ql into two chains by search (m<=5)."""
    # try all subsets A containing a minimal element; check A and comp are chains
    els=list(range(m))
    for r in range(1, m):
        for Asub in itertools.combinations(els, r):
            Aset=set(Asub); Bset=set(els)-Aset
            okA=all(Ql[a][b] or Ql[b][a] for a in Aset for b in Aset if a!=b) or len(Aset)<=1
            okB=all(Ql[a][b] or Ql[b][a] for a in Bset for b in Bset if a!=b) or len(Bset)<=1
            if okA and okB:
                A=sorted(Aset, key=lambda v: sum(1 for u in Aset if Ql[u][v]))
                B=sorted(Bset, key=lambda v: sum(1 for u in Bset if Ql[u][v]))
                return A,B
    return None,None

def corner_test(less,n,z):
    Q=[v for v in range(n) if v!=z]
    Ql,idx=induced(less,n,Q); m=n-1
    A,B=honest_partition(Ql,m)
    if A is None: return None
    # z-pair p-values with chain elements
    pzA=[]; pzB=[]
    for v in A: pzA.append(pval(less,n,z,Q[v]) if not(less[v and Q[v] or 0][z] or less[z][Q[v]]) else None)
    for v in B: pzB.append(pval(less,n,z,Q[v]) if not(less[Q[v]][z] or less[z][Q[v]]) else None)
    # recompute properly: element Q[v]; incomparable to z?
    pzA=[]; pzB=[]
    for v in A:
        u=Q[v]
        pzA.append(None if (less[u][z] or less[z][u]) else pval(less,n,z,u))
    for v in B:
        u=Q[v]
        pzB.append(None if (less[u][z] or less[z][u]) else pval(less,n,z,u))
    # if any z-pair balanced -> not a fail case (caller ensures); here find cuts:
    def cut(pz):
        # i* = last index with p < 1/3 among non-None; requires monotone & all unbalanced
        last_small=None
        for k,p in enumerate(pz):
            if p is None: continue
            if p<third: last_small=k
        return last_small
    ist=cut(pzA); jst=cut(pzB)
    if ist is None: ist=-1
    if jst is None: jst=-1
    res=[]
    # corner cells (ist,ist+1)x(jst,jst+1) mapped to elements; both +1 elements must be incomparable to z
    def cell(di,dj):
        ia, jb = ist+di, jst+dj
        if ia>=len(A) or jb>=len(B): return None
        if ia<0 or jb<0: return None
        ua=Q[A[ia]] if ia<len(A) else None
        ub=Q[B[jb]] if jb<len(B) else None
        if ua is None or ub is None: return None
        if less[ua][z] or less[z][ua] or less[ub][z] or less[z][ub]: return None
        if less[ua][ub] or less[ub][ua]: return ('comp',ua,ub)
        return ('inc',ua,ub,pval(less,n,ua,ub))
    for di in (0,1):
        for dj in (0,1):
            c=cell(di,dj)
            if c and c[0]=='inc':
                p=c[3]
                res.append((di,dj,c[1],c[2],p,min(p,1-p)))
    return res

stats=dict(fail=0, corner_bal=0, corner_not=0, comp_skip=0, no_cell=0)
bad_list=[]
for n in (5,6):
    npairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for mask in range(1<<len(npairs)):
        rel=[npairs[k] for k in range(len(npairs)) if (mask>>k)&1]
        L=closure(n,rel)
        if not is_poset(L,n): continue
        if antichain_width(L,n)<2: continue
        for z in range(n):
            Qz=[v for v in range(n) if v!=z]
            Ql,_=induced(L,n,Qz)
            if antichain_width(Ql,n-1)>2: continue
            r=analyze(L,n,z)
            if r['zpair']: continue
            stats['fail']+=1
            res=corner_test(L,n,z)
            if res is None: stats['no_cell']+=1; continue
            if not res: stats['no_cell']+=1; continue
            if any(minp>=third for (_,_,_,_,_,minp) in res):
                stats['corner_bal']+=1
            else:
                stats['corner_not']+=1
                if len(bad_list)<8: bad_list.append((n,rel,z,res))
print(stats)
for (n,rel,z,res) in bad_list:
    print("CORNER-NOT-BAL n=%d z=%d rel=%s"%(n,z,rel))
    for (di,dj,ua,ub,p,minp) in res:
        print("   cell(+%d,+%d) pair(%d,%d) p=%s min=%s"%(di,dj,ua,ub,p,minp))
