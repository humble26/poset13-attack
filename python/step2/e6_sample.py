import itertools, random
from fractions import Fraction
from pathlib import Path
exec(open(Path(__file__).with_name("e5_ideal.py"), encoding="utf-8").read().split("random.seed(5)")[0])
third=Fraction(1,3)

def rand_width2(m, rng):
    perm=list(range(m)); rng.shuffle(perm)
    h=m//2
    C1=perm[:h]; C2=perm[h:]
    rel={(C1[i],C1[i+1]) for i in range(len(C1)-1)}|{(C2[i],C2[i+1]) for i in range(len(C2)-1)}
    for a in C1:
        for b in C2:
            r=rng.random()
            if r<0.25: rel.add((a,b))
            elif r<0.5: rel.add((b,a))
    L=closure(m,rel)
    if not is_poset(L,m) or antichain_width(L,m)>2: return None
    return L

def add_point(L,m,rng):
    els=list(range(m))
    D=[v for v in els if rng.random()<0.3]
    U=[v for v in els if v not in D and rng.random()<0.3]
    if any(L[b][a] for a in D for b in U): return None
    rel=[(i,j) for i in range(m) for j in range(m) if L[i][j]]
    rel+=[(a,m) for a in D]+[(m,b) for b in U]
    P=closure(m+1,rel)
    if not is_poset(P,m+1): return None
    return P

rng=random.Random(2026)
for n in (7,8,9):
    stats=dict(total=0,target_fail=0,zpair_fail=0,both_fail=0)
    min_d=Fraction(2); found=0; trials=0
    while found<1500 and trials<40000:
        trials+=1
        Q=rand_width2(n-1,rng)
        if Q is None: continue
        P=add_point(Q,n-1,rng)
        if P is None: continue
        if antichain_width(P,n)<2: continue  # chain: conjecture vacuous
        z=n-1
        r=analyze(P,n,z)
        found+=1; stats['total']+=1
        if not r['target']: stats['target_fail']+=1; print("TARGET FAIL n=%d"%n, P and "see rel")
        if not r['zpair']:
            stats['zpair_fail']+=1
            if not r['qpair']: stats['both_fail']+=1
        if r['delta']<min_d: min_d=r['delta']
    print("n=%d found=%d min_delta=%s stats=%s"%(n,found,min_d,stats))
