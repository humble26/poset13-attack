import itertools
from fractions import Fraction
from pathlib import Path
exec(open(Path(__file__).with_name("e5_ideal.py"), encoding="utf-8").read().split("random.seed(5)")[0])

third=Fraction(1,3); twothird=Fraction(2,3)

def corridor_data(less,n,z):
    """chain partition of Q (width<=2) via greedy 2-coloring of incomparability,
       then box (dA,dB,uA,uB), column/row masses of nu over box ideals."""
    Q=[v for v in range(n) if v!=z]
    Ql,idx=induced(less,n,Q); m=n-1
    # greedy chain partition: repeatedly take maximal chain
    inc=[[not(Ql[i][j] or Ql[j][i]) for j in range(m)] for i in range(m)]
    color=[-1]*m  # 0=A,1=B ; same color must be chain: use longest-path coloring
    # assign color by parity of longest chain ending at v? For width-2, 2-coloring chains:
    # simple greedy: process topologically, color[v] = 0 if some u<v with color 1 else 0... 
    # Use: color[v] = 0 if v minimal; else flip of any predecessor.
    order=sorted(range(m), key=lambda v: sum(1 for u in range(m) if Ql[u][v]))
    for v in order:
        preds=[u for u in range(m) if Ql[u][v]]
        color[v]= 1-color[preds[0]] if preds else 0
    A=[v for v in range(m) if color[v]==0]
    Bq=[v for v in range(m) if color[v]==1]
    # chain-sort within colors
    A.sort(key=lambda v: sum(1 for u in range(m) if Ql[u][v]))
    Bq.sort(key=lambda v: sum(1 for u in range(m) if Ql[u][v]))
    s,t=len(A),len(Bq)
    aidx={v:k for k,v in enumerate(A)}; bidx={v:k for k,v in enumerate(Bq)}
    # ideals of Ql as (i,j) grid points, E,F counts
    def E_of(i,j):
        Il=[A[k] for k in range(i)]+[Bq[k] for k in range(j)]
        Il2,_=induced(less,n,[Q[x] for x in Il])
        return count_ext(Il2,len(Il2))
    # total elements mapping: I subset as actual Q-elements
    Dmask=sum(1<<idx[v] for v in Q if less[v][z])
    Umask=sum(1<<idx[v] for v in Q if less[z][v])
    boxpts=[]
    for S in ideals(Ql,m):
        if (S & Dmask)==Dmask and (S & Umask)==0:
            i=sum(1 for k in range(s) if (S>>aidx[A[k]])&1)
            j=sum(1 for k in range(t) if (S>>bidx[Bq[k]])&1)
            Clist=[Q[x] for x in range(m) if not (S>>x)&1]
            Cl,_=induced(less,n,Clist)
            boxpts.append((i,j,E_of(i,j),count_ext(Cl,len(Clist))))
    eP=sum(E*F for (_,_,E,F) in boxpts)
    colmass={}
    rowmass={}
    for (i,j,E,F) in boxpts:
        colmass[i]=colmass.get(i,0)+E*F
        rowmass[j]=rowmass.get(j,0)+E*F
    return dict(eP=eP, colmass=colmass, rowmass=rowmass, s=s,t=t,
                dA=min((i for (i,j,_,_) in boxpts), default=None),
                uA=max((i for (i,j,_,_) in boxpts), default=None),
                dB=min((j for (i,j,_,_) in boxpts), default=None),
                uB=max((j for (i,j,_,_) in boxpts), default=None),
                npts=len(boxpts))

third_ = Fraction(1,3)
stats=dict(n=0, both_heavy=0, only_col=0, only_row=0, neither=0, weird=0)
examples=[]
min_heavy=Fraction(2)

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
            if r['zpair']: continue   # only zpair-fail cases
            stats['n']+=1
            cd=corridor_data(L,n,z)
            eP=cd['eP']
            # vertical strips: V(i0)=mass of columns < i0, i0 in (dA, uA]
            heavy_col=any(third_ < Fraction(cmass,eP) for cmass in cd['colmass'].values())
            heavy_row=any(third_ < Fraction(rmass,eP) for rmass in cd['rowmass'].values())
            for cmass in cd['colmass'].values():
                f=Fraction(cmass,eP)
                if f>third_ and f<min_heavy: min_heavy=f
            if heavy_col and heavy_row: stats['both_heavy']+=1
            elif heavy_col: stats['only_col']+=1
            elif heavy_row: stats['only_row']+=1
            else:
                stats['neither']+=1
                if len(examples)<6: examples.append((n,rel,z,cd))

print("zpair-fail cases:",stats['n'])
print(stats)
print("min heavy-line mass (col or row > 1/3 observed):",min_heavy)
for (n,rel,z,cd) in examples[:6]:
    print("neither-heavy example: n=%d z=%d s=%d t=%d box=(%s,%s,%s,%s) eP=%d"%(
        n,z,cd['s'],cd['t'],cd['dA'],cd['dB'],cd['uA'],cd['uB'],cd['eP']))
    print("   colmass:",dict(sorted(cd['colmass'].items())),"rowmass:",dict(sorted(cd['rowmass'].items())))
