"""Exact optimization in a stated range-certificate model (Python stdlib only).

Not an ASIC area model and not an exact reachable-range analyzer.
For unsigned nonnegative circuits, load h certifies 0 <= z <= h*(p-1).
For signed centered circuits, load h certifies |z| <= h*(p-1)/2.
Every raw operation
must fit headroom H BEFORE any optional canonical reduction at that node.
"""
from dataclasses import dataclass
from itertools import combinations, product
import random


@dataclass(frozen=True)
class Node:
    name: str
    terms: tuple  # (parent id, nonzero integer coefficient)
    output: bool = False


@dataclass(frozen=True)
class Circuit:
    name: str
    inputs: int
    nodes: tuple
    outputs: tuple

    def validate(self):
        for k, v in enumerate(self.nodes, self.inputs):
            assert v.terms and all(0 <= j < k and a != 0 for j, a in v.terms)
        assert all(self.inputs <= j < self.inputs + len(self.nodes) for j in self.outputs)
        assert set(self.outputs) == {self.inputs+i for i,v in enumerate(self.nodes) if v.output}

    def matrix(self):
        rows = [tuple(int(i == j) for i in range(self.inputs)) for j in range(self.inputs)]
        for v in self.nodes:
            rows.append(tuple(sum(a*rows[j][k] for j,a in v.terms) for k in range(self.inputs)))
        return tuple(rows[j] for j in self.outputs)


def circuit(name, n, operations, out_names):
    ids = {f'x{i}':i for i in range(n)}
    nodes = []
    for label, terms in operations:
        assert label not in ids
        nodes.append(Node(label, tuple((ids[j],a) for j,a in terms)))
        ids[label] = n+len(nodes)-1
    outs=[]
    for i,label in enumerate(out_names):
        nodes.append(Node(f'out{i}', ((ids[label],1),), True))
        outs.append(n+len(nodes)-1)
    c=Circuit(name,n,tuple(nodes),tuple(outs))
    c.validate()
    return c


def hl_m4():
    return circuit('HorizenLabs_M4',4,[
        ('s0',[('x0',1),('x1',1)]), ('s1',[('x2',1),('x3',1)]),
        ('d1',[('x1',2)]), ('t2',[('d1',1),('s1',1)]),
        ('d3',[('x3',2)]), ('t3',[('d3',1),('s0',1)]),
        ('d_s1',[('s1',2)]), ('q_s1',[('d_s1',2)]),
        ('t4',[('q_s1',1),('t3',1)]),
        ('d_s0',[('s0',2)]), ('q_s0',[('d_s0',2)]),
        ('t5',[('q_s0',1),('t2',1)]),
        ('t6',[('t3',1),('t5',1)]), ('t7',[('t2',1),('t4',1)]),
    ],['t6','t5','t7','t4'])


def p3_m4():
    return circuit('Plonky3_M4',4,[
        ('s01',[('x0',1),('x1',1)]), ('s23',[('x2',1),('x3',1)]),
        ('s',[('s01',1),('s23',1)]),
        ('a',[('s',1),('x1',1)]), ('b',[('s',1),('x3',1)]),
        ('d0',[('x0',2)]), ('d2',[('x2',2)]),
        ('y3',[('b',1),('d0',1)]), ('y1',[('a',1),('d2',1)]),
        ('y0',[('a',1),('s01',1)]), ('y2',[('b',1),('s23',1)]),
    ],['y0','y1','y2','y3'])


def signed_m4():
    """A norm-five MDS lift; 12 additions/subtractions, peak certificate 5."""
    return circuit('Signed_A49',4,[
        ('s01',[('x0',1),('x1',1)]),
        ('d01',[('x1',1),('x0',-1)]),
        ('s23',[('x2',1),('x3',1)]),
        ('d23',[('x3',1),('x2',-1)]),
        ('a',[('s01',1),('x0',1)]),
        ('b',[('d01',1),('x1',1)]),
        ('c',[('x2',1),('d23',-1)]),
        ('d',[('s23',1),('x3',1)]),
        ('y0',[('a',1),('s23',1)]),
        ('y1',[('b',1),('d23',1)]),
        ('y2',[('c',1),('d01',1)]),
        ('y3',[('d',1),('s01',-1)]),
    ],['y0','y1','y2','y3'])


def legal_choices(v, raw, H):
    if raw > H:
        return ()
    if v.output:
        return ((1, int(raw > 1)),)
    if raw == 1:
        return ((raw, 0),)
    return ((raw, 0), (1, 1))


def assess(c, H, selected):
    selected=set(selected)
    loads=[1]*c.inputs
    details=[]
    cost=0
    for i,v in enumerate(c.nodes,c.inputs):
        raw=sum(loads[j]*abs(a) for j,a in v.terms)
        if raw>H:
            return None
        red=(raw>1) if v.output else i in selected
        post=1 if red else raw
        if red:
            cost+=1
        loads.append(post)
        details.append({'id':i,'name':v.name,'raw_load':raw,'post_load':post,'reduce':red})
    return {'reductions':cost,'loads':loads,'details':details}


def optimize(c,H, dominance=False):
    """Frontier DP; exact for the unit-cost compositional load model.

    Identical live-load vectors merge. Optional Pareto dominance is sound
    because transitions and unit reduction costs are monotone in loads.
    The witness is a bit mask of selected internal nodes.
    """
    c.validate()
    last=[-1]*(c.inputs+len(c.nodes))
    for i,v in enumerate(c.nodes,c.inputs):
        for j,a in v.terms:
            last[j]=max(last[j],i)
    live=[j for j in range(c.inputs) if last[j]>=c.inputs]
    states={tuple(1 for _ in live):(0,0)}
    peak=len(states)
    max_frontier=len(live)
    transitions=0
    for i,v in enumerate(c.nodes,c.inputs):
        index={j:k for k,j in enumerate(live)}
        newlive=[j for j in live if last[j]>i]
        keep=[index[j] for j in newlive]
        retain=last[i]>i
        if retain:
            newlive.append(i)
        nxt={}
        for state,(cost,mask) in states.items():
            raw=sum(abs(a)*state[index[j]] for j,a in v.terms)
            for post,red in legal_choices(v,raw,H):
                transitions+=1
                key=tuple(state[k] for k in keep) + ((post,) if retain else ())
                val=(cost+red, mask | ((1<<i) if red and not v.output else 0))
                if key not in nxt or val<nxt[key]:
                    nxt[key]=val
        if dominance and len(nxt)>1:
            kept={}
            for key,val in sorted(nxt.items(), key=lambda kv:(kv[1][0],sum(kv[0]),kv[0])):
                if not any(oldval[0]<=val[0] and all(a<=b for a,b in zip(old,key))
                           for old,oldval in kept.items()):
                    kept[key]=val
            nxt=kept
        live,states=newlive,nxt
        peak=max(peak,len(states))
        max_frontier=max(max_frontier,len(live))
        if not states:
            return None
    assert live==[] and list(states)==[()]
    cost,mask=states[()]
    selected=[i for i in range(c.inputs,c.inputs+len(c.nodes)) if mask>>i&1]
    result=assess(c,H,selected)
    assert result and result['reductions']==cost
    result.update(selected=selected,peak_states=peak,max_frontier=max_frontier,transitions=transitions)
    return result


def brute(c,H):
    """Independent exhaustive enumeration of optional normalization subsets."""
    optional=[i for i,v in enumerate(c.nodes,c.inputs) if not v.output]
    best=None
    for mask in range(1<<len(optional)):
        pick=[i for b,i in enumerate(optional) if mask>>b&1]
        r=assess(c,H,pick)
        if r is not None and (best is None or r['reductions']<best['reductions']):
            best=r
            best['selected']=pick
    return best


def evaluate(c,p,selected,x,centered=False):
    def canonical(z):
        r=z%p
        return r-p if centered and r>p//2 else r
    vals=list(x)
    for i,v in enumerate(c.nodes,c.inputs):
        raw=sum(a*vals[j] for j,a in v.terms)
        vals.append(canonical(raw) if v.output or i in selected else raw)
    return [vals[j] for j in c.outputs],vals


def vertex_cover_circuit(n,edges):
    ops=[]
    for i in range(n):
        ops.append((f'v{i}',[(f'x{2*i}',1),(f'x{2*i+1}',1)]))
    for i,(u,v) in enumerate(edges):
        ops.append((f'e{i}',[(f'v{u}',1),(f'v{v}',1)]))
    return circuit('VC_gadget',2*n,ops,[f'e{i}' for i in range(len(edges))])


def min_vertex_cover(n,edges):
    for k in range(n+1):
        for pick in combinations(range(n),k):
            s=set(pick)
            if all(u in s or v in s for u,v in edges):
                return k


def det_integer(a):
    """Fraction-free Bareiss determinant, with row pivoting."""
    n=len(a)
    if n==1:
        return a[0][0]
    a=[list(r) for r in a]
    sign=1
    prev=1
    for k in range(n-1):
        if a[k][k]==0:
            row=next((i for i in range(k+1,n) if a[i][k]),None)
            if row is None:return 0
            a[k],a[row]=a[row],a[k]
            sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                val=a[i][j]*pivot-a[i][k]*a[k][j]
                assert val%prev==0
                a[i][j]=val//prev
        for i in range(k+1,n):a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]


def minor_certificate(matrix):
    n=len(matrix)
    return [dict(rows=list(rs),cols=list(cs),det=det_integer([[matrix[r][s] for s in cs] for r in rs]))
            for k in range(1,n+1) for rs in combinations(range(n),k) for cs in combinations(range(n),k)]


def cse_circuit(matrix,seed):
    """Standard positive common-subexpression heuristic, not claimed new.

    A randomized tie-break yields a finite catalog of exactly equivalent DAGs.
    """
    rng=random.Random(seed)
    n=len(matrix[0]); nextid=n
    rows=[{j:a for j,a in enumerate(row) if a} for row in matrix]
    ops=[]
    while any(sum(row.values())>1 for row in rows):
        scores={}
        for row in rows:
            keys=sorted(row)
            for i,a in enumerate(keys):
                for b in keys[i:]:
                    freq=row[a]//2 if a==b else min(row[a],row[b])
                    if freq:scores[a,b]=scores.get((a,b),0)+freq
        best=max(scores.values())
        a,b=rng.choice(sorted(k for k,v in scores.items() if v==best))
        terms=((a,2),) if a==b else ((a,1),(b,1))
        ops.append(Node(f'g{nextid}',terms))
        for row in rows:
            freq=row.get(a,0)//2 if a==b else min(row.get(a,0),row.get(b,0))
            if not freq:continue
            row[a]-=freq;row[b]-=freq
            row[nextid]=freq
            for k in [a,b]:
                if row.get(k)==0:row.pop(k,None)
        nextid+=1
    outputs=[]
    for i,row in enumerate(rows):
        assert sum(row.values())==1
        j=next(iter(row))
        outputs.append(nextid)
        ops.append(Node(f'out{i}',((j,1),),True));nextid+=1
    c=Circuit(f'CSE_{seed}',n,tuple(ops),tuple(outputs))
    c.validate();assert c.matrix()==tuple(tuple(r) for r in matrix)
    return c
