"""Revision-4 certificates for all three signed representatives.

Reachability is exhaustive. Circuit search is only a witness search: no global
arithmetic optimality is inferred for A7 or A17 from the search.
"""
from collections import Counter
from itertools import product, permutations
from pathlib import Path
import json
from classify import REPS
from catalog import l1, canonical, representations, trees
from reduction_model import circuit, optimize, brute
ROOT=Path(__file__).resolve().parents[1]


def a7_eleven():
    return circuit('A7_eleven',4,[
        ('w',[('x2',1),('x3',-1)]),
        ('t',[('x1',1),('x0',-1)]),
        ('a',[('t',1),('x1',1)]),('y1',[('a',1),('w',-1)]),
        ('q',[('x2',1),('w',1)]),('y2',[('t',1),('q',1)]),
        ('b',[('x0',1),('x1',1)]),('c',[('b',1),('x3',1)]),
        ('y3',[('c',1),('w',-1)]),
        ('d',[('x0',1),('q',1)]),('y0',[('d',1),('y3',1)]),
    ],['y0','y1','y2','y3'])


def a17_eleven():
    return circuit('A17_eleven',4,[
        ('w',[('x2',1),('x1',-1)]),
        ('u',[('x0',1),('w',1)]),('v',[('x0',1),('w',-1)]),
        ('a',[('x1',1),('x3',1)]),('b',[('x2',1),('x3',1)]),
        ('y1',[('a',1),('u',-1)]),('y2',[('b',1),('u',1)]),
        ('t',[('y2',1),('x1',1)]),('y0',[('t',1),('v',1)]),
        ('d3',[('x3',2)]),('y3',[('d3',1),('v',1)]),
    ],['y0','y1','y2','y3'])


def first_forms(A):
    universe=[w for w in product(range(-4,5),repeat=4)
              if 2<=l1(w)<=4 and next(x for x in w if x)>0]
    rows=[]
    for w in universe:
        reached={0};transitions=[]
        for mask in range(16):
            if mask not in reached:continue
            done=[j for j in range(4) if mask>>j&1]
            for i in range(4):
                if mask>>i&1:continue
                reps=representations(A[i],[w]+[A[j] for j in done])
                if reps:
                    reached.add(mask|1<<i)
                    transitions.append(dict(mask=mask,output=i,representations=[dict(coefficients=a,inputs=u) for a,u in reps]))
        rows.append(dict(w=w,norm=l1(w),reachable_masks=sorted(reached),transitions=transitions))
    return dict(candidate_count=len(rows),max_reachable_histogram=dict(Counter(max(m.bit_count() for m in r['reachable_masks']) for r in rows)),
                successful_by_norm=dict(Counter(r['norm'] for r in rows if 15 in r['reachable_masks'])),all_first_forms=rows)


def witness_search(A,record,beam=80):
    """Bounded beam search over sign-compatible residual expression trees.

    Only used to locate concrete witnesses. The search grammar/beam cannot
    establish an arithmetic lower bound over unrestricted integer programs.
    """
    w=tuple(record['w']);dim=9
    wc,ws=canonical(w+(0,)*5)
    initial=[(t,[(None,wc,t)]) for t in trees(wc)]
    best=None
    transitions={(r['mask'],r['output']):r['representations'] for r in record['transitions']}
    for order in permutations(range(4)):
        states=initial;mask=0
        for i in order:
            rr=transitions.get((mask,i),[])
            if not rr:states=[];break
            done=[j for j in range(4) if mask>>j&1]
            options=[]
            for rep in rr:
                a,u=rep['coefficients'],rep['inputs']
                v=list(u)+[a[0]]+[0]*4
                for j,b in zip(done,a[1:]):v[5+j]=b
                vc,sgn=canonical(tuple(v))
                for tree in trees(vc):options.append((tree,vc,sgn))
            nxt={}
            for previous,stages in states:
                d=dict(previous)
                for tree,vc,sgn in options:
                    if any(k in d and d[k]!=g for k,g in tree):continue
                    merged=previous|tree
                    if best is not None and len(merged)>best[0]:continue
                    if merged not in nxt:nxt[merged]=stages+[(i,vc,tree)]
            states=sorted(nxt.items(),key=lambda s:(len(s[0]),str(sorted(s[0]))))[:beam]
            mask|=1<<i
            if not states:break
        if states:
            candidate=min(states,key=lambda s:len(s[0]))
            if best is None or len(candidate[0])<best[0]:best=(len(candidate[0]),order,candidate[1])
    if best is None:return None
    count,order,stages=best
    # A generic perturbed evaluation orients wires without zero formal forms.
    # All output rows evaluate positively near the all-one input. Different
    # formal reset atoms get tiny distinct perturbations to avoid degeneracy.
    probe=(10**9,10**9+1,10**9+101,10**9+10201)
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    scale=10**12
    values=[a*scale for a in probe]+[dot(w,probe)*scale+1]+[dot(r,probe)*scale+10**(j+1) for j,r in enumerate(A)]
    ori=lambda v:1 if dot(v,values)>0 else -1
    units=[tuple(int(i==j) for i in range(dim)) for j in range(dim)]
    names={units[i]:f'x{i}' for i in range(4)}
    ops=[];outs={};reset=[]
    for output,root,tree in stages:
        for v,children in sorted(tree,key=lambda q:(l1(q[0]),q[0])):
            if v in names:continue
            assert dot(v,values)!=0
            terms={}
            for u,s in children:
                label=names[u];coef=ori(v)*s*ori(u)
                terms[label]=terms.get(label,0)+coef
            terms={k:a for k,a in terms.items() if a}
            assert terms and any(a>0 for a in terms.values())
            assert all(a in (-1,1,2) for a in terms.values())
            label=f'g{len(ops)}';ops.append((label,list(terms.items())));names[v]=label
        label=names[root]
        if output is None:
            assert ori(root)==ori(units[4])
            names[units[4]]=label
        else:
            names[units[5+output]]=label;outs[output]=label
        reset.append(label)
    c=circuit('residual_witness',4,ops,[outs[i] for i in range(4)])
    assert c.matrix()==A,(c.matrix(),A,w,order)
    r=optimize(c,4)
    if r is None or r['reductions']!=5:return None
    return c,r,dict(w=w,output_order=order,beam=beam,operations=len(ops))


def main():
    all_results={}
    for label,A in REPS.items():
        result=first_forms(A)
        all_results[str(label)]=result
        print(label,result['max_reachable_histogram'],result['successful_by_norm'],flush=True)
    (ROOT/'data/all_class_reachability.json').write_text(json.dumps(all_results,indent=2)+'\n')
    import sys
    if '--search' not in sys.argv:
        return
    best=None
    for rec in all_results['17']['all_first_forms']:
        if 15 not in rec['reachable_masks'] or rec['norm']!=2:continue
        found=witness_search(REPS[17],rec)
        if found:
            c,r,meta=found
            print('A17 candidate',meta,flush=True)
            if best is None or meta['operations']<best[2]['operations']:best=found
    assert best
    c,r,meta=best
    (ROOT/'data/a17_search_witness.json').write_text(json.dumps(dict(metadata=meta,nodes=[dict(name=v.name,terms=v.terms,output=v.output) for v in c.nodes],outputs=c.outputs,schedule=r),indent=2)+'\n')
    print('A17 best',meta,flush=True)
    print([(v.name,v.terms) for v in c.nodes],flush=True)

if __name__=='__main__':main()
