"""Focused verification for symmetric lower bounds and unit-cost tripling.

Existing classification and catalog data are unchanged.  This driver checks the
new primitive model and both five-normalization witnesses, plus the local cone
counts needed for the tripling corollary.
"""
from functools import lru_cache
from pathlib import Path
from itertools import product
import json
from reduction_model import circuit, p3_m4, signed_m4, optimize, brute, evaluate, assess
from catalog import canceling_fifteen, subforms, canonical, l1
ROOT=Path(__file__).resolve().parents[1]


def tripling_fourteen():
    return circuit('A49_unit_tripling_14',4,[
        ('s12',[('x1',1),('x2',1)]),
        ('s123',[('s12',1),('x3',1)]),
        ('w',[('x0',1),('s123',-1)]),
        ('t0',[('x0',3)]),('y0',[('t0',1),('w',-1)]),
        ('d2',[('x2',2)]),('t1',[('x1',1),('d2',-1)]),
        ('y1',[('t1',1),('w',-1)]),
        ('d3',[('x3',2)]),('t2',[('x2',1),('d3',-1)]),
        ('y2',[('t2',1),('w',-1)]),
        ('d1',[('x1',2)]),('t3',[('x3',1),('d1',-1)]),
        ('y3',[('t3',1),('w',-1)]),
    ],['y0','y1','y2','y3'])


@lru_cache(None)
def cone_cost(v,tripling):
    # Sign-compatible expression trees suffice after the tight-load proof.
    # Dedicated doubling/tripling options account for reusing a child wire.
    if l1(v)==1:return 0
    costs=[]
    for k in ([2,3] if tripling else [2]):
        if all(x%k==0 for x in v):costs.append(1+cone_cost(tuple(x//k for x in v),tripling))
    for a in subforms(v):
        b=tuple(x-y for x,y in zip(v,a))
        if any(b):costs.append(1+cone_cost(canonical(a)[0],tripling)+cone_cost(canonical(b)[0],tripling))
    assert costs
    return min(costs)


def main():
    data=json.loads((ROOT/'data/five_reduction_certificate.json').read_text())
    cones=[]
    for row in data['successful_forms']:
        w=tuple(row['w'])
        vectors=[tuple(o['input_coefficients'])+(o['w_coefficient'],) for o in row['outputs']]
        patterns=sorted(tuple(sorted(abs(a) for a in v if a)) for v in vectors)
        assert patterns==[(1,1,2)]*3+[(1,3)]
        sets=[{canonical(v)[0] for v in subforms(t) if l1(v)>=2} for t in vectors]
        pre={canonical(v+(0,))[0] for v in subforms(w) if 2<=l1(v)<4}
        assert all(not(pre&s) for s in sets)
        assert all(not(sets[i]&sets[j]) for i in range(4) for j in range(i))
        costs={}
        for triple in [False,True]:
            initial=cone_cost(canonical(w)[0],triple)
            outputs=[cone_cost(canonical(v)[0],triple) for v in vectors]
            total=initial+sum(outputs)
            assert initial==3 and total==(14 if triple else 15)
            costs['tripling' if triple else 'binary']=dict(initial=initial,outputs=outputs,total=total)
        cones.append(dict(w=w,costs=costs))
    results=[];semantic=0
    for c in [p3_m4(),tripling_fourteen()]:
        assert c.matrix()==(p3_m4().matrix() if c.name=='Plonky3_M4' else signed_m4().matrix())
        n=sum(not v.output for v in c.nodes)
        expected=11 if c.name=='Plonky3_M4' else 14
        assert n==expected
        norms=[l1(r) for r in c.matrix()]
        assert all(x>4 for x in norms)
        best=optimize(c,4);check=brute(c,4)
        assert best['reductions']==check['reductions']==5
        replay=assess(c,4,best['selected'])
        assert replay['reductions']==5
        assert sum(x['reduce'] for x in replay['details'] if x['raw_load']>1)==5
        for x in product(range(-6,7),repeat=4):
            y,values=evaluate(c,13,best['selected'],x,True)
            assert [t%13 for t in y]==[sum(a*b for a,b in zip(r,x))%13 for r in c.matrix()]
            for i,v in enumerate(c.nodes,c.inputs):
                raw=sum(a*values[j] for j,a in v.terms)
                assert abs(raw)<=6*replay['details'][i-c.inputs]['raw_load']
                assert abs(values[i])<=6*replay['loads'][i]
            semantic+=1
        results.append(dict(circuit=c.name,operations=n,row_norms=norms,normalizations=5,
                            schedule=best,nodes=[dict(name=v.name,terms=v.terms,output=v.output) for v in c.nodes],outputs=c.outputs))
    out=dict(capacity=4,cone_certificates=cones,witnesses=results,
             exhaustive_subset_comparisons=2,centered_semantic_inputs=semantic,
             common_normalization_lower_bound=5,
             P3_arithmetic_upper_bound_at_five=11,
             A49_binary_arithmetic_minimum_at_five=15,
             A49_tripling_arithmetic_minimum_at_five=14,
             A49_binary_budget_at_most_14_normalization_optimum=6)
    (ROOT/'data/revision3_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['cone_certificates','witnesses']},indent=2))

if __name__=='__main__':main()
