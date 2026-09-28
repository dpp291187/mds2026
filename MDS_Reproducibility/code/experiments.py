"""Deterministic additional experiments; all numbers are certificate-model results."""
from pathlib import Path
import json, random, platform, time, statistics
from reduction_model import hl_m4,p3_m4,signed_m4,optimize,Circuit,Node,det_integer,minor_certificate
from classify import prime_divisors

ROOT=Path(__file__).resolve().parents[1]

def reorder(c,seed):
    rng=random.Random(seed); old=list(range(c.inputs,c.inputs+len(c.nodes)))
    done=set(range(c.inputs)); order=[]
    while old:
        ready=[j for j in old if all(k in done for k,a in c.nodes[j-c.inputs].terms)]
        j=rng.choice(ready);order.append(j);old.remove(j);done.add(j)
    ids={j:j for j in range(c.inputs)}
    ids.update({j:c.inputs+k for k,j in enumerate(order)})
    ns=tuple(Node(c.nodes[j-c.inputs].name,tuple((ids[k],a) for k,a in c.nodes[j-c.inputs].terms),c.nodes[j-c.inputs].output) for j in order)
    r=Circuit(c.name,c.inputs,ns,tuple(ids[j] for j in c.outputs));r.validate()
    assert r.matrix()==c.matrix()
    return r

def main():
    circuits=[hl_m4(),p3_m4(),signed_m4()]
    grid=[]; orders=[]; width=[]; summary={}
    for c in circuits:
        M=c.matrix(); depth=[0]*c.inputs
        for v in c.nodes:
            depth.append(max(depth[j] for j,a in v.terms)+int(not v.output))
        summary[c.name]=dict(L=max(sum(abs(a) for a in row) for row in M),
                            det=det_integer(M),
                            bad_primes=sorted({p for m in minor_certificate(M) for p in prime_divisors(m['det'])}),
                            depth=max(depth),operations=sum(not v.output for v in c.nodes))
        for H in range(2,17):
            p=optimize(c,H);d=optimize(c,H,True)
            assert p['reductions']==d['reductions']
            grid.append(dict(circuit=c.name,H=H,cost=p['reductions'],plain_peak=p['peak_states'],dominance_peak=d['peak_states'],plain_transitions=p['transitions'],dominance_transitions=d['transitions'],frontier=p['max_frontier']))
        for H in [3,5,8]:
            base=optimize(c,H)['reductions']; records=[]
            for seed in range(100):
                ordered=reorder(c,20260926+seed)
                r=optimize(ordered,H); dom=optimize(ordered,H,True)
                assert r['reductions']==dom['reductions']==base
                record={k:r[k] for k in ['max_frontier','peak_states','transitions']}
                record.update(dominance_peak=dom['peak_states'],dominance_transitions=dom['transitions'])
                records.append(record)
            orders.append(dict(circuit=c.name,H=H,cost=base,seeds=100,records=records))
    # Integer-only width computation: 1 + ceil(log2(q L + 1)).
    for p in [13,65521,2013265921,2147483647,18446744069414584321]:
        q=(p-1)//2
        for L in [5,7,16]:
            width.append(dict(p=p,q=q,L=L,width=1+(q*L).bit_length()))
    out=dict(environment={'python':platform.python_version(),'platform':platform.platform()},grid=grid,order_experiments=orders,widths=width)
    ROOT.joinpath('data/extended_results.json').write_text(json.dumps(out,indent=2)+'\n')
    ROOT.joinpath('data/circuit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('45 capacity cases; 45 dominance equalities; 900 valid-order invariance checks.')
    for e in orders:
        r=e['records'];print(e['circuit'],e['H'],'frontier',min(a['max_frontier'] for a in r),max(a['max_frontier'] for a in r),'states',min(a['peak_states'] for a in r),max(a['peak_states'] for a in r))
    print('widths',width)

if __name__=='__main__':main()
