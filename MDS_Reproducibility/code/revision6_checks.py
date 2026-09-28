"""Focused finite checks for the revised word-size and A49 cone arguments.

Run this file directly in PyCharm, or run `python code/revision6_checks.py`
from the project root. Python standard library only; working-directory
independent. Only data/revision6_checks.json is written.

The BN254 scalar modulus is the specified standard prime. This script certifies
integer widths, minor invertibility and compositional node bounds at that exact
modulus; it does not claim a new primality proof or a measured implementation
speedup. The A49 finite checks support the explicitly stated tight-load cone
lemma, not a general inference from four input occurrences to three gates.
"""
from collections import Counter
from itertools import combinations, product
from math import gcd
from pathlib import Path
import hashlib
import json

from catalog import (canonical, five_reduction_certificate, l1,
                     representations, subforms, trees)
from classify import A49, A7, A17, B
from class_extension import a7_eleven, a17_eleven
from reduction_model import circuit, minor_certificate, optimize
from revision5_checks import OBSERVABLE_A7, observable_a7_eleven
from reduction_model import p3_m4

ROOT = Path(__file__).resolve().parents[1]
BN254_SCALAR_P = 21888242871839275222246405745257275088548364400416034343698204186575808495617


def centered_word_capacity(p, width):
    """Largest integer H with H*(p-1)/2 <= 2^(width-1)-1."""
    assert p > 2 and p % 2 and width >= 2
    q = (p-1)//2
    return ((1 << (width-1))-1)//q


def centered_signed_width(p, load):
    """Least W accommodating the symmetric certificate [-load*q,load*q]."""
    assert load >= 1
    peak = load*((p-1)//2)
    width = 1+peak.bit_length()
    assert peak <= (1 << (width-1))-1
    assert width == 2 or peak > (1 << (width-2))-1
    return width


def word_size_checks():
    p = BN254_SCALAR_P
    q = (p-1)//2
    maximum = (1 << 255)-1
    assert centered_word_capacity(p, 256) == 5
    assert 5*q <= maximum < 6*q
    widths = {str(L): centered_signed_width(p, L) for L in (5, 7, 16)}
    assert widths == {'5':256, '7':257, '16':258}
    small = []
    fields = [('BabyBear', 2013265921), ('Mersenne31', 2147483647),
              ('KoalaBear', 2130706433)]
    for name, prime in fields:
        rows = []
        for width, expected in [(32, 2), (33, 4), (34, 8)]:
            H = centered_word_capacity(prime, width)
            sq = (prime-1)//2
            cap = (1 << (width-1))-1
            assert H == expected and H*sq <= cap < (H+1)*sq
            rows.append(dict(signed_width_bits=width, capacity_H=H,
                             centered_input_peak=sq, signed_positive_limit=cap))
        small.append(dict(name=name, prime=prime, word_sizes=rows))
    return dict(BN254_scalar=dict(prime=p, centered_input_peak=q,
                                 signed_width_bits=256, signed_positive_limit=maximum,
                                 capacity_H=5, norm5_peak=5*q, norm6_peak=6*q,
                                 signed_one_shot_widths=widths,
                                 primality_scope='Specified standard BN254 scalar prime; no new primality certificate claimed'),
                small_fields=small, small_field_capacity_checks=9,
                one_shot_width_scope='Exact widths for one complete row-dot-product range; not a limb implementation cost')


def modulus_and_schedule_checks():
    p = BN254_SCALAR_P
    q = (p-1)//2
    maximum = (1 << 255)-1
    mats = [('A49', A49), ('A7', A7), ('A17', A17), ('B', B),
            ('observable_A7', OBSERVABLE_A7), ('P3', p3_m4().matrix())]
    matrix_rows = []
    for name, A in mats:
        cert = minor_certificate(A)
        assert len(cert) == 69
        assert all(gcd(m['det'], p) == 1 for m in cert)
        matrix_rows.append(dict(name=name, row_norms=[l1(row) for row in A],
                                invertible_minors=len(cert),
                                distinct_minor_magnitudes=sorted({abs(m['det']) for m in cert}),
                                MDS_over_specified_prime=True))
    schedules = []
    for c, A, expected in [(a7_eleven(), A7, 4), (a17_eleven(), A17, 4),
                           (observable_a7_eleven(), OBSERVABLE_A7, 4),
                           (p3_m4(), p3_m4().matrix(), 5)]:
        assert c.matrix() == A
        operations = sum(not n.output for n in c.nodes)
        assert operations == 11
        schedule = optimize(c, 5)
        assert schedule and schedule['reductions'] == expected
        node_bounds = []
        for node, detail in zip(c.nodes, schedule['details']):
            raw_peak = detail['raw_load']*q
            retained_peak = detail['post_load']*q
            assert raw_peak <= maximum and retained_peak <= maximum
            node_bounds.append(dict(name=node.name, output_sink=node.output,
                                    raw_load=detail['raw_load'], raw_absolute_bound=raw_peak,
                                    raw_signed_width=centered_signed_width(p, detail['raw_load']),
                                    retained_load=detail['post_load'], retained_absolute_bound=retained_peak,
                                    retained_signed_width=centered_signed_width(p, detail['post_load']),
                                    normalize=detail['reduce']))
        assert len(node_bounds) == 15
        schedules.append(dict(circuit=c.name, matrix=A, operations=operations,
                              capacity_H=5, normalizations=expected,
                              selected_internal_nodes=schedule['selected'],
                              exact_linear_map_verified=True,
                              every_raw_and_retained_bound_fits_signed_256=True,
                              node_bounds=node_bounds))
    return dict(prime=p, minor_checks=sum(r['invertible_minors'] for r in matrix_rows),
                matrices=matrix_rows, schedules=schedules,
                normalization_count_comparisons=4, raw_node_bound_checks=60,
                retained_node_bound_checks=60,
                arithmetic_count_scope='Four supplied eleven-operation witnesses; not global arithmetic minimality',
                capacity_scope='Centered signed storage with no fused modular/scalar operation assumption')


def sign_canonical_or_zero(vector):
    if not any(vector):
        return tuple(vector)
    return canonical(tuple(vector))[0]


def next_gate_results(available):
    """All binary +/- and doubling results, modulo free wire sign.

    Allowing wire signs for free enlarges the implementation model, so exclusion
    in this set remains a valid lower bound for the paper's gate model. Repeated
    operands are allowed; their sum is doubling and their difference is zero.
    """
    out = set()
    for u in available:
        for v in available:
            for sign in (-1, 1):
                out.add(sign_canonical_or_zero(tuple(a+sign*b for a, b in zip(u, v))))
    return out


def at_most_two_gate_vectors(dimension):
    """Exhaust all programs with at most two gates on independent atoms."""
    atoms = tuple(tuple(int(i == j) for i in range(dimension))
                  for j in range(dimension))
    first = next_gate_results(atoms)
    reached = set(atoms) | first
    for g in first:
        reached.update(next_gate_results(atoms+(g,)))
    return reached, len(first)


def explicit_A49_cone_checks():
    # Recompute the existing complete finite certificate without writing it.
    certificate = five_reduction_certificate()
    records = certificate['all_first_forms']
    assert certificate['universe_size'] == len(records) == 156
    expected_universe = {v for v in product(range(-4, 5), repeat=4)
                         if 2 <= l1(v) <= 4 and next(x for x in v if x) > 0}
    assert {tuple(r['w']) for r in records} == expected_universe
    successful = [tuple(r['w']) for r in records if 15 in r['reachable_masks']]
    assert successful == [(1,-1,-1,-1), (1,-1,1,1),
                          (1,1,-1,1), (1,1,1,-1)]
    hist = Counter(max(mask.bit_count() for mask in r['reachable_masks']) for r in records)
    assert [hist[k] for k in range(5)] == [32, 48, 72, 0, 4]
    short4, first4 = at_most_two_gate_vectors(4)
    short5, first5 = at_most_two_gate_vectors(5)
    formal_rows = []
    pattern_counts = Counter()
    output_pair_checks = 0
    output_pre_checks = 0
    total_trees = 0
    for w in successful:
        assert all(abs(x) == 1 for x in w) and l1(w) == 4
        wc = canonical(w)[0]
        assert wc not in short4
        pre_trees = trees(wc)
        pre_gate_count = min(map(len, pre_trees))
        assert pre_gate_count == 3
        total_trees += len(pre_trees)
        vectors = []
        outputs = []
        for i, target in enumerate(A49):
            others = [j for j in range(4) if j != i]
            reps = representations(target, [w]+[A49[j] for j in others])
            assert len(reps) == 1
            coefficients, inputs = reps[0]
            assert coefficients[1:] == (0, 0, 0)
            assert abs(coefficients[0]) == 1 and l1(inputs) == 3
            assert tuple(inputs[k]+coefficients[0]*w[k] for k in range(4)) == target
            formal = tuple(inputs)+(coefficients[0],)
            vectors.append(formal)
            pattern = tuple(sorted((abs(x) for x in formal if x), reverse=True))
            assert pattern in [(3, 1), (2, 1, 1)] and l1(formal) == 4
            pattern_counts[pattern] += 1
            fc = canonical(formal)[0]
            assert fc not in short5
            family = trees(fc)
            cone_gates = min(map(len, family))
            assert cone_gates == 3
            total_trees += len(family)
            outputs.append(dict(output=i, unique_representation_count=1,
                                reset_w_coefficient=coefficients[0],
                                input_coefficients=inputs,
                                other_output_coefficients=coefficients[1:],
                                formal_vector=formal, absolute_nonzero_pattern=pattern,
                                coefficient_norm=4, generated_with_at_most_two_gates=False,
                                sign_compatible_DAGs=len(family), minimum_cone_gates=cone_gates))
        assert Counter(tuple(row['absolute_nonzero_pattern']) for row in outputs) == {
            (3, 1):1, (2, 1, 1):3}
        output_forms = [{canonical(s)[0] for s in subforms(v) if l1(s) >= 2}
                        for v in vectors]
        pre_forms = {canonical(s+(0,))[0] for s in subforms(w) if 2 <= l1(s) < 4}
        intersections = []
        for i, j in combinations(range(4), 2):
            shared = output_forms[i] & output_forms[j]
            assert not shared
            output_pair_checks += 1
            intersections.append(dict(output_pair=[i, j], nontrivial_shared_forms=0))
        for i in range(4):
            shared = output_forms[i] & pre_forms
            assert not shared
            output_pre_checks += 1
        formal_rows.append(dict(first_normalized_form=w,
                                pre_normalization_gate_minimum=pre_gate_count,
                                pre_normalization_two_gate_exclusion=True,
                                pre_normalization_sign_compatible_DAGs=len(pre_trees),
                                outputs=outputs, output_pair_intersections=intersections,
                                output_proper_preform_intersections=[0]*4,
                                total_disjoint_gate_lower_bound=pre_gate_count+12))
    assert pattern_counts == {(3, 1):4, (2, 1, 1):12}
    # Tight load alone does not give a three-gate lower bound: reuse permits
    # four atom occurrences to be generated with two gates.
    example = circuit('load_four_two_gate_counterexample', 2,
                      [('s', [('x0',1), ('x1',1)]), ('y', [('s',2)])], ['y'])
    assert example.matrix() == ((2, 2),)
    example_operations = sum(not n.output for n in example.nodes)
    example_schedule = optimize(example, 4)
    assert example_operations == 2
    assert [row['raw_load'] for row in example_schedule['details'][:2]] == [2, 4]
    assert example_schedule['selected'] == []
    assert l1(example.matrix()[0]) == 4
    short2, _ = at_most_two_gate_vectors(2)
    assert (2, 2) in short2
    return dict(candidate_forms=156, first_form_histogram=[hist[k] for k in range(5)],
                successful_forms=4, unique_formal_output_representations=16,
                absolute_coefficient_pattern_counts=[
                    dict(pattern=list(pattern), count=count)
                    for pattern, count in sorted(pattern_counts.items(), reverse=True)],
                first_gate_vector_counts={'four_atoms':first4, 'five_atoms':first5},
                at_most_two_gate_vector_counts={'four_atoms':len(short4), 'five_atoms':len(short5)},
                explicit_two_gate_exclusions=20,
                sign_compatible_DAGs_inspected=total_trees,
                output_output_disjointness_checks=output_pair_checks,
                output_proper_first_form_disjointness_checks=output_pre_checks,
                verified_pre_normalization_gates=3, verified_each_output_cone_gates=3,
                verified_total_disjoint_lower_bound=15, records=formal_rows,
                lemma_scope='After the tight-load argument reduces every cone to its unique sign-compatible formal vector; doubling and repeated wire use remain permitted',
                logical_caveat='Load equality alone is insufficient; the coefficient-pattern and no-sharing checks are essential',
                tight_load_counterexample=dict(expression='2*(x+y)', operations=2,
                    sequence=['s=x+y', 'y=2*s'], output_coefficients=[2,2],
                    coefficient_norm=4, raw_loads=[2,4],
                    reason='The same intermediate wire is used twice at the doubling gate'))


def snapshot_existing_data():
    return {str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted((ROOT/'data').rglob('*'))
            if path.is_file() and path.name != 'revision6_checks.json'}


def main():
    before = snapshot_existing_data()
    widths = word_size_checks()
    print('Word-size checks passed: BN254 H(256)=5; L5/L7/L16 widths=256/257/258.', flush=True)
    implementations = modulus_and_schedule_checks()
    print('BN254 minor and schedule checks passed:', implementations['minor_checks'],
          'minors; normalization counts', [r['normalizations'] for r in implementations['schedules']], flush=True)
    cones = explicit_A49_cone_checks()
    print('A49 focused cone checks passed: 156 first forms; four successful; sixteen unique output cones.', flush=True)
    after = snapshot_existing_data()
    assert before == after
    report = dict(version='revision6', word_sizes=widths,
                  BN254_implementations=implementations, A49_focused_lemma=cones,
                  preservation=dict(existing_data_files=len(before),
                                    every_existing_data_hash_unchanged=True))
    destination = ROOT/'data/revision6_checks.json'
    destination.write_text(json.dumps(report, indent=2)+'\n')
    print('Saved', destination.name, '; existing data files unchanged:', len(before), flush=True)


if __name__ == '__main__':
    main()
