"""UNEXECUTED analytical fixtures; parent source/math review required before use.

No main, numerical package or actual artifact loading. These toy full tensors
need not have the Hubbard particle sector and contribute zero campaign ideas.
"""
from fractions import Fraction

def _pack(entries):
    real, imag = [0.0] * 4096, [0.0] * 4096
    for index, re, im in entries:
        real[index], imag[index] = re, im
    return {'shape': [4096], 'order': 'F', 'real': real, 'imag': imag}

def _fraction(record):
    return Fraction(record['numerator'], record['denominator'])

def _ln2_oracle():
    # Independent positive series ln2=sum_{k>=1}1/(k*2^k).
    # Tail after256 terms is at most1/(257*2^256).
    lower = sum((Fraction(1, k * (1 << k)) for k in range(1, 257)), Fraction(0))
    return lower, lower + Fraction(1, 257 * (1 << 256))

def run_analytic_fixtures(core):
    cases = [
        ('product-real', _pack([(0, 1.0, 0.0)]), Fraction(1), 0),
        ('product-complex-phase', _pack([(0, 0.0, 1.0)]), Fraction(1), 0),
        ('Bell-two-unnormalized', _pack([(0, 1.0, 0.0), (65, 1.0, 0.0)]), Fraction(1, 2), 1),
        ('maximal64-unnormalized', _pack([(65 * j, 1.0, 0.0) for j in range(64)]), Fraction(1, 64), 6),
        ('nonuniform-Schmidt-3-4', _pack([(0, 0.75, 0.0), (65, 1.0, 0.0)]), Fraction(337, 625), None),
        ('Bell-local-permutation-phase', _pack([(1 + 64 * 2, -1.0, 0.0), (0 + 64 * 3, 0.0, 1.0)]), Fraction(1, 2), 1),
    ]
    ln2_lower, ln2_upper = _ln2_oracle()
    records = []
    for label, pack, purity, log2_multiple in cases:
        result = core.enclose_raw_site_state(pack)
        assert _fraction(result['normalized_purity_lower']) <= purity <= _fraction(result['normalized_purity_upper'])
        if log2_multiple == 0:
            assert _fraction(result['S2_lower']) <= 0 <= _fraction(result['S2_upper'])
        elif log2_multiple is not None:
            assert _fraction(result['S2_lower']) <= log2_multiple * ln2_lower
            assert log2_multiple * ln2_upper <= _fraction(result['S2_upper'])
        records.append({'label': label, 'analytical_purity': str(purity), 'status': 'PASS_TOY_ANALYTICAL_ONLY', 'scientific_credit': 0})
    conditional = run_conditional_join_fixtures(core)
    return {'status': 'PASS_ANALYTICAL_FIXTURES_ONLY_ZERO_CREDIT', 'records': records,
            'conditional_tensor_join': conditional,
            'actual_scientific_artifacts_used': False, 'execution_authority': False}

def run_conditional_join_fixtures(core):
    """Genuine analytic product tensor/vector, never a fake runtime callback.

    The tensor product of six (1,1,4) vacuum tensors equals the supplied complete
    vector exactly, so E=0 is an analytic certificate. These objects are marked
    synthetic and do not attest an original P13 producer or library execution.
    """
    tensors = [
        {'shape': [1, 1, 4], 'order': 'F', 'real': [1.0, 0.0, 0.0, 0.0],
         'imag': [0.0, 0.0, 0.0, 0.0]}
        for _ in range(6)
    ]
    state = {'axes': ['left', 'right', 'physical'], 'sites': tensors,
             'export_stage': 'AFTER_ACTUAL_USER_ENTROPY_CALL',
             'external_norm_scalar': 1.0,
             'complete_site_amplitudes': _pack([(0, 1.0, 0.0)])}
    certificate = {
        'backend': 'native', 'arm': 'primary', 'time_numerator': 0,
        'time_denominator': 8, 'library_S2_original': 0.0,
        'complete_tensor_export_identity_error_upper': 0.0,
        'fresh_normalized_MPS_purity': {'center': 1.0, 'radius': 0.0},
        'fresh_normalized_MPS_entropy': {'center': 0.0, 'radius': 0.0},
        'TeNPy_omitted_entropy_allowance': 0.0,
        'TeNPy_Schmidt_normalization_entropy_allowance': 0.0,
        'TeNPy_omitted_purity_allowance': 0.0,
    }
    raw = core.check_registered_export(state, 'native', 0.0)
    joined = core.join_original_tensor_identity(raw, certificate, 0.0)
    assert joined['status'] == 'PASS_CONDITIONAL_ORIGINAL_TENSOR_IDENTITY_JOIN_ONLY_PENDING_WHOLE_AUDIT'
    assert _fraction(joined['conditional_tensor_purity_lower']) <= 1 <= _fraction(joined['conditional_tensor_purity_upper'])
    assert _fraction(joined['conditional_tensor_S2_lower']) <= 0 <= _fraction(joined['conditional_tensor_S2_upper'])
    assert joined['source_bound_original_tensor_export_identity_error_recomputed'] is False
    records = [{'label': 'exact-product-MPS-E0-conditional-join', 'status': 'PASS_TOY_ANALYTICAL_ONLY'}]
    negative_cases = [
        ('negative-original-E', {'complete_tensor_export_identity_error_upper': -1.0}),
        ('negative-omitted-entropy-budget', {'TeNPy_omitted_entropy_allowance': -1.0}),
        ('negative-normalization-budget', {'TeNPy_Schmidt_normalization_entropy_allowance': -1.0}),
        ('negative-omitted-purity-budget', {'TeNPy_omitted_purity_allowance': -1.0}),
        ('wrong-backend-key', {'backend': 'tenpy'}),
        ('unregistered-arm-key', {'arm': 'other'}),
        ('out-of-range-time-key', {'time_numerator': 33}),
    ]
    for label, changes in negative_cases:
        bad = certificate.copy()
        bad.update(changes)
        try:
            core.join_original_tensor_identity(raw, bad, 0.0)
        except core.PhysicalHold:
            records.append({'label': label, 'status': 'PASS_TOY_EXPECTED_HOLD'})
        else:
            raise AssertionError('Invalid synthetic tensor certificate accepted: ' + label)
    return {'status': 'PASS_CONDITIONAL_JOIN_ANALYTICAL_AND_NEGATIVE_FIXTURES_ZERO_CREDIT',
            'records': records, 'scientific_credit': 0,
            'actual_scientific_artifacts_or_runtime_callbacks_used': False,
            'different_but_valid_time_or_arm_requires_qualified_caller_binding': True}
