"""SOURCE ONLY: integer full-state purity/log enclosure; no main or admission.

A separately source/runtime/resource-qualified caller must invoke this only after
its own qualification. No model, NumPy, external solver or libm log/sqrt is used.
Raw complete_site_amplitudes are already physical-site tensors, not spin-grouped
sector vectors. No observed phase fitting or charge projection is performed.
"""

STATE_LENGTH = 4096
SIDE = 64
QUANTIZATION_BITS = 96
LOG_BITS = 128
LOG_TERMS = 64
MAX_COMPONENT_ABS = 2
MAX_LOG_EXPONENT_ABS = 16

class PhysicalHold(ValueError):
    pass

def _need(ok, message):
    if not ok:
        raise PhysicalHold(message)

def _ceil_div(n, d):
    _need(type(n) is int and type(d) is int and d > 0, 'Integer positive denominator required')
    return -((-n) // d)

def _binary_ratio(value, absolute_cap):
    """Exact represented binary64 ratio, or a bounded exact integer."""
    _need(type(value) in (int, float), 'Exact int or binary64 input required; no bool')
    try:
        n, d = value.as_integer_ratio() if type(value) is float else (value, 1)
    except (ValueError, OverflowError) as error:
        raise PhysicalHold('Nonfinite binary64 input') from error
    _need(type(n) is int and type(d) is int and 0 < d and d & (d - 1) == 0,
          'Exact dyadic ratio required')
    _need(n.bit_length() <= 1076 and d.bit_length() <= 1075,
          'Fixed binary64 ratio bit cap exceeded')
    _need(abs(n) <= absolute_cap * d, 'Input magnitude exceeds fixed semantic work cap')
    return n, d

def _quantize(value):
    n, d = _binary_ratio(value, MAX_COMPONENT_ABS)
    q = (abs(n) << QUANTIZATION_BITS) // d
    return -q if n < 0 else q

def _integer_sqrt(n):
    """Exact checked integer Newton square root; no floating sqrt assumption."""
    _need(type(n) is int and 0 <= n and n.bit_length() <= 224,
          'Fixed integer norm cap exceeded')
    if n == 0:
        return 0
    x = 1 << ((n.bit_length() + 1) // 2)
    for _ in range(256):
        y = (x + n // x) // 2
        if y >= x:
            _need(x * x <= n < (x + 1) * (x + 1), 'Integer square root identity failed')
            return x
        x = y
    raise PhysicalHold('Predeclared integer square root iteration cap')

def _atanh_fixed(z_lower, z_upper):
    """Enclose 2*atanh(z), z in [0,1/2], by directed integer intervals.

    Unit is 2^-128. After64 positive odd-power terms, the remaining tail is
    <= 2*(1/2)^129/[129*(1-1/4)] < 2^-128. One upper grid unit covers it.
    Every multiplication/division below is explicitly floored/ceiled.
    """
    scale = 1 << LOG_BITS
    _need(type(z_lower) is int and type(z_upper) is int and
          0 <= z_lower <= z_upper <= scale // 2, 'Fixed nonnegative log-series domain required')
    if z_upper == 0:
        return 0, 0
    square_lower = (z_lower * z_lower) // scale
    square_upper = _ceil_div(z_upper * z_upper, scale)
    power_lower, power_upper = z_lower, z_upper
    total_lower = total_upper = 0
    for j in range(LOG_TERMS):
        odd = 2 * j + 1
        total_lower += power_lower // odd
        total_upper += _ceil_div(power_upper, odd)
        power_lower = (power_lower * square_lower) // scale
        power_upper = _ceil_div(power_upper * square_upper, scale)
    return 2 * total_lower, 2 * total_upper + 1

def _positive_log_ratio(n, d):
    """Return dyadic lower/upper log(n/d), without evaluating a float log."""
    _need(type(n) is int and type(d) is int and n > 0 and d > 0 and
          n.bit_length() <= 768 and d.bit_length() <= 768,
          'Positive bounded rational log input required')
    exponent = n.bit_length() - d.bit_length()
    at_least_power = n >= (d << exponent) if exponent >= 0 else (n << (-exponent)) >= d
    if not at_least_power:
        exponent -= 1
    _need(-MAX_LOG_EXPONENT_ABS <= exponent <= MAX_LOG_EXPONENT_ABS,
          'Registered purity logarithm exponent cap exceeded')
    m_num, m_den = (n, d << exponent) if exponent >= 0 else (n << (-exponent), d)
    _need(m_den <= m_num < 2 * m_den, 'Exact logarithm power-two reduction failed')
    zn, zd = m_num - m_den, m_num + m_den
    scale = 1 << LOG_BITS
    ln_m_lower, ln_m_upper = _atanh_fixed((zn << LOG_BITS) // zd, _ceil_div(zn << LOG_BITS, zd))
    ln2_lower, ln2_upper = _atanh_fixed(scale // 3, _ceil_div(scale, 3))
    if exponent >= 0:
        return ln_m_lower + exponent * ln2_lower, ln_m_upper + exponent * ln2_upper
    return ln_m_lower + exponent * ln2_upper, ln_m_upper + exponent * ln2_lower

def _record(n, d):
    _need(type(n) is int and type(d) is int and d > 0 and
          n.bit_length() <= 1536 and d.bit_length() <= 1536,
          'Bounded exact rational output required')
    return {'numerator': n, 'denominator': d}

def enclose_raw_site_state(pack):
    """Independently enclose normalized entropy of this complete raw4096 export.

    Does NOT reconstruct the MPS, certify time evolution, or reuse saved purity,
    entropy, full-state distance or propagated entropy radii. Returned rational
    endpoints are exact, with no inward floating conversion.
    """
    _need(type(pack) is dict and set(pack) == {'shape', 'order', 'real', 'imag'} and
          pack['shape'] == [STATE_LENGTH] and pack['order'] in ('F', 'C') and
          type(pack['real']) is list and type(pack['imag']) is list and
          len(pack['real']) == len(pack['imag']) == STATE_LENGTH,
          'Complete fixed physical-site4096 complex vector pack required')
    # A one-dimensional pack has the same F/C order. Frozen producers formed it
    # from site tensors with first site fastest, so index = row +64*column.
    qr = [_quantize(x) for x in pack['real']]
    qi = [_quantize(x) for x in pack['imag']]
    rows_r = [[qr[r + SIDE * c] for c in range(SIDE)] for r in range(SIDE)]
    rows_i = [[qi[r + SIDE * c] for c in range(SIDE)] for r in range(SIDE)]
    norm_integer = purity_integer = 0
    for a in range(SIDE):
        ar, ai = rows_r[a], rows_i[a]
        diagonal = sum(ar[k] * ar[k] + ai[k] * ai[k] for k in range(SIDE))
        norm_integer += diagonal
        purity_integer += diagonal * diagonal
        for b in range(a + 1, SIDE):
            br, bi = rows_r[b], rows_i[b]
            re = im = 0
            for k in range(SIDE):
                re += ar[k] * br[k] + ai[k] * bi[k]
                im += ai[k] * br[k] - ar[k] * bi[k]
            purity_integer += 2 * (re * re + im * im)
    _need(norm_integer > 0 and norm_integer.bit_length() <= 224,
          'Quantized state is zero or exceeds fixed norm work cap')
    norm_squared = norm_integer * norm_integer
    _need(0 < purity_integer <= norm_squared and SIDE * purity_integer >= norm_squared,
          'Exact normalized Gram purity physical identities failed')
    s = _integer_sqrt(norm_integer)
    _need(s > 128, 'Raw quantization cannot certify a positive normalizing norm')
    scale_squared = 1 << (2 * QUANTIZATION_BITS)
    norm_lower_num, norm_upper_num = (s - 128) ** 2, (s + 129) ** 2
    # Per-real-component truncation is <2^-96. Over8192 components, raw norm
    # error <=128*2^-96. Normalized ray distance <=256/s, purity error<=1024/s.
    purity_den = norm_squared * s
    purity_lower_num = purity_integer * s - 1024 * norm_squared
    purity_upper_num = purity_integer * s + 1024 * norm_squared
    _need(purity_lower_num > 0, 'Independent purity lower bound nonpositive; no clipping')
    log_lower, _ = _positive_log_ratio(purity_lower_num, purity_den)
    _, log_upper = _positive_log_ratio(purity_upper_num, purity_den)
    log_scale = 1 << LOG_BITS
    return {
        'method': 'EXACT_INTEGER_QUANTIZED_FULL4096_GRAM_PURITY_AND_DIRECTED_FIXED_POINT_LOG_V1',
        'scientific_credit': 0,
        'normalized_raw_export_only_not_MPS_or_evolution_certificate': True,
        'quantization_bits': QUANTIZATION_BITS,
        'log_bits': LOG_BITS,
        'log_series_terms': LOG_TERMS,
        'raw_norm_lower': _record(s - 128, 1 << QUANTIZATION_BITS),
        'raw_norm_upper': _record(s + 129, 1 << QUANTIZATION_BITS),
        'norm_squared_lower': _record(norm_lower_num, scale_squared),
        'norm_squared_upper': _record(norm_upper_num, scale_squared),
        'normalized_purity_lower': _record(purity_lower_num, purity_den),
        'normalized_purity_upper': _record(purity_upper_num, purity_den),
        'S2_lower': _record(-log_upper, log_scale),
        'S2_upper': _record(-log_lower, log_scale),
        'normalized_ray_quantization_distance_upper': _record(256, s),
        'purity_quantization_error_upper': _record(1024, s),
        'no_libm_log_sqrt_or_floating_Gram_arithmetic': True,
        'no_charge_projection_or_observed_phase_alignment': True,
    }

def check_registered_export(state, backend, library_S2):
    """Additional raw-export norm/library check using unchanged frozen limits.

    The separate original whole-artifact audit still requires all11692 frozen
    gates, tensor-export identity, source/runtime/queue/resource proofs. This
    independent calculation does not replace any original acceptance gate.
    """
    _need(backend in ('native', 'tenpy') and type(state) is dict and
          state.get('axes') == ['left', 'right', 'physical'] and
          type(state.get('sites')) is list and len(state['sites']) == 6 and
          state.get('export_stage') == ('AFTER_ACTUAL_USER_ENTROPY_CALL' if backend == 'native' else 'AFTER_ACTUAL_TENPY_ENTROPY_CALL'),
          'Exact original post-entropy physical-site export stage required')
    scalar_n, scalar_d = _binary_ratio(state['external_norm_scalar'], 2)
    _need(scalar_n == scalar_d, 'Frozen producer external norm scalar must already be1')
    result = enclose_raw_site_state(state['complete_site_amplitudes'])
    nl, nu = result['norm_squared_lower'], result['norm_squared_upper']
    _need(nl['denominator'] == nu['denominator'], 'Norm denominator identity failed')
    norm_num = max(abs(nl['numerator'] - nl['denominator']), abs(nu['numerator'] - nu['denominator']))
    norm_den = nl['denominator']
    limit_n, limit_d = _binary_ratio(1e-12, 1)
    _need(norm_num * limit_d <= limit_n * norm_den,
          'Additional conservative raw-export norm envelope exceeds unchanged1e-12 gate; HOLD')
    ln, ld = _binary_ratio(library_S2, 16)
    lo, hi = result['S2_lower'], result['S2_upper']
    _need(lo['denominator'] == hi['denominator'], 'Entropy denominator identity failed')
    entropy_num = max(abs(ln * lo['denominator'] - lo['numerator'] * ld),
                      abs(ln * hi['denominator'] - hi['numerator'] * ld))
    entropy_den = ld * lo['denominator']
    limit_n, limit_d = _binary_ratio(1e-10, 1)
    _need(entropy_num * limit_d <= limit_n * entropy_den,
          'Independent raw-export entropy/library discrepancy exceeds unchanged1e-10 gate; HOLD')
    result.update(
        status='PASS_ADDITIONAL_RAW_EXPORT_PHYSICAL_CHECK_ONLY_PENDING_WHOLE_AUDIT',
        backend=backend,
        raw_export_norm_squared_discrepancy_upper=_record(norm_num, norm_den),
        raw_export_S2_vs_library_discrepancy_upper=_record(entropy_num, entropy_den),
        original_MPS_tensor_reconstruction_and_evolution_bounds_recomputed=False,
        source_bound_original_tensor_identity_and_whole_audit_still_required=True,
    )
    return result

def _checked_record(value, cap=768):
    _need(type(value) is dict and set(value) == {'numerator', 'denominator'},
          'Exact internal rational endpoint record required')
    n, d = value['numerator'], value['denominator']
    _need(type(n) is int and type(d) is int and d > 0 and
          n.bit_length() <= cap and d.bit_length() <= cap,
          'Internal rational endpoint bit cap exceeded')
    return n, d

def _original_ball_interval(value):
    _need(type(value) is dict and set(value) == {'center', 'radius'},
          'Exact original serialized scalar ball required')
    cn, cd = _binary_ratio(value['center'], 16)
    rn, rd = _binary_ratio(value['radius'], 1)
    _need(rn >= 0, 'Negative original scalar ball radius')
    return (cn * rd - rn * cd, cd * rd), (cn * rd + rn * cd, cd * rd)

def _overlaps_ratio_interval(lower, upper, old_ball):
    old_lower, old_upper = _original_ball_interval(old_ball)
    ln, ld = lower
    un, ud = upper
    return ln * old_upper[1] <= old_upper[0] * ld and old_lower[0] * ud <= un * old_lower[1]

def _nonnegative_fixed_upper(value):
    n, d = _binary_ratio(value, 1)
    _need(n >= 0, 'Negative source-bound error/omission/normalization budget')
    return _ceil_div(n << LOG_BITS, d), (n, d)

def join_original_tensor_identity(raw_result, original_state_certificate, library_S2):
    """Conditional original tensor/export certificate join, not a new contraction.

    The qualified caller must bind this exact raw sample to this exact original
    per-state certificate. This method independently propagates its retained
    export-identity bound; it does NOT recompute that bound or the MPS tensors.
    Original ball overlap is a consistency check, not a proof of its radius.
    """
    _need(type(raw_result) is dict and raw_result.get('status') ==
          'PASS_ADDITIONAL_RAW_EXPORT_PHYSICAL_CHECK_ONLY_PENDING_WHOLE_AUDIT' and
          raw_result.get('method') == 'EXACT_INTEGER_QUANTIZED_FULL4096_GRAM_PURITY_AND_DIRECTED_FIXED_POINT_LOG_V1' and
          raw_result.get('scientific_credit') == 0,
          'Complete independently checked raw-export result required before conditional join')
    cert = original_state_certificate
    _need(type(cert) is dict and cert.get('backend') == raw_result.get('backend') and
          cert.get('backend') in ('native', 'tenpy') and cert.get('arm') in ('coarse', 'primary') and
          type(cert.get('time_numerator')) is int and 0 <= cert['time_numerator'] <= 32 and
          type(cert.get('time_denominator')) is int and cert['time_denominator'] == 8,
          'Original exact full-state certificate identity required')
    ln, ld = _binary_ratio(library_S2, 16)
    cn, cd = _binary_ratio(cert['library_S2_original'], 16)
    _need(ln * cd == cn * ld, 'Original library/certificate scalar association differs')
    a_num, a_den = _checked_record(raw_result['raw_norm_lower'])
    _need(a_num > 0 and a_den == 1 << QUANTIZATION_BITS,
          'Positive independent raw-export norm lower bound required')
    p_lower_num, p_lower_den = _checked_record(raw_result['normalized_purity_lower'])
    p_upper_num, p_upper_den = _checked_record(raw_result['normalized_purity_upper'])
    _need(p_lower_num > 0 and p_lower_den == p_upper_den and p_lower_num <= p_upper_num,
          'Independent positive raw purity interval required')
    error_upper, (error_n, error_d) = _nonnegative_fixed_upper(cert['complete_tensor_export_identity_error_upper'])
    limit_n, limit_d = _binary_ratio(1e-10, 1)
    _need(error_n * limit_d <= limit_n * error_d,
          'Original tensor/export identity bound exceeds unchanged1e-10 gate')
    log_scale = 1 << LOG_BITS
    _need(error_upper * a_den < log_scale * a_num,
          'Original tensor norm not certified positive through raw-export identity')
    # E is rounded UP only to the fixed2^-128 grid. This controls denominator
    # growth, including binary64 subnormal E, and never shrinks an original error.
    allowance_num, allowance_den = 8 * error_upper * a_den, log_scale * a_num
    joined_den = p_lower_den * allowance_den
    joined_lower_num = p_lower_num * allowance_den - allowance_num * p_lower_den
    joined_upper_num = p_upper_num * allowance_den + allowance_num * p_upper_den
    _need(joined_lower_num > 0, 'Conditional tensor purity lower bound nonpositive; no clipping')
    log_lower, _ = _positive_log_ratio(joined_lower_num, joined_den)
    _, log_upper = _positive_log_ratio(joined_upper_num, joined_den)
    entropy_lower, entropy_upper = -log_upper, -log_lower
    _need(_overlaps_ratio_interval((joined_lower_num, joined_den), (joined_upper_num, joined_den),
                                  cert['fresh_normalized_MPS_purity']),
          'Conditional tensor/raw purity interval inconsistent with original fresh tensor-purity ball')
    _need(_overlaps_ratio_interval((entropy_lower, log_scale), (entropy_upper, log_scale),
                                  cert['fresh_normalized_MPS_entropy']),
          'Conditional tensor/raw entropy interval inconsistent with original fresh tensor-entropy ball')
    omission, _ = _nonnegative_fixed_upper(cert['TeNPy_omitted_entropy_allowance'])
    normalization, _ = _nonnegative_fixed_upper(cert['TeNPy_Schmidt_normalization_entropy_allowance'])
    _nonnegative_fixed_upper(cert['TeNPy_omitted_purity_allowance'])
    if cert['backend'] == 'native':
        _need(omission == normalization == 0 and cert['TeNPy_omitted_purity_allowance'] == 0,
              'Native certificate acquired a TeNPy omission/normalization budget')
    discrepancy_num = max(abs(ln * log_scale - entropy_lower * ld),
                          abs(ln * log_scale - entropy_upper * ld))
    discrepancy_num += (omission + normalization) * ld
    discrepancy_den = ld * log_scale
    _need(discrepancy_num * limit_d <= limit_n * discrepancy_den,
          'Additional conditional tensor/library discrepancy plus source budgets exceeds unchanged1e-10 gate; HOLD')
    return {
        'status': 'PASS_CONDITIONAL_ORIGINAL_TENSOR_IDENTITY_JOIN_ONLY_PENDING_WHOLE_AUDIT',
        'scientific_credit': 0,
        'state_key_from_original_certificate': [cert['backend'], cert['arm'], cert['time_numerator']],
        'raw_sample_to_original_certificate_binding_is_qualified_caller_obligation': True,
        'source_bound_original_tensor_export_identity_error_recomputed': False,
        'original_MPS_tensors_or_eigensystem_evolution_recomputed': False,
        'source_error_bound_rounded_up_grid_bits': LOG_BITS,
        'source_error_bound_upward_grid': _record(error_upper, log_scale),
        'normalized_tensor_export_purity_allowance_upper': _record(allowance_num, allowance_den),
        'conditional_tensor_purity_lower': _record(joined_lower_num, joined_den),
        'conditional_tensor_purity_upper': _record(joined_upper_num, joined_den),
        'conditional_tensor_S2_lower': _record(entropy_lower, log_scale),
        'conditional_tensor_S2_upper': _record(entropy_upper, log_scale),
        'original_fresh_purity_and_entropy_ball_overlap_consistency': True,
        'original_ball_radii_independently_recomputed': False,
        'source_omission_and_normalization_budgets_rounded_up': True,
        'conditional_tensor_library_discrepancy_plus_source_budgets_upper': _record(discrepancy_num, discrepancy_den),
        'all_original_frozen_gates_and_source_runtime_whole_audit_still_required': True,
    }
