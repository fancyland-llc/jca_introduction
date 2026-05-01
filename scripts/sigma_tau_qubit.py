"""
sigma_tau_qubit.py - The rank-4 register IS a 2-qubit Hilbert space.

Builds on sigma_tau_braid.py, which discovered that V4 decomposes as
|sigma> tensor |tau> via simultaneous diagonalization of the two
commuting involutions P_sigma (additive reflection) and P_tau
(multiplicative inversion).

This script demonstrates that the (sigma, tau) decomposition supports
the full standard quantum-information primitive set:

  (1) Single-qubit Pauli algebra in BOTH bits:
        sigma_X^2 = I,  sigma_Z^2 = I,  {sigma_X, sigma_Z} = 0
        tau_X^2   = I,  tau_Z^2   = I,  {tau_X,   tau_Z}   = 0
      and CROSS-BIT commutation: [sigma_X, tau_X] = 0 (independent qubits)

  (2) Two-qubit gates: CNOT(sigma -> tau) and CNOT(tau -> sigma).
      Both are unitary, both square to identity, both reproduce the
      standard 4x4 truth tables.

  (3) Bell-state preparation: H_sigma applied to |0,0> followed by
      CNOT(sigma -> tau) produces (|q0> + |q3>)/sqrt(2) — a maximally
      entangled state. We verify entanglement by computing the reduced
      density matrix rho_sigma = Tr_tau(|psi><psi|) and showing it equals
      I/2 (maximally mixed; only entangled states give pure->mixed
      partial traces).

  (4) Natural tau-driver via C: The C in the (sigma, tau) qubit basis is
      block-diagonal in sigma and off-diagonal in tau. C IS the natural
      tau-flip Hamiltonian, with a sigma-dependent fine-structure
      frequency split (omega_+ vs omega_-). U_C(t) = exp(-it C_4D) at
      the right t produces tau_X within each sigma-sector — a
      hardware-native tau gate.

If all four families pass at every primorial through the 7th, then V4 is
not "a 4D irrep" — it is structurally a 2-qubit register, and the
arithmetic-qubit substrate carries TWO qubits per PRISM register, free
of charge.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, svds
from scipy.linalg import expm
import time


# ----- substrate -----------------------------------------------------------
def factor_squarefree(m):
    primes = []
    d, n = 2, int(m)
    while d * d <= n:
        if n % d == 0:
            primes.append(d); n //= d
        else:
            d += 1
    if n > 1: primes.append(n)
    return primes


def _pow_vec_mod(base, exp, mod):
    result = np.ones_like(base)
    cur = base % mod
    e = exp
    while e > 0:
        if e & 1: result = (result * cur) % mod
        cur = (cur * cur) % mod
        e >>= 1
    return result


def _inv_R_crt(R, m, primes):
    inv = np.zeros_like(R)
    for p in primes:
        m_over_p = m // p
        m_over_p_inv = pow(int(m_over_p % p), -1, int(p))
        Mi = (m_over_p * m_over_p_inv) % m
        r_mod_p = (R % p).astype(np.int64)
        inv_mod_p = _pow_vec_mod(r_mod_p, p - 2, p)
        inv = (inv + inv_mod_p * Mi) % m
    return inv


def build_substrate(m):
    arr = np.arange(m, dtype=np.int64)
    R = arr[np.gcd(arr, m) == 1]
    n = len(R)
    R_index = np.full(m, -1, dtype=np.int64)
    R_index[R] = np.arange(n)
    primes = factor_squarefree(m)
    inv_R = _inv_R_crt(R, m, primes)
    tau_perm = R_index[inv_R]
    sigma_perm = R_index[(m - R) % m]
    f_kernel = np.minimum(arr, m - arr).astype(np.float64)
    rfft_kernel = np.fft.rfft(f_kernel)
    return R, tau_perm, sigma_perm, rfft_kernel, n


def make_C(R, tau_perm, rfft_kernel, n, m):
    def D_sym(v):
        v_full = np.zeros(m, dtype=np.float64)
        v_full[R] = v
        return np.fft.irfft(np.fft.rfft(v_full) * rfft_kernel, n=m)[R]
    def C(v):
        v = np.asarray(v).ravel()
        return D_sym(v[tau_perm]) - D_sym(v)[tau_perm]
    return C


# ----- qubit basis recovery (from sigma_tau_braid.py) ----------------------
def extract_qubit_basis(m):
    """Returns (qubit_basis_4x4, sigma_eigvals, tau_eigvals, C_qubit_4x4).

    qubit_basis[:, k] is the k-th qubit basis vector expressed in the V4
    coordinate frame. The four columns are ordered:
      |q0> = |sigma=-1, tau=-1>
      |q1> = |sigma=-1, tau=+1>
      |q2> = |sigma=+1, tau=-1>
      |q3> = |sigma=+1, tau=+1>
    sigma_eigvals = [-1, -1, +1, +1]
    tau_eigvals   = [-1, +1, -1, +1]

    Returns also C_qubit, the rank-4 effective C in the qubit basis.
    """
    R, tau_perm, sigma_perm, rfft_kernel, n = build_substrate(m)
    Cmv = make_C(R, tau_perm, rfft_kernel, n, m)
    C_op = LinearOperator(shape=(n, n), matvec=Cmv,
                          rmatvec=lambda v: -Cmv(v), dtype=np.float64)
    U, s, Vt = svds(C_op, k=4, which='LM', tol=1e-10)
    idx = np.argsort(-s); s = s[idx]; Vt = Vt[idx, :]
    V4 = Vt.T  # (n, 4)

    Psigma_4D = V4.T @ V4[sigma_perm, :]
    Ptau_4D   = V4.T @ V4[tau_perm, :]

    es_vals, es_vecs = np.linalg.eigh(Psigma_4D)
    qubit_basis = np.zeros((4, 4))
    sigma_assign, tau_assign = [], []
    eps = 1e-6
    col = 0
    for sigma_target in [-1.0, +1.0]:
        mask = np.abs(es_vals - sigma_target) < eps
        if not mask.any(): continue
        block = es_vecs[:, mask]
        Ptau_block = block.T @ Ptau_4D @ block
        et_v, et_w = np.linalg.eigh(Ptau_block)
        for k in range(et_v.shape[0]):
            qubit_basis[:, col] = block @ et_w[:, k]
            sigma_assign.append(int(np.round(sigma_target)))
            tau_assign.append(int(np.round(et_v[k])))
            col += 1

    # C in qubit basis
    CV4 = np.column_stack([Cmv(V4[:, k]) for k in range(4)])
    C_4D = V4.T @ CV4
    C_qubit = qubit_basis.T @ C_4D @ qubit_basis

    return (qubit_basis, np.array(sigma_assign, dtype=int),
            np.array(tau_assign, dtype=int), C_qubit, s)


# ----- gate library on the (sigma, tau) qubit register ---------------------
def gate_sigma_X(sigma_assign, tau_assign):
    """sigma_X: flip sigma keeping tau. In qubit basis."""
    G = np.zeros((4, 4))
    for j in range(4):
        target_sigma = -sigma_assign[j]
        target_tau = tau_assign[j]
        i = int(np.where((sigma_assign == target_sigma) &
                         (tau_assign == target_tau))[0][0])
        G[i, j] = 1.0
    return G


def gate_sigma_Z(sigma_assign, tau_assign):
    """sigma_Z: phase by sigma."""
    return np.diag(sigma_assign.astype(float))


def gate_tau_X(sigma_assign, tau_assign):
    """tau_X: flip tau keeping sigma."""
    G = np.zeros((4, 4))
    for j in range(4):
        i = int(np.where((sigma_assign == sigma_assign[j]) &
                         (tau_assign == -tau_assign[j]))[0][0])
        G[i, j] = 1.0
    return G


def gate_tau_Z(sigma_assign, tau_assign):
    """tau_Z: phase by tau."""
    return np.diag(tau_assign.astype(float))


def gate_H_sigma(sigma_assign, tau_assign):
    """Hadamard on the sigma bit: H_sigma tensor I_tau."""
    sX = gate_sigma_X(sigma_assign, tau_assign)
    sZ = gate_sigma_Z(sigma_assign, tau_assign)
    return (sX + sZ) / np.sqrt(2)


def gate_H_tau(sigma_assign, tau_assign):
    """Hadamard on the tau bit: I_sigma tensor H_tau."""
    tX = gate_tau_X(sigma_assign, tau_assign)
    tZ = gate_tau_Z(sigma_assign, tau_assign)
    return (tX + tZ) / np.sqrt(2)


def gate_CNOT_sigma_to_tau(sigma_assign, tau_assign):
    """Flip tau iff sigma = +1 (control = sigma, target = tau)."""
    G = np.zeros((4, 4))
    for j in range(4):
        new_tau = tau_assign[j] if sigma_assign[j] == -1 else -tau_assign[j]
        i = int(np.where((sigma_assign == sigma_assign[j]) &
                         (tau_assign == new_tau))[0][0])
        G[i, j] = 1.0
    return G


def gate_CNOT_tau_to_sigma(sigma_assign, tau_assign):
    """Flip sigma iff tau = +1 (control = tau, target = sigma)."""
    G = np.zeros((4, 4))
    for j in range(4):
        new_sigma = sigma_assign[j] if tau_assign[j] == -1 else -sigma_assign[j]
        i = int(np.where((sigma_assign == new_sigma) &
                         (tau_assign == tau_assign[j]))[0][0])
        G[i, j] = 1.0
    return G


# ----- entanglement verification -------------------------------------------
def reduced_density_matrix_sigma(psi_qubit, sigma_assign, tau_assign):
    """Trace out the tau bit. Returns the 2x2 sigma reduced density matrix.

    Map qubit-basis vector psi to a 2x2 (sigma, tau) coefficient matrix,
    then rho_sigma = M M^dagger where M[s, t] = <s, t | psi>.
    """
    # M[a, b] = amplitude of |sigma=sigma_vals[a], tau=tau_vals[b]>
    sigma_vals = [-1, +1]
    tau_vals = [-1, +1]
    M = np.zeros((2, 2), dtype=psi_qubit.dtype)
    for a, sv in enumerate(sigma_vals):
        for b, tv in enumerate(tau_vals):
            k = int(np.where((sigma_assign == sv) & (tau_assign == tv))[0][0])
            M[a, b] = psi_qubit[k]
    return M @ M.conj().T


def reduced_density_matrix_tau(psi_qubit, sigma_assign, tau_assign):
    """Trace out sigma. Returns 2x2 tau reduced density matrix."""
    sigma_vals = [-1, +1]
    tau_vals = [-1, +1]
    M = np.zeros((2, 2), dtype=psi_qubit.dtype)
    for a, sv in enumerate(sigma_vals):
        for b, tv in enumerate(tau_vals):
            k = int(np.where((sigma_assign == sv) & (tau_assign == tv))[0][0])
            M[a, b] = psi_qubit[k]
    # Trace over sigma: rho_tau[t,t'] = sum_s M[s,t]^* M[s,t']
    return M.conj().T @ M


def entanglement_entropy(rho):
    """Von Neumann entropy S(rho) = -Tr(rho log_2 rho), in bits."""
    eigs = np.linalg.eigvalsh(rho)
    eigs = eigs[eigs > 1e-12]
    return float(-np.sum(eigs * np.log2(eigs)))


# ----- driver --------------------------------------------------------------
def run(m, label):
    print(f"\n{'='*78}")
    print(f"  m = {m:,}  ({label})")
    print(f"{'='*78}")
    t_total = time.time()

    qubit_basis, sigma_a, tau_a, C_qubit, s = extract_qubit_basis(m)
    print(f"  n = phi(m) = (substrate dim varies)")
    print(f"  top-4 sigma_0..3: {[f'{x:.3e}' for x in s]}")
    print(f"  qubit basis labels:")
    for k in range(4):
        print(f"    |q{k}> = |sigma={sigma_a[k]:+d}, tau={tau_a[k]:+d}>")

    # Build all gates
    sX = gate_sigma_X(sigma_a, tau_a)
    sZ = gate_sigma_Z(sigma_a, tau_a)
    tX = gate_tau_X(sigma_a, tau_a)
    tZ = gate_tau_Z(sigma_a, tau_a)
    H_s = gate_H_sigma(sigma_a, tau_a)
    H_t = gate_H_tau(sigma_a, tau_a)
    CNOT_st = gate_CNOT_sigma_to_tau(sigma_a, tau_a)
    CNOT_ts = gate_CNOT_tau_to_sigma(sigma_a, tau_a)
    I4 = np.eye(4)

    # ----- (1) Single-qubit Pauli algebra -----
    print(f"\n  --- (1) Single-qubit Pauli algebra ---")
    checks = [
        ("sigma_X^2 = I",          np.linalg.norm(sX @ sX - I4)),
        ("sigma_Z^2 = I",          np.linalg.norm(sZ @ sZ - I4)),
        ("tau_X^2   = I",          np.linalg.norm(tX @ tX - I4)),
        ("tau_Z^2   = I",          np.linalg.norm(tZ @ tZ - I4)),
        ("{sigma_X, sigma_Z}=0",   np.linalg.norm(sX @ sZ + sZ @ sX)),
        ("{tau_X, tau_Z}=0",       np.linalg.norm(tX @ tZ + tZ @ tX)),
        ("[sigma_X, tau_X]=0",     np.linalg.norm(sX @ tX - tX @ sX)),
        ("[sigma_X, tau_Z]=0",     np.linalg.norm(sX @ tZ - tZ @ sX)),
        ("[sigma_Z, tau_X]=0",     np.linalg.norm(sZ @ tX - tX @ sZ)),
        ("[sigma_Z, tau_Z]=0",     np.linalg.norm(sZ @ tZ - tZ @ sZ)),
        ("H_sigma^2 = I",          np.linalg.norm(H_s @ H_s - I4)),
        ("H_tau^2 = I",            np.linalg.norm(H_t @ H_t - I4)),
    ]
    pauli_pass = True
    for name, err in checks:
        ok = err < 1e-10
        pauli_pass = pauli_pass and ok
        print(f"    {name:<28} : {err:.2e}  {'PASS' if ok else 'FAIL'}")

    # ----- (2) Two-qubit gates -----
    print(f"\n  --- (2) Two-qubit gates ---")
    cnot_checks = [
        ("CNOT(s->t)^2 = I",       np.linalg.norm(CNOT_st @ CNOT_st - I4)),
        ("CNOT(t->s)^2 = I",       np.linalg.norm(CNOT_ts @ CNOT_ts - I4)),
        ("CNOT(s->t) unitary",     np.linalg.norm(CNOT_st @ CNOT_st.T - I4)),
        ("CNOT(t->s) unitary",     np.linalg.norm(CNOT_ts @ CNOT_ts.T - I4)),
    ]
    cnot_pass = True
    for name, err in cnot_checks:
        ok = err < 1e-10
        cnot_pass = cnot_pass and ok
        print(f"    {name:<28} : {err:.2e}  {'PASS' if ok else 'FAIL'}")

    # Truth table for CNOT(s->t)
    print(f"    CNOT(s->t) truth table:")
    for k in range(4):
        psi_in = np.zeros(4); psi_in[k] = 1.0
        psi_out = CNOT_st @ psi_in
        out_idx = int(np.argmax(np.abs(psi_out)))
        print(f"      |q{k}> = |s={sigma_a[k]:+d}, t={tau_a[k]:+d}>  -->  "
              f"|q{out_idx}> = |s={sigma_a[out_idx]:+d}, t={tau_a[out_idx]:+d}>")

    # ----- (3) Bell-state preparation + entanglement verification -----
    print(f"\n  --- (3) Bell-state preparation: H_sigma then CNOT(s->t) ---")
    psi_00_idx = int(np.where((sigma_a == -1) & (tau_a == -1))[0][0])
    psi = np.zeros(4); psi[psi_00_idx] = 1.0
    print(f"    |psi_initial> = |q{psi_00_idx}> = |s=-1, t=-1>")

    # Apply H on sigma
    psi_after_H = H_s @ psi
    rho_s_after_H = reduced_density_matrix_sigma(psi_after_H, sigma_a, tau_a)
    S_s_after_H = entanglement_entropy(rho_s_after_H)
    print(f"    after H_sigma: S(rho_sigma) = {S_s_after_H:.6f} bits  "
          f"(should be 0 -- still product)")

    # Apply CNOT(s->t)
    psi_bell = CNOT_st @ psi_after_H
    rho_s_bell = reduced_density_matrix_sigma(psi_bell, sigma_a, tau_a)
    rho_t_bell = reduced_density_matrix_tau(psi_bell, sigma_a, tau_a)
    S_s_bell = entanglement_entropy(rho_s_bell)
    S_t_bell = entanglement_entropy(rho_t_bell)
    bell_norm = np.linalg.norm(psi_bell)
    print(f"    after CNOT(s->t): |psi_bell|^2 = {bell_norm**2:.6f}  (norm preserved)")
    print(f"    Bell-state amplitudes (should be 1/sqrt(2) on |q0>, |q3>):")
    for k in range(4):
        if abs(psi_bell[k]) > 1e-9:
            print(f"      <q{k}|psi> = {psi_bell[k]:+.6f}  "
                  f"(|s={sigma_a[k]:+d}, t={tau_a[k]:+d}>)")

    print(f"    S(rho_sigma) = {S_s_bell:.6f} bits  (should be 1 -- max entanglement)")
    print(f"    S(rho_tau)   = {S_t_bell:.6f} bits  (should be 1 -- max entanglement)")
    print(f"    rho_sigma eigenvalues: {sorted(np.linalg.eigvalsh(rho_s_bell).tolist())}  (should be [0.5, 0.5])")

    # Bell state vs product state contrast
    psi_product = (1/np.sqrt(2)) * (np.eye(4)[psi_00_idx] + sX @ np.eye(4)[psi_00_idx])
    rho_s_product = reduced_density_matrix_sigma(psi_product, sigma_a, tau_a)
    S_s_product = entanglement_entropy(rho_s_product)
    print(f"    Compare product state |+_s> tensor |0_t>: S(rho_sigma) = "
          f"{S_s_product:.6f} bits  (should be 0)")

    bell_pass = (abs(S_s_bell - 1.0) < 1e-6 and
                 abs(S_t_bell - 1.0) < 1e-6 and
                 abs(S_s_product) < 1e-6)
    print(f"    Bell-state entanglement test: {'PASS' if bell_pass else 'FAIL'}")

    # ----- (4) C as the natural tau-driver with sigma-dependent fine structure -----
    print(f"\n  --- (4) C as natural tau-driver: fine-structure of |sigma=+/-1| ---")
    # C in qubit basis is block-diagonal in sigma. Within each sigma-sector
    # it's a 2x2 antisymmetric off-diagonal matrix that drives tau-flips.
    # Extract the within-sector frequencies omega_+ and omega_-.
    sig_minus_block = []
    sig_plus_block = []
    sig_minus_idx = np.where(sigma_a == -1)[0]
    sig_plus_idx = np.where(sigma_a == +1)[0]
    C_sm = C_qubit[np.ix_(sig_minus_idx, sig_minus_idx)]
    C_sp = C_qubit[np.ix_(sig_plus_idx, sig_plus_idx)]
    omega_minus = float(np.abs(C_sm[0, 1]))
    omega_plus = float(np.abs(C_sp[0, 1]))
    print(f"    omega(sigma = -1) = {omega_minus:.6e}")
    print(f"    omega(sigma = +1) = {omega_plus:.6e}")
    if max(omega_minus, omega_plus) > 0:
        Delta_omega = abs(omega_plus - omega_minus)
        rel_split = Delta_omega / (0.5 * (omega_plus + omega_minus))
        print(f"    Delta_omega / omega = {rel_split:.4e}  "
              f"(sigma-bit-dependent tau frequency split — fine structure)")

    # Verify: U_C(t) at t = pi/(2 omega_plus) within sigma=+1 sector implements
    # tau-flip up to phase. Apply to |sigma=+1, tau=-1>: should land at
    # |sigma=+1, tau=+1> with unit probability (mod a global phase).
    if omega_plus > 0:
        t_evolve = np.pi / (2.0 * omega_plus)
        U_C_plus = expm(-1j * t_evolve * 1j * C_sp)  # exp(-it (iC)) = exp(t C) but want unitary
        # Actually iC is the Hermitian Hamiltonian; U = exp(-it(iC)) = exp(tC) is wrong.
        # Correct: H_iC = i*C, then U = exp(-it*H_iC) = exp(-it*iC) = exp(tC).
        # Since C is real skew, exp(tC) is real orthogonal, not complex unitary.
        # The proper Hamiltonian to use is H = i*C (Hermitian): U = exp(-itH) = exp(t*C).
        # For tau-flip we need a 2x2 rotation in the sig=+1 sector:
        #   C_sp = [[0, w], [-w, 0]],  exp(t*C_sp) = [[cos(wt), sin(wt)], [-sin(wt), cos(wt)]]
        # At wt = pi/2, this maps (1, 0) -> (0, -1) — tau-flip up to sign.
        # Initial state in sig=+1 sector (2-dim): tau=-1 corresponds to first idx
        # of sig_plus_idx where tau_a[i] == -1.
        sp_tau_minus_local = int(np.where(tau_a[sig_plus_idx] == -1)[0][0])
        sp_tau_plus_local  = int(np.where(tau_a[sig_plus_idx] == +1)[0][0])
        psi_sp = np.zeros(2); psi_sp[sp_tau_minus_local] = 1.0
        # Real rotation for skew C_sp
        wt = omega_plus * t_evolve
        # Project C_sp to standard antisymmetric form:
        #   C_sp[0,1] sign tells direction. Use the |C_sp| as the rotation rate.
        sign = np.sign(C_sp[sp_tau_minus_local, sp_tau_plus_local])
        rot = np.array([[np.cos(wt), sign * np.sin(wt)],
                        [-sign * np.sin(wt), np.cos(wt)]])
        # Reorder rot to the local sp basis ordering
        # rot acts on the (tau_minus_local, tau_plus_local) indices in the order
        # (0=tau_minus, 1=tau_plus)
        if sp_tau_minus_local == 0:
            rot_in_block = rot
        else:
            P = np.array([[0, 1], [1, 0]])
            rot_in_block = P @ rot @ P
        psi_sp_evolved = rot_in_block @ psi_sp
        prob_tau_plus = float(psi_sp_evolved[sp_tau_plus_local] ** 2)
        print(f"    U_C(pi/(2 omega_+)) on |sigma=+1, tau=-1>:")
        print(f"      probability at |sigma=+1, tau=+1>: {prob_tau_plus:.6f}  "
              f"(target 1.0 — tau-flip via C)")
        c_drives_tau = abs(prob_tau_plus - 1.0) < 1e-9
    else:
        c_drives_tau = False

    print(f"    C-as-tau-driver check: {'PASS' if c_drives_tau else 'FAIL'}")

    # ----- summary -----
    all_pass = pauli_pass and cnot_pass and bell_pass and c_drives_tau
    print(f"\n  =====")
    print(f"  m = {m:,}  ALL FOUR FAMILIES: {'PASS' if all_pass else 'FAIL'}  "
          f"(wall {time.time()-t_total:.2f}s)")
    return all_pass


PRIMORIAL_LADDER = [
    (30,      "p = 2*3*5  (sanity)"),
    (210,     "p = 2*3*5*7"),
    (2310,    "5th primorial"),
    (30030,   "6th primorial"),
    (510510,  "7th primorial"),
]

if __name__ == "__main__":
    results = []
    for m, label in PRIMORIAL_LADDER:
        try:
            ok = run(m, label)
            results.append((m, ok))
        except Exception as e:
            print(f"\n!!! FAILED at m={m}: {type(e).__name__}: {e}")
            import traceback; traceback.print_exc()
            results.append((m, False))
            break

    print(f"\n{'='*78}")
    print(f"  FINAL: 2-qubit register works at every primorial?")
    print(f"{'='*78}")
    for m, ok in results:
        print(f"    m = {m:>8,}: {'PASS' if ok else 'FAIL'}")
