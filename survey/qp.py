"""Active-set solver for bounded, equality-constrained convex quadratic programs.

Minimize  0.5 z' D z - q' z  subject to  A z = b, lo <= z <= hi,
where D is a strictly positive diagonal.  Equality rows may be rank
deficient: each reduced system first selects a maximal independent row
subset.  A Phase-I program with artificial slack variables decides
feasibility, so iteration stagnation or a singular matrix is never taken
as proof of infeasibility.
"""
from itertools import combinations as _combinations, product


class MaxIterations(Exception):
    """Raised when the active-set sweep exhausts its iteration budget."""


class InfeasibleError(Exception):
    """Raised when Phase I cannot satisfy the equality system within bounds."""


def _scale(matrix):
    best = 0.0
    for row in matrix:
        for value in row:
            magnitude = abs(value)
            if magnitude > best:
                best = magnitude
    return best or 1.0


def solve_linear(matrix, rhs):
    """Gauss-Jordan elimination with partial pivoting for a dense system."""
    size = len(matrix)
    work = [list(matrix[i]) + [rhs[i]] for i in range(size)]
    tolerance = 1e-14 * max(size, 10) * _scale(matrix)
    for column in range(size):
        pivot = max(range(column, size), key=lambda r: abs(work[r][column]))
        if abs(work[pivot][column]) <= tolerance:
            raise MaxIterations("unusable pivot in linear solve")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_row = work[column]
        factor_row = pivot_row[column]
        for row_index in range(size):
            if row_index == column or work[row_index][column] == 0.0:
                continue
            factor = work[row_index][column] / factor_row
            target = work[row_index]
            for j in range(column, size + 1):
                target[j] -= factor * pivot_row[j]
    return [work[i][size] / work[i][i] for i in range(size)]


def independent_rows(rows, column_count):
    """Original row indices of a maximal linearly independent subset."""
    if not rows:
        return []
    work = [list(row) for row in rows]
    tags = list(range(len(rows)))
    tolerance = 1e-12 * _scale(rows)
    rank = 0
    chosen = []
    for column in range(column_count):
        if rank >= len(work):
            break
        pivot = max(range(rank, len(work)), key=lambda r: abs(work[r][column]))
        if abs(work[pivot][column]) <= tolerance:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        tags[rank], tags[pivot] = tags[pivot], tags[rank]
        pivot_row = work[rank]
        factor_row = pivot_row[column]
        for row_index in range(len(work)):
            if row_index == rank or work[row_index][column] == 0.0:
                continue
            factor = work[row_index][column] / factor_row
            target = work[row_index]
            for j in range(column, column_count):
                target[j] -= factor * pivot_row[j]
        chosen.append(tags[rank])
        rank += 1
    return chosen


def _direction(D, q, A, b, kept, movable, fixed, n, z):
    """Newton step toward the equality-constrained target.

    Only indices in ``movable`` may change; indices in ``fixed`` are held
    at their current value and their columns enter the equality RHS.
    """
    k = len(kept)
    if movable:
        total = len(movable) + k
        kkt = [[0.0] * total for _ in range(total)]
        for position, variable in enumerate(movable):
            kkt[position][position] = D[variable]
        for t, row_index in enumerate(kept):
            for position, variable in enumerate(movable):
                value = A[row_index][variable]
                # Canonical target KKT: D z* + A'mu = q, A z* = b.
                kkt[position][len(movable) + t] = value
                kkt[len(movable) + t][position] = value
        rhs = [q[i] for i in movable]
        rhs.extend(b[row_index] - sum(A[row_index][i] * z[i] for i in fixed)
                   for row_index in kept)
        solution = solve_linear(kkt, rhs)
        target_movable = solution[:len(movable)]
        # lambda = -mu, so stationarity reads grad - A'lambda = 0.
        multipliers = [-value for value in solution[len(movable):]]
    else:
        target_movable = []
        multipliers = []
    step = [0.0] * n
    for position, variable in enumerate(movable):
        step[variable] = target_movable[position] - z[variable]
    return step, multipliers


def _reduced_system(A, b, m, movable, n=None):
    sub_rows = [[A[r][i] for i in movable] for r in range(m)]
    return independent_rows(sub_rows, len(movable)) if movable else []


def _active_set(D, q, A, b, lo, hi, start, initial_status=None):
    """Phase II sweep. ``start`` must satisfy the equalities and bounds."""
    n = len(D)
    m = len(A)
    z = [float(v) for v in start]
    if initial_status is None:
        status = [0] * n
    else:
        status = list(initial_status)
        for i in range(n):
            if status[i] == -1:
                z[i] = lo[i]
            elif status[i] == 1:
                z[i] = hi[i]
    limit = 400 * (n + 1) + 200
    for _ in range(limit):
        free = [i for i in range(n) if status[i] == 0]
        pinned = [i for i in range(n) if status[i] != 0]
        kept = _reduced_system(A, b, m, free)
        step, multipliers = _direction(D, q, A, b, kept, free, pinned, n, z)
        step_scale = max((abs(v) for v in step), default=0.0)
        z_scale = max(1.0, max((abs(v) for v in z), default=0.0))
        if step_scale <= 1e-12 * z_scale:
            gradient = [D[i] * z[i] - q[i] for i in range(n)]
            release = _select_release(
                D, q, A, b, m, status, gradient, multipliers, kept,
                free, pinned, n, z, z_scale)
            if release is None:
                return z
            status[release] = 0
            continue
        alpha = 1.0
        blocking = None
        for variable in free:
            move = step[variable]
            if move < 0.0:
                candidate = (lo[variable] - z[variable]) / move
            elif move > 0.0:
                candidate = (hi[variable] - z[variable]) / move
            else:
                continue
            if candidate < alpha:
                alpha = candidate
                blocking = variable
        for variable in free:
            z[variable] += alpha * step[variable]
        if blocking is not None:
            if step[blocking] < 0.0:
                z[blocking] = lo[blocking]
                status[blocking] = -1
            else:
                z[blocking] = hi[blocking]
                status[blocking] = 1
    raise MaxIterations("active-set sweep did not converge")


def _select_release(D, q, A, b, m, status, gradient, multipliers, kept,
                    free, pinned, n, z, z_scale):
    """Pick a pinned variable whose freeing gives a feasible descent step.

    Strong KKT multiplier signs are tried before the degenerate near-zero
    residuals produced by redundant equality rows. Every candidate is
    confirmed by recomputing the projected direction so a release only
    happens when the freed variable can actually move.
    """
    candidates = []
    for i in pinned:
        transposed = 0.0
        for t, row_index in enumerate(kept):
            transposed += multipliers[t] * A[row_index][i]
        candidates.append((i, status[i], gradient[i] - transposed))
    for i, pin, residual in candidates:
        if pin == -1 and residual < -1e-11 * z_scale:
            if _feasible_descent(D, q, A, b, m, status, i, pin, n, z,
                                 z_scale):
                return i
        elif pin == 1 and residual > 1e-11 * z_scale:
            if _feasible_descent(D, q, A, b, m, status, i, pin, n, z,
                                 z_scale):
                return i
    for i, pin, residual in candidates:
        if abs(residual) <= 1e-11 * z_scale:
            if _feasible_descent(D, q, A, b, m, status, i, pin, n, z,
                                 z_scale):
                return i
    return None


def _feasible_descent(D, q, A, b, m, status, release, pin, n, z, z_scale):
    """True when freeing ``release`` yields a feasible descending step."""
    trial = list(status)
    trial[release] = 0
    movable = [i for i in range(n) if trial[i] == 0]
    fixed = [i for i in range(n) if trial[i] != 0]
    kept = _reduced_system(A, b, m, movable)
    step, _ = _direction(D, q, A, b, kept, movable, fixed, n, z)
    move = step[release]
    if pin == -1:
        return move > 1e-11 * z_scale
    return move < -1e-11 * z_scale


def _combinatorial_feasible(A, b, lo, hi):
    """Find a bounded feasible point by enumerating bound configurations.

    For every choice of which variables sit at a bound, the remaining
    free variables solve the equality system directly. Used as a Phase-I
    fallback when the smooth active set stalls at a degenerate vertex.
    """
    n = len(lo)
    m = len(A)
    independent = independent_rows(A, n)
    rank = len(independent)
    if rank > n:
        return None

    def search(free, fixed):
        if not free:
            if all(abs(b[r] - sum(
                    A[r][i] * (lo[i] if fixed[i] == -1 else hi[i])
                    for i in range(n) if fixed[i] != 0)) <= 1e-8
                    for r in independent):
                return [lo[i] if fixed[i] == -1 else hi[i]
                        for i in range(n)]
            return None
        sub = [[A[r][i] for i in free] for r in independent]
        kept = independent_rows(sub, len(free))
        if not kept:
            return None
        adjusted = [b[r] - sum(A[r][i] * (lo[i] if fixed[i] == -1 else hi[i])
                               for i in range(n) if fixed[i] != 0)
                    for r in independent]
        block = [[sub[t][p] for p in range(len(free))] for t in kept]
        # Under-determined systems get the minimum-norm solution through
        # the square KKT form [I, A'; A, 0][w, mu] = [0, b].
        f_count = len(free)
        k = len(kept)
        kkt = [[0.0] * (f_count + k) for _ in range(f_count + k)]
        for p in range(f_count):
            kkt[p][p] = 1.0
        for t in range(k):
            for p in range(f_count):
                value = block[t][p]
                kkt[p][f_count + t] = value
                kkt[f_count + t][p] = value
        rhs = [0.0] * f_count + adjusted
        try:
            solution = solve_linear(kkt, rhs)
        except MaxIterations:
            return None
        values = solution[:f_count]
        point = [0.0] * n
        for i in range(n):
            if fixed[i] == -1:
                point[i] = lo[i]
            elif fixed[i] == 1:
                point[i] = hi[i]
        for p, i in enumerate(free):
            value = values[p]
            if value < lo[i] - 1e-9 or value > hi[i] + 1e-9:
                return None
            point[i] = min(hi[i], max(lo[i], value))
        # Check every original equality, including rows that were linearly
        # dependent on the free-variable sub-matrix.
        for r in range(m):
            if abs(sum(A[r][i] * point[i] for i in range(n)) - b[r]) > 1e-8:
                return None
        return point

    # The combinatorial search is only tractable on small supports, where
    # the smooth active set can stall at a tight degenerate vertex. Large
    # problems always have enough free variables for the active set.
    if n > 12:
        return None
    # Prefer the interior (all free), then pin small variable subsets.
    # With at most four control equations, a vertex is determined by no
    # more than ``rank`` pins, so larger sets cannot be necessary.
    fixed = [0] * n
    result = search(list(range(n)), fixed)
    if result is not None:
        return result
    work = list(range(n))
    budget = 20000
    for count in range(1, min(rank, n) + 1):
        for combo in _combinations(work, count):
            for signs in product((-1, 1), repeat=count):
                budget -= 1
                if budget < 0:
                    return None
                trial = [0] * n
                for index, sign in zip(combo, signs):
                    trial[index] = sign
                free = [i for i in range(n) if trial[i] == 0]
                result = search(free, trial)
                if result is not None:
                    return result
    return None


def phase_one(A, b, lo, hi, denominators):
    """Find a bounded point satisfying A z = b, or raise InfeasibleError."""
    n = len(lo)
    m = len(A)
    artificial = 1.0
    for j in range(m):
        bound_sum = 0.0
        for i in range(n):
            bound_sum += max(abs(A[j][i] * lo[i]), abs(A[j][i] * hi[i]))
        artificial = max(artificial, 2.0 * abs(b[j]),
                         bound_sum + abs(b[j]) + 1.0)
    size = n + m
    extended = [list(A[j]) + [1.0 if j == t else 0.0 for t in range(m)]
                for j in range(m)]
    lo_extended = list(lo) + [-artificial] * m
    hi_extended = list(hi) + [artificial] * m
    # Pin every original variable at its lower bound and choose each
    # artificial slack so the extended equality holds exactly at start.
    start = list(lo)
    for j in range(m):
        slack = b[j] - sum(A[j][i] * lo[i] for i in range(n))
        start.append(slack)
    # Phase I minimizes the squared artificial slacks. Original variables
    # carry only a negligible strict-convexity regularizer so the program
    # stays uniquely solvable while feasibility dominates the objective.
    epsilon = 1e-10
    diagonal = [epsilon] * n + [1.0] * m
    initial_status = [-1] * n + [0] * m
    try:
        point = _active_set(diagonal, [0.0] * size, extended, b,
                            lo_extended, hi_extended, start, initial_status)
    except MaxIterations:
        point = None
    if point is not None:
        candidate = point[:n]
        converged = all(
            abs(sum(A[j][i] * candidate[i] for i in range(n)) - b[j])
            <= 1e-7 * denominators[j] for j in range(m))
    else:
        converged = False
        candidate = None
    if not converged:
        # Degenerate vertex fallback: enumerate bound configurations.
        fallback = _combinatorial_feasible(A, b, lo, hi)
        if fallback is None:
            raise InfeasibleError(
                "phase I could not meet the equality system")
        candidate = fallback
    return candidate


def bounded_minimum(D, q, A, b, lo, hi, denominators):
    """Return the unique optimizer, or raise MaxIterations after restarts."""
    lo = [float(v) for v in lo]
    hi = [float(v) for v in hi]
    candidate = phase_one(A, b, lo, hi, denominators)
    for _attempt in range(4):
        try:
            return _active_set(D, q, A, b, lo, hi, candidate)
        except MaxIterations:
            # Restart from the same feasible point with an empty working set.
            candidate = [min(hi[i], max(lo[i], candidate[i]))
                         for i in range(len(lo))]
    raise MaxIterations("repeated active-set restarts failed")
