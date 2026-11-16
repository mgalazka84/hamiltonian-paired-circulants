"""Reproduce the finite checks in 'Hamiltonian-paired circulant graphs of valency four'.

Run 'python verify.py census' or 'python verify.py witness 13 1 5'.
The Hamilton-cycle search uses only the adjacency list of the input graph.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
from math import gcd


Edge = tuple[int, int]
Cycle = tuple[int, ...]


def edge(v: int, w: int) -> Edge:
    return (v, w) if v < w else (w, v)


def magnitude(value: int, n: int) -> int:
    residue = value % n
    return min(residue, n - residue)


def check_input(n: int, a: int, b: int) -> None:
    if n < 5 or not (0 < 2 * a < n and 0 < 2 * b < n and a != b):
        raise ValueError("Use n >= 5 and two distinct integers strictly between 0 and n/2.")
    if gcd(n, gcd(a, b)) != 1:
        raise ValueError("The two steps must generate the whole cyclic group.")


def factor(n: int, step: int) -> set[Edge]:
    return {edge(x, (x + step) % n) for x in range(n)}


def switch(n: int, a: int, b: int, t: int, red: set[Edge], blue: set[Edge]) -> None:
    x, y, z, w = (t + offset for offset in (0, a, a + b, b))
    x, y, z, w = x % n, y % n, z % n, w % n
    horizontal = {edge(x, y), edge(w, z)}
    vertical = {edge(x, w), edge(y, z)}
    if not horizontal <= red or not vertical <= blue:
        raise AssertionError("The proposed square does not have the required colours.")
    red.difference_update(horizontal)
    blue.difference_update(vertical)
    red.update(vertical)
    blue.update(horizontal)


def read_cycles(n: int, edges: set[Edge]) -> tuple[Cycle, ...]:
    neighbors = [[] for _ in range(n)]
    for v, w in edges:
        neighbors[v].append(w)
        neighbors[w].append(v)
    if any(len(row) != 2 for row in neighbors):
        raise AssertionError("The proposed factor is not spanning and 2-regular.")

    seen = set()
    cycles = []
    for start in range(n):
        if start in seen:
            continue
        sequence = [start]
        seen.add(start)
        previous, current = start, min(neighbors[start])
        while current != start:
            if current in seen:
                raise AssertionError("A cycle repeats a vertex before closing.")
            seen.add(current)
            sequence.append(current)
            left, right = neighbors[current]
            previous, current = current, right if left == previous else left
        cycles.append(tuple(sequence))
    return tuple(cycles)


@dataclass
class Witness:
    steps: tuple[int, int]
    square_origins: tuple[int, ...]
    hamilton_cycle: Cycle
    remaining_cycles: tuple[Cycle, ...]


def construct(n: int, a: int, b: int) -> Witness | None:
    """Return a checked counterexample in the user's original vertex labels."""
    check_input(n, a, b)
    original_edges = factor(n, a) | factor(n, b)
    d, f = gcd(a, n), gcd(b, n)
    if d > f:
        a, b, d, f = b, a, f, d

    if d != f:
        origins = tuple((i * b + (i % 2) * a) % n for i in range(d - 1))
    else:
        # Connectedness makes both gcds equal to one.
        k = magnitude(b * pow(a, -1, n), n)
        if (n, k) == (5, 2) or (n % 4 == 0 and 2 * k + 2 == n):
            return None
        inverse = pow(k, -1, n)
        h = magnitude(inverse, n)
        for t in range(2, n // 2 + 1):
            if t == k:
                continue
            splits_red = t < k
            splits_blue = 2 <= magnitude(inverse * t, n) < h
            if splits_red != splits_blue:
                break
        else:
            raise AssertionError("No distinguishing shift was found.")
        # Multiplication by a returns from the normal form to the input labels.
        b = (a * k) % n
        origins = (0, (a * t) % n)

    red, blue = factor(n, a), factor(n, b)
    for t in origins:
        switch(n, a, b, t, red, blue)
    if red & blue or red | blue != original_edges:
        raise AssertionError("The two factors do not partition the input edges.")
    first, second = read_cycles(n, red), read_cycles(n, blue)
    if len(first) != 1:
        first, second = second, first
    if len(first) != 1 or len(second) < 2:
        raise AssertionError("The constructed partition is not a counterexample.")
    return Witness((a, b), origins, first[0], second)


def exhaustive_paired(n: int, a: int, b: int) -> bool:
    """Search all Hamilton cycles until a disconnected complement is found."""
    adj = [sorted({(x + s) % n for s in (a, -a, b, -b)}) for x in range(n)]
    path = [0]
    saw_hamilton_cycle = False

    def complement_disconnected() -> bool:
        remaining = [set(row) for row in adj]
        for v, w in zip(path, path[1:] + [0]):
            remaining[v].remove(w)
            remaining[w].remove(v)
        reached, stack = {0}, [0]
        while stack:
            v = stack.pop()
            for w in remaining[v] - reached:
                reached.add(w)
                stack.append(w)
        return len(reached) < n

    def has_witness(v: int, visited: int) -> bool:
        nonlocal saw_hamilton_cycle
        if len(path) == n:
            if 0 not in adj[v] or path[1] >= v:
                return False
            saw_hamilton_cycle = True
            return complement_disconnected()
        for w in adj[v]:
            if visited & (1 << w):
                continue
            path.append(w)
            found = has_witness(w, visited | (1 << w))
            path.pop()
            if found:
                return True
        return False

    found = has_witness(0, 1)
    return saw_hamilton_cycle and not found


def census(max_order: int, exhaustive_order: int) -> None:
    if max_order < 5 or not 0 <= exhaustive_order <= max_order:
        raise ValueError("Require max-order >= 5 and 0 <= exhaustive-order <= max-order.")
    total = positives = small_total = small_positives = 0
    for n in range(5, max_order + 1):
        steps = range(1, (n - 1) // 2 + 1)
        for a, b in combinations(steps, 2):
            if gcd(n, gcd(a, b)) != 1:
                continue
            paired = construct(n, a, b) is None
            total += 1
            positives += paired
            if n <= exhaustive_order:
                independent = exhaustive_paired(n, a, b)
                if independent != paired:
                    raise AssertionError(f"Disagreement for C_{n}({a},{b})")
                small_total += 1
                small_positives += independent
    print(f"Through order {max_order}: {total} parameter triples, {positives} "
          f"in the listed families, {total - positives} negative certificates checked.")
    if exhaustive_order >= 5:
        print(f"Direct search through order {exhaustive_order}: {small_total} triples, "
              f"{small_positives} paired, {small_total - small_positives} non-paired.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    count = commands.add_parser("census", help="run the finite checks")
    count.add_argument("--max-order", type=int, default=80)
    count.add_argument("--exhaustive-order", type=int, default=26)
    example = commands.add_parser("witness", help="print a negative certificate, if one exists")
    for arg in ("n", "a", "b"):
        example.add_argument(arg, type=int)
    args = parser.parse_args()
    try:
        if args.command == "census":
            census(args.max_order, args.exhaustive_order)
        else:
            result = construct(args.n, args.a, args.b)
            if result is None:
                print("Hamiltonian-paired by the classification.")
            else:
                print("Square steps:", result.steps, "Origins:", result.square_origins)
                print("Hamilton cycle:", result.hamilton_cycle)
                print("Complement cycles:", result.remaining_cycles)
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
