"""Minimal SMILES reader used when RDKit is not installed: elemental composition, formula and
monoisotopic / average mass. It is not RDKit; it covers common organic SMILES (organic subset,
bracket atoms, rings, branches, aromatic atoms, charges)."""
import re

MONO = {"H": 1.00782503207, "B": 11.0093054, "C": 12.0, "N": 14.0030740048, "O": 15.99491461956,
        "F": 18.99840322, "Na": 22.9897692809, "Si": 27.9769265325, "P": 30.97376163, "S": 31.97207100,
        "Cl": 34.96885268, "K": 38.96370668, "Br": 78.9183371, "I": 126.904473}
AVG = {"H": 1.008, "B": 10.811, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Na": 22.990,
       "Si": 28.086, "P": 30.974, "S": 32.065, "Cl": 35.453, "K": 39.098, "Br": 79.904, "I": 126.904}
VALENCES = {"B": [3], "C": [4], "N": [3, 5], "O": [2], "P": [3, 5], "S": [2, 4, 6],
            "F": [1], "Cl": [1], "Br": [1], "I": [1], "Si": [4]}
ELECTRON = 0.000548579909

_TOKEN = re.compile(r"(\[[^\]]+\]|Br|Cl|[BCNOPSFI]|[bcnops]|%\d\d|\d|[=#$:/\\\-.()])")
_BRACKET = re.compile(r"^\[(\d+)?([A-Z][a-z]?|[a-z]{1,2})(@{1,2}(?:TH\d|AL\d|SP\d|TB\d+|OH\d+)?)?(H\d*)?([+-]+\d*)?(:\d+)?\]$")
_BOND = {"-": 1.0, "=": 2.0, "#": 3.0, "$": 4.0, ":": 1.5, "/": 1.0, "\\": 1.0}


class Atom:
    def __init__(self, idx, symbol, aromatic, bracket, hcount, charge):
        self.idx, self.symbol, self.aromatic, self.bracket = idx, symbol, aromatic, bracket
        self.hcount, self.charge = hcount, charge
        self.bonds = []          # bond orders (aromatic = 1.5)
        self.radicals = 0
        self.no_implicit = bracket

    def GetSymbol(self): return self.symbol
    def GetIdx(self): return self.idx
    def GetNumRadicalElectrons(self): return self.radicals
    def SetNumRadicalElectrons(self, n): self.radicals = n
    def SetNoImplicit(self, flag): self.no_implicit = flag
    def GetFormalCharge(self): return self.charge
    def GetIsAromatic(self): return self.aromatic

    def bond_sum(self):
        # an aromatic bond counts 1, plus 1 for taking part in the aromatic system
        n_arom = sum(1 for b in self.bonds if b == 1.5)
        s = sum(1 for b in self.bonds if b == 1.5) + sum(b for b in self.bonds if b != 1.5)
        return s + (1 if n_arom else 0)

    def _valence_for(self, used):
        vals = VALENCES.get(self.symbol, [0])
        if self.charge:   # isoelectronic: N+ ~ C (4), O- ~ F (1), O+ ~ N (3)...
            iso = {("N", 1): [4], ("N", -1): [2], ("O", 1): [3], ("O", -1): [1], ("C", -1): [3], ("C", 1): [3],
                   ("S", 1): [3, 5], ("S", -1): [1, 3, 5], ("P", 1): [4]}
            vals = iso.get((self.symbol, self.charge), vals)
        for v in vals:
            if v >= used - 1e-9:
                return v
        return vals[-1]

    def total_h(self):
        if self.aromatic and not self.bracket:
            # organic-subset aromatic atoms: c and b take one H if a valence is left after the ring bonds
            # (each aromatic bond counts 1, plus 1 for the pi system); n, o, s, p take none (pyrrole-type
            # nitrogen must be written [nH], as in RDKit)
            if self.symbol in ("C", "B"):
                n_arom = sum(1 for b in self.bonds if b == 1.5)
                used = n_arom + sum(b for b in self.bonds if b != 1.5) + 1
                return max(0, int(round(VALENCES[self.symbol][0] - used)))
            return 0
        if not self.no_implicit:      # organic subset, or a "corrected" radical
            used = self.bond_sum() + (self.hcount if self.bracket else 0)
            return (self.hcount if self.bracket else 0) + max(0, int(round(self._valence_for(used) - used)))
        return self.hcount


class Mol:
    def __init__(self, atoms):
        self.atoms = atoms

    def GetAtoms(self): return list(self.atoms)
    def GetNumAtoms(self): return len(self.atoms)


def _parse_bracket(tok):
    m = _BRACKET.match(tok)
    if not m:
        raise ValueError(tok)
    sym = m.group(2)
    aromatic = sym[0].islower()
    sym = sym.capitalize() if aromatic else sym
    h = m.group(4)
    hcount = 0 if not h else (1 if h == "H" else int(h[1:]))
    ch = m.group(5) or ""
    if ch:
        sign = 1 if ch[0] == "+" else -1
        digits = ch.lstrip("+-")
        charge = sign * (int(digits) if digits else len(ch))
    else:
        charge = 0
    return sym, aromatic, hcount, charge


def parse(smi):
    try:
        toks = _TOKEN.findall(smi)
        if "".join(toks) != smi:
            return None
        atoms, stack, rings = [], [], {}
        prev, bond = None, None
        for t in toks:
            if t == "(":
                stack.append(prev); continue
            if t == ")":
                prev = stack.pop(); continue
            if t == ".":
                prev = None; continue
            if t in _BOND:
                bond = _BOND[t]; continue
            if t[0] == "%" or t.isdigit():
                n = int(t.lstrip("%"))
                if n in rings:
                    other, b0 = rings.pop(n)
                    b = bond or b0 or (1.5 if (prev.aromatic and other.aromatic) else 1.0)
                    prev.bonds.append(b); other.bonds.append(b)
                else:
                    rings[n] = (prev, bond)
                bond = None
                continue
            if t[0] == "[":
                sym, aromatic, h, ch = _parse_bracket(t)
                a = Atom(len(atoms), sym, aromatic, True, h, ch)
            else:
                aromatic = t[0].islower()
                a = Atom(len(atoms), t.capitalize() if aromatic else t, aromatic, False, 0, 0)
            atoms.append(a)
            if prev is not None:
                b = bond or (1.5 if (prev.aromatic and a.aromatic) else 1.0)
                prev.bonds.append(b); a.bonds.append(b)
            prev, bond = a, None
        if rings or stack:
            return None
        # radicals: bracket atom, not aromatic, with incomplete valence
        for a in atoms:
            if a.bracket and not a.aromatic and a.symbol in VALENCES:
                used = a.bond_sum() + a.hcount
                v = a._valence_for(used)
                a.radicals = max(0, int(round(v - used)))
        return Mol(atoms)
    except Exception:
        return None


def _sanitize(m):
    return m


def composition(m):
    comp = {}
    for a in m.atoms:
        comp[a.symbol] = comp.get(a.symbol, 0) + 1
        h = a.total_h()
        if h:
            comp["H"] = comp.get("H", 0) + h
    return comp


def charge(m):
    return sum(a.charge for a in m.atoms)


def exact_mass(m) -> float:
    return sum(MONO[e] * n for e, n in composition(m).items()) - charge(m) * ELECTRON


def average_mass(m) -> float:
    return sum(AVG[e] * n for e, n in composition(m).items())


def formula(m) -> str:
    comp = composition(m)
    order = (["C", "H"] if "C" in comp else []) + sorted(e for e in comp if e not in ("C", "H") or "C" not in comp)
    s = "".join(f"{e}{comp[e] if comp[e] > 1 else ''}" for e in order)
    q = charge(m)
    if q:
        s += ("+" if q > 0 else "-") + (str(abs(q)) if abs(q) > 1 else "")
    return s
