"""OECD (2021) PFAS definition gate + OECD/UNEP (2018) structure-category assignment.

Definition gate (ENV/CBC/MONO(2021)25):
  a substance is a PFAS if it contains at least one FULLY FLUORINATED, SATURATED
  methyl (-CF3) or methylene (-CF2-) carbon with no H/Cl/Br/I attached to it.
  >CF- alone does not qualify; =CF2 (unsaturated) does not qualify.

Categories follow the 2018 OECD/UNEP structure series, assigned by parsing the
molecule as CnF2n+1-X-R: the fluorinated carbon moiety, the linkage X, and the
functional group R.
"""
import sys
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

SRC, OUT = sys.argv[1], sys.argv[2]
HAL = {9: "F", 17: "Cl", 35: "Br", 53: "I"}

def pf_carbons(mol):
    """Indices of fully fluorinated saturated CF3/CF2 carbons (OECD gate)."""
    out = set()
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C" or a.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        if a.GetTotalNumHs() != 0:
            continue
        nb = [n.GetSymbol() for n in a.GetNeighbors()]
        nF = nb.count("F")
        if any(x in ("Cl", "Br", "I") for x in nb):
            continue                      # fully halogenated but not fully fluorinated
        heavy = len([n for n in a.GetNeighbors() if n.GetSymbol() != "F"])
        if (nF == 3 and heavy == 1) or (nF == 2 and heavy == 2):
            out.add(a.GetIdx())
    return out

def f_carbons(mol):
    """Carbons bearing F with no H and no Cl/Br/I - used for CONNECTIVITY, not the gate.
    Includes >CF- (e.g. the GenX alpha carbon, perfluorocycloalkane ring carbons),
    which the OECD gate itself excludes."""
    out = set()
    for a in mol.GetAtoms():
        if a.GetSymbol() != "C" or a.GetTotalNumHs() != 0:
            continue
        nb = [n.GetSymbol() for n in a.GetNeighbors()]
        if any(x in ("Cl", "Br", "I") for x in nb) or "F" not in nb:
            continue
        out.add(a.GetIdx())
    return out

def attached(mol, fc, syms, exclude_arom=True):
    """Is any fluorinated carbon bonded to an atom of one of `syms`?"""
    for i in fc:
        for n in mol.GetAtomWithIdx(i).GetNeighbors():
            if n.GetSymbol() in syms and not (exclude_arom and n.GetIsAromatic()):
                return n
    return None

def chain_stats(mol, pfc):
    """Longest contiguous run of perfluorinated carbons, and whether ether-interrupted."""
    if not pfc: return 0, False
    adj = {i: set() for i in pfc}
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in pfc and j in pfc:
            adj[i].add(j); adj[j].add(i)
    best, seen = 0, set()
    for s in pfc:
        if s in seen: continue
        comp, stack = set(), [s]
        while stack:
            u = stack.pop()
            if u in comp: continue
            comp.add(u); stack.extend(adj[u] - comp)
        seen |= comp; best = max(best, len(comp))
    return best, None

Q = {k: Chem.MolFromSmarts(v) for k, v in {
    # head groups attached directly to a perfluorinated carbon
    "pf_carboxyl":  "[CX4;!$([CX4][H])](F)(F)[CX3](=O)[OX2H1,OX1-,OX2]",
    "pf_sulfonyl":  "[CX4](F)(F)[SX4](=O)(=O)",
    "pf_sulfonamide":"[CX4](F)(F)[SX4](=O)(=O)[NX3]",
    "pf_phosphorus":"[CX4](F)(F)[PX4]",
    "ft_n2":        "[CX4](F)(F)[CH2][CH2]",
    "ft_n1":        "[CX4](F)(F)[CH2;!$([CH2][CH2]);!$([CH2][CH])]",
    "pf_ether":     "[CX4](F)(F)[OX2][CX4](F)(F)",
    "aromatic_sc":  "c~*~[CX4](F)(F)F",
    "aromatic_dir": "c[CX4](F)(F)F",
}.items()}
assert all(v is not None for v in Q.values())

def has(m, k): return m.HasSubstructMatch(Q[k])

def classify(smi):
    m = Chem.MolFromSmiles(str(smi))
    if m is None:
        return dict(oecd_pfas=False, oecd_class="unparsable", n_pfC=0, pf_chain=0, ether=False)
    pfc = pf_carbons(m)
    fc = f_carbons(m)
    n, _ = chain_stats(m, pfc)
    # ether: O bridging two fluorinated carbons, or a fluorinated carbon and an aryl O
    ether = False
    for o in m.GetAtoms():
        if o.GetSymbol() != "O" or o.GetDegree() != 2:
            continue
        nb = list(o.GetNeighbors())
        nf = sum(1 for x in nb if x.GetIdx() in fc)
        if nf >= 2 or (nf == 1 and any(x.GetIsAromatic() for x in nb)):
            ether = True; break
    if not pfc:
        return dict(oecd_pfas=False, oecd_class="not PFAS (OECD 2021)",
                    n_pfC=0, pf_chain=0, ether=ether)

    # ---- category cascade: head group reached from a fluorinated carbon ----
    sulfonyl   = attached(m, fc, ("S",))
    phosphorus = attached(m, fc, ("P",))
    def _is_acid(c):
        """Carboxylic acid or carboxylate: C(=O) bearing a second O with H or a negative charge."""
        dbl = any(b.GetBondType() == Chem.BondType.DOUBLE and
                  b.GetOtherAtom(c).GetSymbol() == "O" for b in c.GetBonds())
        if not dbl:
            return False
        for b in c.GetBonds():
            o = b.GetOtherAtom(c)
            if (o.GetSymbol() == "O" and b.GetBondType() == Chem.BondType.SINGLE
                    and (o.GetTotalNumHs() > 0 or o.GetFormalCharge() < 0)):
                return True
        return False

    def _carbonyl_from(i, depth):
        """Walk out from a fluorinated carbon, optionally crossing one sp2 carbon.

        An alpha,beta-unsaturated perfluoroalkenoic acid carries its acid carbon two
        bonds from the nearest fluorinated carbon, across a =CH- that bears hydrogen
        and is therefore not itself fluorinated; a one-bond walk misses it entirely.
        """
        for nb in m.GetAtomWithIdx(i).GetNeighbors():
            if nb.GetSymbol() != "C" or nb.GetIsAromatic():
                continue
            if _is_acid(nb):
                return nb, True
            if nb.GetHybridization() == Chem.HybridizationType.SP2 and depth > 0:
                hit, acid = _carbonyl_from(nb.GetIdx(), depth - 1)
                if hit is not None:
                    return hit, acid
            # a carbonyl that is not an acid: ketone, ester or acyl halide
            if any(b.GetBondType() == Chem.BondType.DOUBLE and
                   b.GetOtherAtom(nb).GetSymbol() == "O" for b in nb.GetBonds()):
                return nb, False
        return None, False

    carbonylC, carbonyl_is_acid = None, False
    for i in fc:
        hit, acid = _carbonyl_from(i, 1)
        if hit is not None and (carbonylC is None or acid):
            carbonylC, carbonyl_is_acid = hit, acid
            if acid:
                break
    is_sulfonamide = bool(sulfonyl) and any(
        x.GetSymbol() == "N" for x in sulfonyl.GetNeighbors()) if sulfonyl else False
    # Fluorotelomer-based substances are CnF2n+1-CH2CH2-R products of telomerisation;
    # the 2018 OECD/UNEP series enumerates n = 4, 6, 8, 10. A C2F5 or C3F7 terminus on
    # an otherwise non-fluorinated scaffold satisfies the SMARTS but is not a telomer
    # product, so the chain-length floor is set at n >= 4 and shorter matches fall
    # through to assignment by head group.
    telomer = (n >= 4) and (has(m, "ft_n2") or has(m, "ft_n1"))
    # 700 requires a fluorinated side chain ON an aromatic ring, not merely the
    # presence of a ring somewhere in the molecule.
    aromatic = False
    for i in fc:
        stack, seen = [(i, 0)], {i}
        while stack:
            j, dist = stack.pop()
            if dist > 2:
                continue
            for nb in m.GetAtomWithIdx(j).GetNeighbors():
                if nb.GetIsAromatic():
                    aromatic = True; stack = []; break
                if nb.GetIdx() not in seen:
                    seen.add(nb.GetIdx()); stack.append((nb.GetIdx(), dist + 1))
            if aromatic:
                break
        if aromatic:
            break

    if telomer:
        c = "400 fluorotelomer"
    elif ether and ((carbonylC is not None and carbonyl_is_acid) or sulfonyl is not None):
        c = "500 perfluoroalkyl ether acid"
    elif is_sulfonamide:
        c = "200b perfluoroalkane sulfonamide"
    elif sulfonyl is not None:
        c = "200a perfluoroalkane sulfonic"
    elif carbonylC is not None and carbonyl_is_acid:
        c = "100 perfluoroalkyl carboxylic"
    elif phosphorus is not None:
        c = "300 perfluoroalkyl phosphorus"
    elif ether:
        c = "500 perfluoroalkyl ether other"
    elif aromatic:
        c = "700 side-chain fluorinated aromatic"
    else:
        c = "600 other perfluoroalkyl"
    return dict(oecd_pfas=True, oecd_class=c, n_pfC=len(pfc), pf_chain=n, ether=ether)

d = pd.read_csv(SRC)
res = pd.DataFrame([classify(s) for s in d.canonical_smiles])
out = pd.concat([d.reset_index(drop=True), res], axis=1)
out.to_csv(OUT, index=False)

print(f"n = {len(out)}   OECD-2021 PFAS: {out.oecd_pfas.sum()}   not PFAS: {(~out.oecd_pfas).sum()}\n")
print("=== OECD structure categories ===")
t = (out.groupby("oecd_class")
        .agg(n=("logTox","size"), median_logTox=("logTox","median"),
             median_chain=("pf_chain","median"), median_pfC=("n_pfC","median"))
        .sort_values("n", ascending=False))
print(t.to_string(float_format=lambda x: f"{x:6.2f}"))

print("\n=== NOT PFAS under OECD 2021 (dropped from the PFAS set) ===")
for _, r in out[~out.oecd_pfas].iterrows():
    print(f"   {str(r['name'])[:70]}")
