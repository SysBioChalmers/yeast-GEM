"""Single-gene-knockout vs the Stanford yeast deletion collection.

Port of ``code/modelTests/essentialGenes.m``. Constrains the model to
the Kennedy synthetic complete medium, runs single-gene deletion via
cobrapy, then compares the predicted essentiality against the curated
reference lists under ``data/essentialGenes/``.
"""
from __future__ import annotations

from dataclasses import dataclass

import cobra
from cobra.flux_analysis import single_gene_deletion

from yeastgem.io import REPO_PATH, read_yeast_model

# Kennedy-medium exchange-reaction bound presets, copied from the
# ``complete_Y7`` local function in essentialGenes.m.
_COMPLETE_Y7_CONSTRAINED_UPTAKE = (
    "r_1604", "r_1639", "r_1873", "r_1879", "r_1880", "r_1881", "r_1671",
    "r_1883", "r_1757", "r_1891", "r_1889", "r_1810", "r_1993", "r_1893",
    "r_1897", "r_1947", "r_1899", "r_1900", "r_1902", "r_1967",
    "r_1903", "r_1548", "r_1904", "r_2028", "r_2038", "r_1906", "r_2067",
    "r_1911", "r_1912", "r_1913", "r_2090", "r_1914", "r_2106",
)
_COMPLETE_Y7_GLUCOSE_EX = "r_1714"
_COMPLETE_Y7_UNCONSTRAINED_UPTAKE = (
    "r_1672", "r_1654", "r_1992", "r_2005", "r_2060", "r_1861", "r_1832",
    "r_2100", "r_4593", "r_4595", "r_4596", "r_4597", "r_2049", "r_4594",
    "r_4600", "r_2020",
)
_KO_TOL = 1e-6


@dataclass(frozen=True)
class EssentialGeneResult:
    """Confusion-matrix breakdown of essential-gene predictions.

    All four attribute lists are sorted, deduplicated, and intersected
    with the verified-ORF set (so the metrics are comparable across
    different model versions).
    """

    accuracy: float
    sensitivity: float
    specificity: float
    mcc: float
    tp: list[str]
    tn: list[str]
    fp: list[str]
    fn: list[str]


def essential_genes(
    model: cobra.Model | None = None,
    *,
    write_output: bool = False,
) -> EssentialGeneResult:
    """Predict essential genes + compare against Stanford deletion lists.

    Parameters
    ----------
    model
        Model to test. Defaults to a fresh :func:`read_yeast_model` load.
    write_output
        When True, save a markdown report to
        ``data/testResults/essentialGenes.md``.
    """
    if model is None:
        model = read_yeast_model()
    model = _apply_complete_Y7(model.copy())

    inviable_orfs = _load_orfs("inviable_orfs.txt")
    verified_orfs = _load_orfs("verified_orfs.txt")
    model_genes = {g.id for g in model.genes}

    exp_inviable = (model_genes & inviable_orfs) & verified_orfs
    exp_viable = (model_genes - inviable_orfs) & verified_orfs

    # Wild-type FBA pins the denominator for the growth ratio.
    wild_type = model.optimize()
    if wild_type.status != "optimal" or wild_type.objective_value <= 0:
        raise RuntimeError(
            f"Wild-type FBA returned {wild_type.status}/"
            f"obj={wild_type.objective_value}; cannot run deletion benchmark."
        )
    knockout = single_gene_deletion(model)
    gr_ratio = _knockout_growth_ratio(
        knockout, model_genes, wild_type.objective_value,
    )

    mod_viable = {gid for gid, ratio in gr_ratio.items() if ratio >= _KO_TOL}
    mod_inviable = model_genes - mod_viable
    mod_viable &= verified_orfs
    mod_inviable &= verified_orfs

    tp = sorted(exp_viable & mod_viable)
    tn = sorted(exp_inviable & mod_inviable)
    fp = sorted(exp_inviable & mod_viable)
    fn = sorted(exp_viable & mod_inviable)
    n_tp, n_tn, n_fp, n_fn = len(tp), len(tn), len(fp), len(fn)
    total = n_tp + n_tn + n_fp + n_fn
    accuracy = (n_tp + n_tn) / total if total else 0.0
    sensitivity = 100 * n_tp / (n_tp + n_fn) if (n_tp + n_fn) else 0.0
    specificity = 100 * n_tn / (n_tn + n_fp) if (n_tn + n_fp) else 0.0
    denom_mcc = ((n_tp + n_fp) * (n_tp + n_fn)
                 * (n_tn + n_fp) * (n_tn + n_fn))
    mcc = (n_tp * n_tn - n_fp * n_fn) / (denom_mcc ** 0.5) if denom_mcc else 0.0

    result = EssentialGeneResult(
        accuracy=accuracy,
        sensitivity=sensitivity,
        specificity=specificity,
        mcc=mcc,
        tp=tp, tn=tn, fp=fp, fn=fn,
    )

    if write_output:
        _write_report(result)
    return result


# --- helpers ----------------------------------------------------------

def _apply_complete_Y7(model: cobra.Model) -> cobra.Model:
    """Constrain to the Kennedy synthetic complete medium."""
    # Reset every exchange reaction to (0, 1000).
    for rxn in model.exchanges:
        rxn.lower_bound = 0
        rxn.upper_bound = 1000
    for rxn_id in _COMPLETE_Y7_CONSTRAINED_UPTAKE:
        _try_set_lb(model, rxn_id, -0.5)
    _try_set_lb(model, _COMPLETE_Y7_GLUCOSE_EX, -20)
    for rxn_id in _COMPLETE_Y7_UNCONSTRAINED_UPTAKE:
        _try_set_lb(model, rxn_id, -1000)
    return model


def _try_set_lb(model: cobra.Model, rxn_id: str, lb: float) -> None:
    try:
        model.reactions.get_by_id(rxn_id).lower_bound = lb
    except KeyError:
        pass  # quietly skip rxns missing from this model version


def _load_orfs(filename: str) -> set[str]:
    path = REPO_PATH / "data" / "essentialGenes" / filename
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def _knockout_growth_ratio(
    knockout, model_genes: set[str], wild_type_growth: float,
) -> dict[str, float]:
    """Extract per-gene growth ratio from cobrapy's deletion DataFrame.

    cobrapy's ``single_gene_deletion`` returns a 0-indexed DataFrame
    with the deleted gene id(s) in the ``ids`` column (a frozenset)
    and the resulting absolute growth in ``growth``. We divide by the
    wild-type growth to get the ratio, matching MATLAB
    ``findGeneDeletions``'s ``grRatioMuts``.
    """
    out: dict[str, float] = {}
    for _, row in knockout.iterrows():
        ids = row["ids"]
        if not isinstance(ids, (frozenset, set, tuple, list)) or len(ids) != 1:
            continue
        (gid,) = ids
        if gid not in model_genes:
            continue
        growth_val = row["growth"]
        # cobrapy sets growth to NaN for infeasible deletions; treat
        # those as zero growth (= inviable).
        if growth_val != growth_val:  # NaN check
            growth_val = 0.0
        out[gid] = float(growth_val) / float(wild_type_growth)
    # Any gene not in the deletion table defaults to wild-type viable
    # (ratio=1) — matches MATLAB findGeneDeletions, which leaves
    # unhit indices at the initialised "1" value.
    for gid in model_genes:
        out.setdefault(gid, 1.0)
    return out


def _write_report(result: EssentialGeneResult) -> None:
    out_dir = REPO_PATH / "data" / "testResults"
    out_dir.mkdir(parents=True, exist_ok=True)
    md = []
    md.append("## False non-essential genes")
    md.extend(result.fp)
    md.append("## False essential genes")
    md.extend(result.fn)
    md.append("## True non-essential genes")
    md.extend(result.tp)
    md.append("## True essential genes")
    md.extend(result.tn)
    (out_dir / "essentialGenes.md").write_text("\n".join(md) + "\n")
