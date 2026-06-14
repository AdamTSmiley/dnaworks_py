"""
Melting temperature (Tm) calculation for oligonucleotide overlaps.

Thin wrapper around BioPython's Bio.SeqUtils.MeltingTemp, configured
to match Benchling's "Modified Breslauer" method by default.

Benchling's "Modified Breslauer" is:
  - Breslauer 1986 NN parameters (BioPython DNA_NN1)
  - SantaLucia 1998 salt correction (saltcorr=5)
  - 50 mM Na+, 0 mM Mg2+, 250 nM total primer
  - No initiation correction (Breslauer convention)
"""

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

_NN_TABLES = {
    "breslauer86": mt.DNA_NN1,
    "santalucia97": mt.DNA_NN3,
}


def tm(
    seq: str,
    method: str = "breslauer86",
    na: float = 50.0,
    mg: float = 0.0,
    dnac: float = 250.0,
    saltcorr: int = 5,
) -> float:
    """
    Calculate Tm for a DNA sequence.

    Parameters
    seq : str
        DNA sequence (ACGT only, case-insensitive).
    method : str
        'breslauer86' (default, matches Benchling) or 'santalucia97'.
    na : float
        Sodium concentration in mM. Default 50.
    mg : float
        Magnesium concentration in mM. Default 0.
    dnac : float
        Total primer concentration in nM. Default 250.
    saltcorr : int
        Salt correction method (BioPython convention). Default 5
        (SantaLucia 1998). See BioPython docs for options 1-7.

    Returns
    float
        Melting temperature in °C. Returns 0.0 for sequences <= 7 nt.
    """
    seq = seq.upper().strip()
    if len(seq) <= 7:
        return 0.0

    nn_table = _NN_TABLES.get(method)
    if nn_table is None:
        raise ValueError(
            f"Unknown method '{method}'. Use 'breslauer86' or 'santalucia97'."
        )

    return mt.Tm_NN(
        Seq(seq),
        nn_table=nn_table,
        saltcorr=saltcorr,
        Na=na,
        Mg=mg,
        dnac1=dnac / 2.0,
        dnac2=dnac / 2.0,
    )
