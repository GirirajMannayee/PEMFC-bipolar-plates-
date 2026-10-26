from adc12pem import verify


def test_all_manuscript_numbers_reproduce():
    df = verify.run_checks()
    bad = df[~df.passed]
    assert bad.empty, "\n" + bad.to_string(index=False)


def test_scorecard_best_is_s4():
    from adc12pem.io import load
    _, norm, comp = verify.scorecard(*(load(n) for n in ["table06_conductivity_hardness", "table07_eis_25C",
                                                           "table08_polarisation_lpr", "table09_icr_vs_pressure"]))
    assert comp.idxmax() == "S4" and abs(comp["S4"] - 1) < 1e-12
