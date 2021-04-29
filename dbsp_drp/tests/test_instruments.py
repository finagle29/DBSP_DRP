import pytest

from pypeit.spectrographs.util import load_spectrograph

from dbsp_drp import instruments

def instrument_tester(ins_class: instruments.Instrument):
    ins = ins_class()
    assert ins.arms == len(ins.archived_sensfuncs)
    assert ins.arms == len(ins.arm_names_pypeit)
    assert ins.arms == len(ins.arm_prefixes)
    assert ins.arms == len(ins.arm_telluric)
    assert ins.arms == len(ins.coadd_threshholds)
    assert ins.arms == len(ins.frac_pos_threshholds)
    assert ins.arms == len(ins.pixel_pos_threshholds)
    assert ins.arms == len(ins.pypeit_name_to_arm)

    for i in range(ins.arms):
        assert ins.arm_prefixes[i] in ins.archived_sensfuncs[i]
        load_spectrograph(ins.arm_names_pypeit[i])
        assert ins.arm_prefixes[i] in ins.arm_names_pypeit[i]
        assert ins.pypeit_name_to_arm[ins.arm_names_pypeit[i]] == ins.arm_prefixes[i]

def test_instrument():
    with pytest.raises(TypeError):
        instrument_tester(instruments.Instrument)

def test_dbsp():
    instrument_tester(instruments.DBSP)
