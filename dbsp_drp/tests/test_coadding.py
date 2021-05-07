from dbsp_drp import coadding

def test_group_coadds():
    fname_to_spats = {
        'a1': [100, 200, 300],
        'a2': [100, 200, 300],
        'a3': [101, 199],
        'a4': [99, 201, 303]
    }

    correct = [
        {
            'spats': [99, 100, 100, 101],
            'fnames': ['a4', 'a1', 'a2', 'a3']
        },
        {
            'spats': [199, 200, 200, 201],
            'fnames': ['a3', 'a1', 'a2', 'a4']
        },
        {
            'spats': [300, 300],
            'fnames': ['a1', 'a2']
        },
        {
            'spats': [303],
            'fnames': ['a4']
        }
    ]

    grouped_coadds = coadding.group_coadds(fname_to_spats)

    assert correct == grouped_coadds
