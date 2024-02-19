from load_data import load_data

def test_load_data():
    df = load_data('rgi_strict', 'human gut')
    assert df.shape[0] == 6729
    df = load_data('rgi_strict', 'dog gut')
    assert df.shape[0] == 129
    df = load_data('rgi', 'dog gut')
    assert df.shape[0] == 129

    # deeparg requires special casing so add test for it
    df = load_data('deeparg', 'dog gut')
    assert df.shape[0] == 129
