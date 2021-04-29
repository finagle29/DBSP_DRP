from dbsp_drp import reduction

def test_parse_param_file(tmp_path):
    red_lines = ['a = 1.0\n','b = [3, 4] \n','c=five\n']
    param_fname = tmp_path / 'test.params'
    with open(param_fname, 'w') as f:
        f.write('[blue]\n')
        f.write('a=0\n')
        f.write('[green]\n')
        f.write('a = 1\n')
        f.write('[red]\n')
        f.writelines(red_lines)
        f.write('[orange]\n')
        f.write('a = blue\n')

    with open(param_fname) as f:
        print(f.readlines())

    result = reduction.parse_pypeit_parameter_file(param_fname, 'red', ['red', 'orange', 'green', 'blue'])
    assert result == red_lines
