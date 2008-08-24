from oldowan.mitomotifs.str2sites import str2sites

def test_empty_string():
    assert [] == str2sites('')

def test_single_substitution():
    assert 16129 == str2sites('16129A')[0].position
    assert 'A' == str2sites('16129A')[0].value

def test_single_substitution_add16k():
    assert 16129 == str2sites('129A',add16k=True)[0].position
    assert 'A' == str2sites('129A',add16k=True)[0].value
    
def test_transition_without_value():
    assert 16129 == str2sites('16129')[0].position
    assert 'C' == str2sites('16129')[0].value
    
def test_transition_without_value_add16k():
    assert 16129 == str2sites('129',True)[0].position
    assert 'C' == str2sites('129',True)[0].value
    
def test_deletion_in_minus_format():
    assert '-' == str2sites('16129-')[0].value
    
def test_deletion_in_minus_format_add16k():
    assert '-' == str2sites('129-',True)[0].value
    
def test_deletion_in_standard_format():
    assert '-' == str2sites('16129d')[0].value

def test_deletion_in_standard_format_add16k():
    assert '-' == str2sites('129d',True)[0].value

def test_deletion_in_del_format():
    assert '-' == str2sites('16129del')[0].value
    
def test_deletion_in_del_and_value_format():
    assert '-' == str2sites('16129delA')[0].value
    
def test_insertion_in_standard_format():
    assert 1 == str2sites('16129.1A')[0].insert

def test_insertion_in_standard_format_add16k():
    assert 1 == str2sites('129.1A',True)[0].insert

def test_insertion_in_ins_format():
    assert 1 == str2sites('16129insA')[0].insert

def test_insertion_in_ins_format_add16k():
    assert 1 == str2sites('129insA',True)[0].insert

def test_insertion_in_plus_format():
    assert 1 == str2sites('16129+A')[0].insert
    
def test_insertion_in_plus_format_add16k():
    assert 1 == str2sites('129+A',True)[0].insert

def test_two_insertions_in_ins_format():
    assert len(str2sites('16129insAA')) == 2

def test_two_insertions_in_ins_format_add16k():
    assert len(str2sites('129insAA',True)) == 2

def test_two_insertions_in_plus_format():
    assert len(str2sites('16129+AA')) == 2
    
def test_two_insertions_in_plus_format_add16k():
    assert len(str2sites('129+AA',True)) == 2

def test_two_insertions_in_one_position():
    assert len(str2sites('16129.1A 16129.2A')) == 2
    assert 16129 == str2sites('16129.1A 16129.2A')[0].position
    assert str2sites('16129.1A 16129.2A')[0].insert in [1,2]
    
def test_two_insertions_in_one_position_add16k():
    assert len(str2sites('129.1A 129.2A', True)) == 2
    assert 16129 == str2sites('129.1A 129.2A', True)[0].position
    assert str2sites('129.1A 129.2A', True)[0].insert in [1,2]

def test_spaces():
    assert len(str2sites('16293G 16311C')) == 2
    assert 'G' == str2sites('16293G 16311C')[0].value
    assert 'C' == str2sites('16293G 16311C')[1].value
    
def test_commas_spaces_and_semicolons():
    assert 'G' == str2sites('16293G, 16311C; 16129')[0].value
    assert 'C' == str2sites('16293G, 16311C; 16129')[1].value
    
def test_out_of_order_insertions():
    assert len(str2sites('1.2A 1.1A')) == 2

def test_10_insertions():
    s = '1.1A 1.2A 1.3A 1.4A 1.5A 1.6A 1.7A 1.8A 1.9A 1.10A'
    assert len(str2sites(s)) == 10

