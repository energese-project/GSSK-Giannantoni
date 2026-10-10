#!/usr/bin/env python3
"""Regenerates the T-ROB-01 corpus (NFR-ROB-001): python3 tests/mop_fuzz/generate.py

ok_* cases load (ok_refused_* load, and the kernel then refuses them: T-OUT-03);
every other case is a load error. The files are committed; this script only
records how they were made."""
import os
D=os.path.dirname(os.path.abspath(__file__))
base_nodes='[{"id":"src","type":"source","initial_value":1.0},{"id":"a","type":"storage","current_level":1.0},{"id":"b","type":"storage","current_level":1.0},{"id":"c","type":"storage","current_level":1.0},{"id":"out","type":"sink"}]'
base_edges='[{"source":"src","target":"a","weight":0.1},{"source":"a","target":"b","weight":0.1},{"source":"b","target":"c","weight":0.1},{"source":"c","target":"out","weight":0.1}]'
def seed(mop):
    m = '' if mop is None else ',"mop":'+mop
    return '{"system_name":"fuzz","nodes":'+base_nodes+',"edges":'+base_edges+',"simulation_params":{"t_val":1.5,"derivative_order":2,"generative_mode":false}'+m+'}\n'
aff='{"from":"a","to":"b","a":1,"b":0.25,"p":1}'
many_samples='['+','.join('[%d,%g,%g]'%(i,1+0.01*i,0.001*i) for i in range(1000))+']'
deep='['*200+']'*200
long_id='x'*5000
cases = [
 ('ok_none', None),
 ('ok_empty_beta', '{"k":1,"beta":[]}'),
 ('ok_affine', '{"k":1,"beta":['+aff+']}'),
 ('ok_network', '{"k":3,"beta":"network"}'),
 ('ok_fraction', '{"k":{"num":3,"den":6},"beta":['+aff+']}'),
 ('ok_all_six_couples', '{"k":2,"beta":['+','.join('{"from":"%s","to":"%s","a":[1,0.1],"b":[0.2,0],"p":0.5}'%(f,t) for f in 'abc' for t in 'abc' if f!=t)+']}'),
 ('ok_1000_samples', '{"k":1,"beta":[{"from":"a","to":"c","samples":'+many_samples+'}]}'),
 ('ok_negative_zero', '{"k":1,"beta":[{"from":"a","to":"b","a":1,"b":-0.0,"p":1}]}'),
 ('ok_full', '{"k":1,"reference":["b","c"],"beta":['+aff+'],"second_equation":{"alpha12_0":[0.3,0.1],"c1":1,"c2":0.5},"eqs":{"psi1":[1,1,1],"psi2":1,"epsilon":[0,0,0],"A":1}}'),
 ('ok_refused_domain', '{"k":1,"beta":[{"from":"a","to":"b","a":1,"b":-1,"p":1}]}'),
 ('ok_refused_second', '{"k":1,"beta":[],"second_equation":{"alpha12_0":0.3,"c1":1,"c2":-1}}'),
 ('ok_refused_eqs_x11', '{"k":1,"beta":[],"eqs":{"psi1":[1,1,1],"psi2":1,"epsilon":[0,0.1,0.2],"A":1}}'),
 ('bad_mop_null', 'null'),
 ('bad_mop_true', 'true'),
 ('bad_mop_string', '"mop"'),
 ('bad_mop_empty', '{}'),
 ('bad_duplicate_key_k', '{"k":1,"k":2,"beta":[]}'),
 ('bad_duplicate_key_in_couple', '{"k":1,"beta":[{"from":"a","from":"c","to":"b","a":1,"b":0,"p":1}]}'),
 ('bad_k_negative', '{"k":-1,"beta":[]}'),
 ('bad_k_huge', '{"k":4294967296,"beta":[]}'),
 ('bad_k_den_negative', '{"k":{"num":1,"den":-2},"beta":[]}'),
 ('bad_k_num_fraction', '{"k":{"num":1.5,"den":2},"beta":[]}'),
 ('bad_k_array', '{"k":[1],"beta":[]}'),
 ('bad_k_null', '{"k":null,"beta":[]}'),
 ('bad_beta_deep_nesting', '{"k":1,"beta":'+deep+'}'),
 ('bad_beta_null', '{"k":1,"beta":null}'),
 ('bad_couple_empty_ids', '{"k":1,"beta":[{"from":"","to":"","a":1,"b":0,"p":1}]}'),
 ('bad_couple_long_id', '{"k":1,"beta":[{"from":"'+long_id+'","to":"b","a":1,"b":0,"p":1}]}'),
 ('bad_couple_unicode_id', '{"k":1,"beta":[{"from":"\\u00e1","to":"b","a":1,"b":0,"p":1}]}'),
 ('bad_couple_module_free_sink', '{"k":1,"beta":[{"from":"a","to":"out","a":1,"b":0,"p":1}]}'),
 ('bad_a_three_elements', '{"k":1,"beta":[{"from":"a","to":"b","a":[1,2,3],"b":0,"p":1}]}'),
 ('bad_a_nested', '{"k":1,"beta":[{"from":"a","to":"b","a":[[1],0],"b":0,"p":1}]}'),
 ('bad_p_overflow', '{"k":1,"beta":[{"from":"a","to":"b","a":1,"b":0,"p":1e999}]}'),
 ('bad_b_minus_overflow', '{"k":1,"beta":[{"from":"a","to":"b","a":1,"b":-1e999,"p":1}]}'),
 ('bad_samples_empty', '{"k":1,"beta":[{"from":"a","to":"b","samples":[]}]}'),
 ('bad_samples_object', '{"k":1,"beta":[{"from":"a","to":"b","samples":{"t":0}}]}'),
 ('bad_samples_strings', '{"k":1,"beta":[{"from":"a","to":"b","samples":[["0","1","0"],["1","1","0"]]}]}'),
 ('bad_samples_four_wide', '{"k":1,"beta":[{"from":"a","to":"b","samples":[[0,1,0,9],[1,1,0,9]]}]}'),
 ('bad_samples_decreasing', '{"k":1,"beta":[{"from":"a","to":"b","samples":[[0,1,0],[2,1,0],[1,1,0]]}]}'),
 ('bad_samples_inf_value', '{"k":1,"beta":[{"from":"a","to":"b","samples":[[0,1,0],[1,1e999,0]]}]}'),
 ('bad_reference_numbers', '{"k":1,"reference":[1,2],"beta":[]}'),
 ('bad_reference_three', '{"k":1,"reference":["a","b","c"],"beta":[]}'),
 ('bad_second_string', '{"k":1,"beta":[],"second_equation":"yes"}'),
 ('bad_second_alpha_three', '{"k":1,"beta":[],"second_equation":{"alpha12_0":[1,2,3],"c1":1,"c2":0}}'),
 ('bad_eqs_epsilon_four', '{"k":1,"beta":[],"eqs":{"psi1":[1,1,1],"psi2":1,"epsilon":[0,0,0,0],"A":1}}'),
 ('bad_eqs_psi2_array', '{"k":1,"beta":[],"eqs":{"psi1":[1,1,1],"psi2":[1],"epsilon":[0,0,0],"A":1}}'),
 ('bad_truncated_json', '{"k":1,"beta":[{"from":"a"'),
]
names=[]
for name, mop in cases:
    fn=name+'.json'
    txt = seed(mop)
    if name=='bad_truncated_json': txt = txt[:-3]
    open(os.path.join(D,fn),'w').write(txt)
    names.append(fn)
open(os.path.join(D,'INDEX'),'w').write('# T-ROB-01 corpus (NFR-ROB-001). ok_* must load; every other case must be a load error.\n'+'\n'.join(names)+'\n')
print(len(names))
