//This file is autogenerated by the python_src/extract_test_cases.py script, do not modify this manually, as it will be lost.
#include "unity_fixture.h"

TEST_GROUP_RUNNER(string_handling)
{
    RUN_TEST_CASE(string_handling,simple_append_num_to_name);
    RUN_TEST_CASE(string_handling,append_num_to_complex_name);
    RUN_TEST_CASE(string_handling,append_longer_num);
}


TEST_GROUP_RUNNER(file_creation_w1_4)
{
    RUN_TEST_CASE(file_creation_w1_4,create_not_existing);
    RUN_TEST_CASE(file_creation_w1_4,fail_to_create_existing_if_mode_is_fail_if_exists);
    RUN_TEST_CASE(file_creation_w1_4,override_existing_if_mode_says_so);
    RUN_TEST_CASE(file_creation_w1_4,ignore_existing_and_do_nothing_if_mode_is_no_op_on_existing);
    RUN_TEST_CASE(file_creation_w1_4,can_create_with_noop_if_exists);
    RUN_TEST_CASE(file_creation_w1_4,can_create_with_noop_if_not_exists);
    RUN_TEST_CASE(file_creation_w1_4,cannot_create_with_readonly);
    RUN_TEST_CASE(file_creation_w1_4,has_proper_group_structure);
    RUN_TEST_CASE(file_creation_w1_4,evade_if_exists_simply_creates_file_if_non_existant);
    RUN_TEST_CASE(file_creation_w1_4,evade_if_exists_evades_one);
}


TEST_GROUP_RUNNER(read_and_write_single_w1_3_6)
{
    RUN_TEST_CASE(read_and_write_single_w1_3_6,write_native_int);
    RUN_TEST_CASE(read_and_write_single_w1_3_6,write_and_read_native_int);
    RUN_TEST_CASE(read_and_write_single_w1_3_6,do_not_write_to_noop_file);
    RUN_TEST_CASE(read_and_write_single_w1_3_6,fail_to_read_with_wrong_sized_communicator);
    RUN_TEST_CASE(read_and_write_single_w1_3_6,write_and_read_native_double);
    RUN_TEST_CASE(read_and_write_single_w1_3_6,write_and_read_various_types);
    RUN_TEST_CASE(read_and_write_single_w1_3_6,rw_weird_characters);
}


TEST_GROUP_RUNNER(read_write_array_w1_4)
{
    RUN_TEST_CASE(read_write_array_w1_4,read_int);
    RUN_TEST_CASE(read_write_array_w1_4,read_uint16);
    RUN_TEST_CASE(read_write_array_w1_4,read_double);
    RUN_TEST_CASE(read_write_array_w1_4,read_float);
    RUN_TEST_CASE(read_write_array_w1_4,read_uint64);
    RUN_TEST_CASE(read_write_array_w1_4,write_unequal_amounts);
    RUN_TEST_CASE(read_write_array_w1_4,write_only_rank_0);
}


TEST_GROUP_RUNNER(convenience_api_w1_3)
{
    RUN_TEST_CASE(convenience_api_w1_3,sample_use);
}


TEST_GROUP_RUNNER(communicate_density_fields_w2)
{
    RUN_TEST_CASE(communicate_density_fields_w2,2domains_1rankperdomain_3types);
}


TEST_GROUP_RUNNER(communicate_density_fields_w8)
{
    RUN_TEST_CASE(communicate_density_fields_w8,2domains_4ranksperdomain_2types);
    RUN_TEST_CASE(communicate_density_fields_w8,8domains);
    RUN_TEST_CASE(communicate_density_fields_w8,4domains_2rankseach_typesaredifferent);
    RUN_TEST_CASE(communicate_density_fields_w8,ranks_in_domain_are_different);
}


TEST_GROUP_RUNNER(calc_non_bonded_energy_w2)
{
    RUN_TEST_CASE(calc_non_bonded_energy_w2,example);
}


