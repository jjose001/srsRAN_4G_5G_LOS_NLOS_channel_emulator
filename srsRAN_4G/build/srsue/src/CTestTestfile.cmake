# CMake generated Testfile for 
# Source directory: /home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/src
# Build directory: /home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/build/srsue/src
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ue_rf_failure "srsue" "/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/ue.conf.example" "--rf.device_name=zmq")
set_tests_properties(ue_rf_failure PROPERTIES  _BACKTRACE_TRIPLES "/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/src/CMakeLists.txt;64;add_test;/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/src/CMakeLists.txt;0;")
add_test(ue_rf_failure_max_channels "srsue" "/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/ue.conf.example" "--rf.device_name=zmq" "--rf.nof_antennas=4" "--rat.eutra.nof_carriers=5")
set_tests_properties(ue_rf_failure_max_channels PROPERTIES  _BACKTRACE_TRIPLES "/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/src/CMakeLists.txt;65;add_test;/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/src/CMakeLists.txt;0;")
add_test(ue_rf_failure_exceeds_channels "srsue" "/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/ue.conf.example" "--rf.device_name=zmq" "--rf.nof_antennas=5" "--rat.eutra.nof_carriers=5")
set_tests_properties(ue_rf_failure_exceeds_channels PROPERTIES  _BACKTRACE_TRIPLES "/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/src/CMakeLists.txt;66;add_test;/home/joao/Escritorio/GITHUB/srsRAN_INVESTIGACION_JOAO/srsRAN_4G/srsue/src/CMakeLists.txt;0;")
subdirs("phy")
subdirs("stack")
subdirs("test")
