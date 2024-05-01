#!/usr/bin/env bats

source $QP_ROOT/tests/bats/common.bats.sh
source $QP_ROOT/quantum_package.rc


function run_Ne() {
  rm -rf Ne_tc_scf
  echo Ne > Ne.xyz
  qp create_ezfio -b cc-pcvdz Ne.xyz -o Ne_tc_scf
  qp run scf

  qp set tc_keywords tc_integ_type numeric
  qp set jastrow env_type Sum_Gauss
  qp set hamiltonian mu_erf 0.87
  qp set jastrow j1e_type None
  qp set jastrow env_coef "[1.]"
  qp set jastrow env_expo "[1.5]"

  qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out
  eref=-128.552134
  energy="$(qp get tc_scf bitc_energy)"
  eq $energy $eref 2e-4
}


@test "Ne" {
 run_Ne
}

function run_C() {
  rm -rf C_tc_scf
  echo C  > C.xyz
  qp create_ezfio -b cc-pcvdz C.xyz -o C_tc_scf -m 3
  qp run scf

  qp set tc_keywords tc_integ_type numeric
  qp set jastrow env_type Sum_Gauss
  qp set hamiltonian mu_erf 0.87
  qp set jastrow j1e_type None
  qp set jastrow env_coef "[1.]"
  qp set jastrow env_expo "[1.5]"

  qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out
  eref=-37.691254356408791
  energy="$(qp get tc_scf bitc_energy)"
  eq $energy $eref 2e-4
}


@test "C" {
 run_C
}


function run_O() {
  rm -rf O_tc_scf
  echo O  > O.xyz
  qp create_ezfio -b cc-pcvdz O.xyz -o O_tc_scf -m 3
  qp run scf

  qp set tc_keywords tc_integ_type numeric
  qp set jastrow env_type Sum_Gauss
  qp set jastrow j1e_type None
  qp set jastrow env_coef "[1.]"
  qp set jastrow env_expo "[1.5]"
  qp set hamiltonian mu_erf 0.87

  qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out
  eref=-74.814687229354590
  energy="$(qp get tc_scf bitc_energy)"
  eq $energy $eref 2e-4
}


@test "O" {
 run_O
}



function run_ch2() {
  rm -rf ch2_tc_scf
  cp ${QP_ROOT}/tests/input/ch2.xyz .
  qp create_ezfio -b "C:cc-pcvdz|H:cc-pvdz" ch2.xyz -o ch2_tc_scf
  qp run scf

  qp set tc_keywords tc_integ_type numeric
  qp set jastrow env_type Sum_Gauss
  qp set jastrow j1e_type None
  qp set jastrow env_coef "[1., 1., 1.]"
  qp set jastrow env_expo '[1.5,10000,10000]'
  qp set hamiltonian mu_erf 0.87

  qp run tc_scf | tee ${EZFIO_FILE}.tc_scf.out
  eref=-38.903247818077737
  energy="$(qp get tc_scf bitc_energy)"
  eq $energy $eref 2e-4
}


@test "ch2" {
 run_ch2
}

