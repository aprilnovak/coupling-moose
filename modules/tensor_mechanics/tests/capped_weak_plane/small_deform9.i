# apply a shear deformation to observe shear hardening.
# Shear gives q_trial = 2*Sqrt(20), p_trial=0
# The solution given by MOOSE correctly satisfies the equations
# 0 = f = q + p*tan(phi) - coh
# 0 = p - p_trial + ga * Ezzzz * dg/dp
# 0 = q - q_trial + ga * Ezxzx * dg/dq
# Here dg/dp = tan(psi), and dg/dq = 1
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -0.5
  zmax = 0.5
[]


[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
[]


[BCs]
  [./bottomx]
    type = PresetBC
    variable = disp_x
    boundary = back
    value = 0.0
  [../]
  [./bottomy]
    type = PresetBC
    variable = disp_y
    boundary = back
    value = 0.0
  [../]
  [./bottomz]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]

  [./topx]
    type = FunctionPresetBC
    variable = disp_x
    boundary = front
    function = 't'
  [../]
  [./topy]
    type = FunctionPresetBC
    variable = disp_y
    boundary = front
    function = '2*t'
  [../]
  [./topz]
    type = FunctionPresetBC
    variable = disp_z
    boundary = front
    function = '0'
  [../]
[]


[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strainp_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strainp_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strainp_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strainp_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strainp_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./strainp_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./straint_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./straint_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./straint_xz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./straint_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./straint_yz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./straint_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_shear]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_tensile]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_compressive]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./intnl_shear]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./intnl_tensile]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./iter]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ls]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  [../]
  [./stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  [../]
  [./stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
  [./strainp_xx]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = strainp_xx
    index_i = 0
    index_j = 0
  [../]
  [./strainp_xy]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = strainp_xy
    index_i = 0
    index_j = 1
  [../]
  [./strainp_xz]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = strainp_xz
    index_i = 0
    index_j = 2
  [../]
  [./strainp_yy]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = strainp_yy
    index_i = 1
    index_j = 1
  [../]
  [./strainp_yz]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = strainp_yz
    index_i = 1
    index_j = 2
  [../]
  [./strainp_zz]
    type = RankTwoAux
    rank_two_tensor = plastic_strain
    variable = strainp_zz
    index_i = 2
    index_j = 2
  [../]
  [./straint_xx]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = straint_xx
    index_i = 0
    index_j = 0
  [../]
  [./straint_xy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = straint_xy
    index_i = 0
    index_j = 1
  [../]
  [./straint_xz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = straint_xz
    index_i = 0
    index_j = 2
  [../]
  [./straint_yy]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = straint_yy
    index_i = 1
    index_j = 1
  [../]
  [./straint_yz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = straint_yz
    index_i = 1
    index_j = 2
  [../]
  [./straint_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = straint_zz
    index_i = 2
    index_j = 2
  [../]
  [./f_shear]
    type = MaterialStdVectorAux
    property = plastic_yield_function
    index = 0
    variable = f_shear
  [../]
  [./f_tensile]
    type = MaterialStdVectorAux
    property = plastic_yield_function
    index = 1
    variable = f_tensile
  [../]
  [./f_compressive]
    type = MaterialStdVectorAux
    property = plastic_yield_function
    index = 2
    variable = f_compressive
  [../]
  [./intnl_shear]
    type = MaterialStdVectorAux
    property = plastic_internal_parameter
    index = 0
    variable = intnl_shear
  [../]
  [./intnl_tensile]
    type = MaterialStdVectorAux
    property = plastic_internal_parameter
    index = 1
    variable = intnl_tensile
  [../]
  [./iter]
    type = MaterialRealAux
    property = plastic_NR_iterations
    variable = iter
  [../]
  [./ls]
    type = MaterialRealAux
    property = plastic_linesearch_needed
    variable = ls
  [../]
[]

[Postprocessors]
  [./stress_xx]
    type = PointValue
    point = '0 0 0'
    variable = stress_xx
  [../]
  [./stress_xy]
    type = PointValue
    point = '0 0 0'
    variable = stress_xy
  [../]
  [./stress_xz]
    type = PointValue
    point = '0 0 0'
    variable = stress_xz
  [../]
  [./stress_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
  [./stress_yz]
    type = PointValue
    point = '0 0 0'
    variable = stress_yz
  [../]
  [./stress_zz]
    type = PointValue
    point = '0 0 0'
    variable = stress_zz
  [../]
  [./strainp_xx]
    type = PointValue
    point = '0 0 0'
    variable = strainp_xx
  [../]
  [./strainp_xy]
    type = PointValue
    point = '0 0 0'
    variable = strainp_xy
  [../]
  [./strainp_xz]
    type = PointValue
    point = '0 0 0'
    variable = strainp_xz
  [../]
  [./strainp_yy]
    type = PointValue
    point = '0 0 0'
    variable = strainp_yy
  [../]
  [./strainp_yz]
    type = PointValue
    point = '0 0 0'
    variable = strainp_yz
  [../]
  [./strainp_zz]
    type = PointValue
    point = '0 0 0'
    variable = strainp_zz
  [../]
  [./straint_xx]
    type = PointValue
    point = '0 0 0'
    variable = straint_xx
  [../]
  [./straint_xy]
    type = PointValue
    point = '0 0 0'
    variable = straint_xy
  [../]
  [./straint_xz]
    type = PointValue
    point = '0 0 0'
    variable = straint_xz
  [../]
  [./straint_yy]
    type = PointValue
    point = '0 0 0'
    variable = straint_yy
  [../]
  [./straint_yz]
    type = PointValue
    point = '0 0 0'
    variable = straint_yz
  [../]
  [./straint_zz]
    type = PointValue
    point = '0 0 0'
    variable = straint_zz
  [../]
  [./f_shear]
    type = PointValue
    point = '0 0 0'
    variable = f_shear
  [../]
  [./f_tensile]
    type = PointValue
    point = '0 0 0'
    variable = f_tensile
  [../]
  [./f_compressive]
    type = PointValue
    point = '0 0 0'
    variable = f_compressive
  [../]
  [./intnl_shear]
    type = PointValue
    point = '0 0 0'
    variable = intnl_shear
  [../]
  [./intnl_tensile]
    type = PointValue
    point = '0 0 0'
    variable = intnl_tensile
  [../]
  [./iter]
    type = PointValue
    point = '0 0 0'
    variable = iter
  [../]
  [./ls]
    type = PointValue
    point = '0 0 0'
    variable = ls
  [../]
[]

[UserObjects]
  [./coh]
    type = TensorMechanicsHardeningExponential
    value_0 = 1
    value_residual = 2
    rate = 1
  [../]
  [./tanphi]
    type = TensorMechanicsHardeningExponential
    value_0 = 1.0
    value_residual = 0.5
    rate = 2
  [../]
  [./tanpsi]
    type = TensorMechanicsHardeningExponential
    value_0 = 0.1
    value_residual = 0.05
    rate = 1
  [../]
  [./t_strength]
    type = TensorMechanicsHardeningConstant
    value = 1E8
  [../]
  [./c_strength]
    type = TensorMechanicsHardeningConstant
    value = 1E8
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric_isotropic
    C_ijkl = '4 4'
  [../]
  [./strain]
    type = ComputeIncrementalSmallStrain
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./admissible]
    type = ComputeMultipleInelasticStress
    inelastic_models = stress
    perform_finite_strain_rotations = false
  [../]
  [./stress]
    type = CappedWeakPlaneStressUpdate
    cohesion = coh
    tan_friction_angle = tanphi
    tan_dilation_angle = tanpsi
    tensile_strength = t_strength
    compressive_strength = c_strength
    max_NR_iterations = 20
    tip_smoother = 0
    smoothing_tol = 1
    yield_function_tol = 1E-3
  [../]
[]


[Executioner]
  end_time = 1
  dt = 1
  type = Transient
[]


[Outputs]
  file_base = small_deform9
  [./csv]
    type = CSV
  [../]
[]
