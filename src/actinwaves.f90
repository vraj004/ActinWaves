!> \file
!> $Id: DiffusionExample.f90 1528 2010-09-21 01:32:29Z chrispbradley $
!> \author Chris Bradley
!> \brief This is an example program to solve a diffusion equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example ClassicalField/ReactionDiffusion/ReactionDiffusionNoSource1D/src/ReactionDiffusionNoSource1DExample.f90
!! Example program to solve a diffusion equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/ReactionDiffusion/ReactionDiffusionNoSource1D/history.html
!<

!> Main program
PROGRAM ActinWaves

  USE OpenCMISS
  USE OpenCMISS_Iron
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Program Parameters
  REAL(CMISSRP), PARAMETER ::  store_coeff = 1.0_CMISSRP
  !Boundary Markers
  INTEGER(CMISSIntg), PARAMETER ::  MEMBRANE_MARKER=1
  INTEGER(CMISSIntg), PARAMETER ::  CYTOSOL_MARKER=0

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3

  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8

  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveMaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveEquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveEquationsSetFieldUserNumber=23
  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveSourceFieldUserNumber=12

  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveMaterialsFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveEquationsSetUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveEquationsSetFieldUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveSourceFieldUserNumber=16

  INTEGER(CMISSIntg), PARAMETER :: FActinFieldUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: stwoFieldUserNumber=26
  INTEGER(CMISSIntg), PARAMETER :: soneFieldUserNumber=27
  INTEGER(CMISSIntg), PARAMETER :: kzeroFieldUserNumber=28

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=21
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=22
  
  !cmfe type variables used for this program

  !mesh/domain related fields
  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_MeshElementsType) :: MeshElements
  TYPE(cmfe_NodesType) :: Nodes
  !simulation equations fields
  TYPE(cmfe_EquationsType) :: NPF_ActiveEquations, NPF_InActiveEquations
  TYPE(cmfe_EquationsSetType) :: NPF_ActiveEquationsSet, NPF_InActiveEquationsSet
  TYPE(cmfe_FieldType) :: NPF_ActiveEquationsSetField, NPF_InActiveEquationsSetField
  TYPE(cmfe_FieldType) :: NPF_ActiveField, NPF_InActiveField, FActinField,stwoField,soneField
  TYPE(cmfe_FieldType) :: kzeroField,NPF_ActiveMaterialsField, NPF_InActiveMaterialsField
  TYPE(cmfe_FieldType) :: NPF_ActiveSourceField, NPF_InActiveSourceField


  !equations solving related fields
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_SolverType) :: Solver, LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

  !CellML related cmiss fields
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField


#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Other program variables
  REAL(CMISSRP),ALLOCATABLE,DIMENSION(:,:) :: NodeCoords
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: ElemMap
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: NodeNums

  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENTS,CONDITION
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS,NODE_NUMBER
  INTEGER(CMISSIntg),DIMENSION(2) :: BCNODES
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER :: node,st,outputfreq,NUMBER_OF_NODES,NUMBER_OF_ATTRIBUTES,NUMBER_OF_COORDS, &
    & BOUNDARY_MARKER,nodedomain,NODES_PER_ELE,ELE_ATTRIBUTES,element,i,GeometricMeshComponent
  REAL(CMISSRP) :: init_NPFActive, init_NPFInActive, Dx_NPFActive,Dy_NPFActive,Dx_NPFInActive,Dy_NPFInActive
  REAL(CMISSRP) :: startT,endT,Tstep,ODE_TIME_STEP,nodex,nodey, VALUE, init_FActin,kzero,stwo,sone
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex,CellMLIndex,constantModelIndex
  INTEGER(CMISSIntg) :: Err
  LOGICAL :: EXPORT_FIELD
  CHARACTER(250) :: CELLID,NODEFILE,ELEMFILE,ActinPolymSignalModel  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif


!_________________________________________________________________________________________________
  !Problem INPUTS. PARAMETERS FROM Holmes et.al. 2012
  !MESH FILES
  open(unit=9,file='inputs.txt',status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening inputs file',st
    STOP
  ELSE
    PRINT *,'inputs file opened correctly'
    READ(9,*) !file info
    READ(9,*) !comment on next set of variables - reading in file info
    READ(9,*) NODEFILE,ELEMFILE
    READ(9,*)
    READ(9,*) ActinPolymSignalModel
    READ(9,*)
    READ(9,*) init_NPFActive,Dx_NPFActive,Dy_NPFActive
    READ(9,*)
    READ(9,*) init_NPFInActive,Dx_NPFInActive,Dy_NPFInActive
    READ(9,*)
    READ(9,*) init_FActin
    READ(9,*)
    READ(9,*) stwo,sone,kzero
    READ(9,*)
    READ(9,*) startT,endT,Tstep,ODE_TIME_STEP
    READ(9,*) 
    READ(9,*) outputfreq
  ENDIF
  CLOSE(9)
  !Write the params out to screen for double checking

  WRITE(*,*) 'Node file:',NODEFILE
  WRITE(*,*) 'Element file:',ELEMFILE
  WRITE(*,*) 'CellML Model File:', ActinPolymSignalModel
  !cell initial conditions
  WRITE(*,*) 'Initial [NPF_Active]i = ',init_NPFActive !dependent field (cytosolic NPFActive).
  WRITE(*,*) 'NPF_Active Diff Coeff in x = ',Dx_NPFActive
  WRITE(*,*) 'NPF_Active Diff Coeff in y = ',Dy_NPFActive
  WRITE(*,*) 'Initial [NPF_InActive]i = ',init_NPFInActive !dependent field (cytosolic NPFInActive).
  WRITE(*,*) 'NPF_InActive Diff Coeff in x = ',Dx_NPFInActive
  WRITE(*,*) 'NPF_InActive Diff Coeff in y = ',Dy_NPFInActive

  WRITE(*,*) 'Tstart=',startT
  WRITE(*,*) 'Tend=',endT
  WRITE(*,*) 'Tstep=',Tstep
  WRITE(*,*) 'ODE_Tstep=',ODE_TIME_STEP
  WRITE(*,*) 'Output Frequency=',outputfreq


  
  EXPORT_FIELD=.TRUE.
!_________________________________________________________________________________________________
  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  !Set all diganostic levels on for testing
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 2D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)
  

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)
  
!_________________________________________________________________________________________________
  !Start the creation of a biilinear-simplex basis

  CALL cmfe_Basis_Initialise(Basis,Err)
  CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
  !set the basis to bilinear simplex
  CALL cmfe_Basis_TypeSet(Basis,CMFE_BASIS_SIMPLEX_TYPE,Err)
  CALL cmfe_Basis_NumberOfXiSet(Basis,2,Err)
  CALL cmfe_Basis_InterpolationXiSet(Basis,(/CMFE_Basis_Linear_Simplex_Interpolation, & 
  & CMFE_Basis_Linear_Simplex_Interpolation/), Err)
  !Finish the creation of the basis
  CALL cmfe_Basis_CreateFinish(Basis,Err)

  !Time to create a mesh
  !Read in nodes.
  open(unit=10,file=NODEFILE,status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening node file',st
    STOP
  ELSE
    PRINT *,'Node file opened correctly'
    READ(10,*) NUMBER_OF_NODES, NUMBER_OF_COORDS, NUMBER_OF_ATTRIBUTES, BOUNDARY_MARKER
    ALLOCATE(NodeNums(NUMBER_OF_NODES,2))
    ALLOCATE(NodeCoords(NUMBER_OF_NODES,NUMBER_OF_COORDS))
    DO i = 1,NUMBER_OF_NODES
      READ(10,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2),NodeNums(i,2)
    ENDDO
  ENDIF
  CLOSE(10)
  PRINT *, 'Total Nodes',NUMBER_OF_NODES
  !Read in elements
  OPEN(unit=11,file=ELEMFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening element file',st
    STOP
  ELSE
    PRINT *,'Element file opened successfully'
    READ(11,*) NUMBER_OF_ELEMENTS,NODES_PER_ELE,ELE_ATTRIBUTES
    ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,4))
    DO i = 1,NUMBER_OF_ELEMENTS
      READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4)
    ENDDO
  ENDIF 
  CLOSE(11)
  PRINT *,'Total Elements',NUMBER_OF_ELEMENTS
  CALL MPI_BCAST(NUMBER_OF_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_COORDS,Mesh,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,1,Err)
  
  CALL cmfe_MeshElements_Initialise(MeshElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,1,Basis,MeshElements,Err)
  DO i = 1,NUMBER_OF_ELEMENTS
    element = ElemMap(i,1)
    CALL cmfe_MeshElements_NodesSet(MeshElements,element,(/ElemMap(i,2),ElemMap(i,3), &
     &   ElemMap(i,4)/),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(MeshElements,Err)
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,cmfe_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components. We have 2 field components in 1 mesh component
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,2,1,Err)

  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Set the geometric field values

  DO i = 1,NUMBER_OF_NODES
    node = NodeNums(i,1)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      nodex = NodeCoords(i,1)
      nodey = NodeCoords(i,2)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,1,nodex,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,2,nodey,Err)
     ENDIF
    ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"Cell_Geom","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"Cell_Geom","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF 
!______________________________________________________________________________________________________________

  !Create the cellml reaction with split reaction diffusion equations_set - 1 for each species
!###################
  !NPF_Active equations
!###################
  CALL cmfe_EquationsSet_Initialise(NPF_ActiveEquationsSet,Err)
  CALL cmfe_Field_Initialise(NPF_ActiveEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(NPF_ActiveEquationsSetUserNumber,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & NPF_ActiveEquationsSetFieldUserNumber,NPF_ActiveEquationsSetField,NPF_ActiveEquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(NPF_ActiveEquationsSet,Err)


  !Create the equations set dependent field variables for NPF_Active
  CALL cmfe_Field_Initialise(NPF_ActiveField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(NPF_ActiveEquationsSet, &
    & NPF_ActiveFieldUserNumber,NPF_ActiveField,Err)
  CALL cmfe_Field_VariableLabelSet(NPF_ActiveField,CMFE_FIELD_U_VARIABLE_TYPE,"NPF_Active Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(NPF_ActiveEquationsSet,Err)
  !Initialise NPF_Active dependent field
  CALL cmfe_Field_ComponentValuesInitialise(NPF_ActiveField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_NPFActive,Err)

  !Create the equations set material field variables - NPF_Active
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL cmfe_Field_Initialise(NPF_ActiveMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(NPF_ActiveEquationsSet, &
    & NPF_ActiveMaterialsFieldUserNumber,NPF_ActiveMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(NPF_ActiveMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"NPF_Active Materials Field",Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(NPF_ActiveEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(NPF_ActiveMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,Dx_NPFActive,Err) !ca diff coeff in x
  CALL cmfe_Field_ComponentValuesInitialise(NPF_ActiveMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,Dy_NPFActive,Err) !ca diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(NPF_ActiveMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,store_coeff,Err) ! storage coefficient

 
  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  CALL cmfe_Field_Initialise(NPF_ActiveSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(NPF_ActiveEquationsSet, &
    & NPF_ActiveSourceFieldUserNumber,NPF_ActiveSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(NPF_ActiveSourceField,CMFE_FIELD_U_VARIABLE_TYPE,"NPFActive Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(NPF_ActiveEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. 
  CALL cmfe_Field_ComponentValuesInitialise(NPF_ActiveSourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSRP,Err)


!###################
  !NPF_InActive equations
!###################
  CALL cmfe_EquationsSet_Initialise(NPF_InActiveEquationsSet,Err)
  CALL cmfe_Field_Initialise(NPF_InActiveEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(NPF_InActiveEquationsSetUserNumber,Region,GeometricField, &
    & [CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS, &
    & CMFE_EQUATIONS_SET_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_EQUATIONS_SET_CELLML_REAC_SPLIT_REAC_DIFF_SUBTYPE], &
    & NPF_InActiveEquationsSetFieldUserNumber,NPF_InActiveEquationsSetField,NPF_InActiveEquationsSet,Err)
  !Finish creating the equations set
  CALL cmfe_EquationsSet_CreateFinish(NPF_InActiveEquationsSet,Err)


  !Create the equations set dependent field variables for NPF_InActive
  CALL cmfe_Field_Initialise(NPF_InActiveField,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(NPF_InActiveEquationsSet, &
    & NPF_InActiveFieldUserNumber,NPF_InActiveField,Err)
  CALL cmfe_Field_VariableLabelSet(NPF_InActiveField,CMFE_FIELD_U_VARIABLE_TYPE,"NPF_InActive Field",Err)
  !Finish the equations set dependent field variables
  CALL cmfe_EquationsSet_DependentCreateFinish(NPF_InActiveEquationsSet,Err)
  !Initialise NPF_InActive dependent field
  CALL cmfe_Field_ComponentValuesInitialise(NPF_InActiveField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,init_NPFInActive,Err)

  !Create the equations set material field variables - NPF_InActive
  !by default 2 comps for reac diff i.e. diff coeff in 1 direction set constant spatially = 1, and storage coeff set to 1
  CALL cmfe_Field_Initialise(NPF_InActiveMaterialsField,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(NPF_InActiveEquationsSet, &
    & NPF_InActiveMaterialsFieldUserNumber,NPF_InActiveMaterialsField,Err)
  CALL cmfe_Field_VariableLabelSet(NPF_InActiveMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,"NPF_InActive Materials Field",Err)
  !Finish the equations set materials field variables
  CALL cmfe_EquationsSet_MaterialsCreateFinish(NPF_InActiveEquationsSet,Err)
  CALL cmfe_Field_ComponentValuesInitialise(NPF_InActiveMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,Dx_NPFInActive,Err) !NPF_InActive diff coeff in x
  CALL cmfe_Field_ComponentValuesInitialise(NPF_InActiveMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,Dy_NPFInActive,Err) !NPF_InActive diff coeff in y
  CALL cmfe_Field_ComponentValuesInitialise(NPF_InActiveMaterialsField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 3,store_coeff,Err) ! storage coefficient

 
  !Set up source field for reaction diffusion equation set. Note that for the split problem subtype, the source field is not used at all.
  CALL cmfe_Field_Initialise(NPF_InActiveSourceField,Err)
  CALL cmfe_EquationsSet_SourceCreateStart(NPF_InActiveEquationsSet, &
    & NPF_InActiveSourceFieldUserNumber,NPF_InActiveSourceField,Err)
  CALL cmfe_Field_VariableLabelSet(NPF_InActiveSourceField,CMFE_FIELD_U_VARIABLE_TYPE,"NPFInActive Source Field",Err)
  !Finish the equations set source field variables
  CALL cmfe_EquationsSet_SourceCreateFinish(NPF_InActiveEquationsSet,Err)
  !Initialising the iCaField to zero everywhere. 
  CALL cmfe_Field_ComponentValuesInitialise(NPF_InActiveSourceField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,0.0_CMISSRP,Err)

!###################
  !Set up an FActin Field
!###################
  CALL cmfe_Field_Initialise(FActinField,Err)
  CALL cmfe_Field_CreateStart(FActinFieldUserNumber,Region,FActinField,Err)
  CALL cmfe_Field_TypeSet(FActinField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FActinField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FActinField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FActinField,1,Err)
  CALL cmfe_Field_VariableTypesSet(FActinField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL cmfe_Field_DataTypeSet(FActinField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DimensionSet(FActinField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FActinField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(FActinField,CMFE_FIELD_U_VARIABLE_TYPE,"FActin Field",Err)
  CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,GeometricMeshComponent,Err)
  !Default to the geometric interpolation setup
  CALL cmfe_Field_ComponentMeshComponentSet(FActinField,CMFE_FIELD_U_VARIABLE_TYPE,1,GeometricMeshComponent,Err)            
  !Specify the interpolation to be same as geometric interpolation
  CALL cmfe_Field_ComponentInterpolationSet(FActinField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_CreateFinish(FActinField,Err)  
  CALL cmfe_Field_ComponentValuesInitialise(FActinField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,init_FActin,Err)

!###################
  !Set up a s_2 Field
!###################
  CALL cmfe_Field_Initialise(stwoField,Err)
  CALL cmfe_Field_CreateStart(stwoFieldUserNumber,Region,stwoField,Err)
  CALL cmfe_Field_TypeSet(stwoField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(stwoField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(stwoField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(stwoField,1,Err)
  CALL cmfe_Field_VariableTypesSet(stwoField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL cmfe_Field_DataTypeSet(stwoField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DimensionSet(stwoField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(stwoField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(stwoField,CMFE_FIELD_U_VARIABLE_TYPE,"stwo Field",Err)
  CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,GeometricMeshComponent,Err)
  !Default to the geometric interpolation setup
  CALL cmfe_Field_ComponentMeshComponentSet(stwoField,CMFE_FIELD_U_VARIABLE_TYPE,1,GeometricMeshComponent,Err)            
  !Specify the interpolation to be same as geometric interpolation
  CALL cmfe_Field_ComponentInterpolationSet(stwoField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_CreateFinish(stwoField,Err)  
  CALL cmfe_Field_ComponentValuesInitialise(stwoField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,stwo,Err)
!###################
  !Set up a s_1 Field
!###################
  CALL cmfe_Field_Initialise(soneField,Err)
  CALL cmfe_Field_CreateStart(soneFieldUserNumber,Region,soneField,Err)
  CALL cmfe_Field_TypeSet(soneField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(soneField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(soneField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(soneField,1,Err)
  CALL cmfe_Field_VariableTypesSet(soneField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL cmfe_Field_DataTypeSet(soneField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DimensionSet(soneField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(soneField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(soneField,CMFE_FIELD_U_VARIABLE_TYPE,"sone Field",Err)
  CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,GeometricMeshComponent,Err)
  !Default to the geometric interpolation setup
  CALL cmfe_Field_ComponentMeshComponentSet(soneField,CMFE_FIELD_U_VARIABLE_TYPE,1,GeometricMeshComponent,Err)            
  !Specify the interpolation to be same as geometric interpolation
  CALL cmfe_Field_ComponentInterpolationSet(soneField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_CreateFinish(soneField,Err)  
  CALL cmfe_Field_ComponentValuesInitialise(soneField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,sone,Err)
!###################
  !Set up k_0 Field
!###################
  CALL cmfe_Field_Initialise(kzeroField,Err)
  CALL cmfe_Field_CreateStart(kzeroFieldUserNumber,Region,kzeroField,Err)
  CALL cmfe_Field_TypeSet(kzeroField,CMFE_FIELD_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(kzeroField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(kzeroField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(kzeroField,1,Err)
  CALL cmfe_Field_VariableTypesSet(kzeroField,[CMFE_FIELD_U_VARIABLE_TYPE],Err)
  CALL cmfe_Field_DataTypeSet(kzeroField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DP_TYPE,Err)
  CALL cmfe_Field_DimensionSet(kzeroField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_SCALAR_DIMENSION_TYPE,Err)
  CALL cmfe_Field_NumberOfComponentsSet(kzeroField,CMFE_FIELD_U_VARIABLE_TYPE,1,Err)
  CALL cmfe_Field_VariableLabelSet(kzeroField,CMFE_FIELD_U_VARIABLE_TYPE,"kzero Field",Err)
  CALL cmfe_Field_ComponentMeshComponentGet(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,1,GeometricMeshComponent,Err)
  !Default to the geometric interpolation setup
  CALL cmfe_Field_ComponentMeshComponentSet(kzeroField,CMFE_FIELD_U_VARIABLE_TYPE,1,GeometricMeshComponent,Err)            
  !Specify the interpolation to be same as geometric interpolation
  CALL cmfe_Field_ComponentInterpolationSet(kzeroField,CMFE_FIELD_U_VARIABLE_TYPE,1, &
    & CMFE_FIELD_NODE_BASED_INTERPOLATION,Err)
  CALL cmfe_Field_CreateFinish(kzeroField,Err)  
  CALL cmfe_Field_ComponentValuesInitialise(kzeroField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,kzero,Err)

!______________________________________________________________________________________________________________  
  !Start to set up CellML Fields

  !Create the CellML environment
  CALL cmfe_CellML_Initialise(CellML,Err)
  CALL cmfe_CellML_CreateStart(CellMLUserNumber,Region,CellML,Err)

  !Import a toy constant source model from a file
  CALL cmfe_CellML_ModelImport(CellML,ActinPolymSignalModel,constantModelIndex,Err)

  CALL cmfe_CellML_VariableSetAsKnown(CellML,constantModelIndex,"cell/s_2",Err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,constantModelIndex,"cell/s_1",Err)
  CALL cmfe_CellML_VariableSetAsKnown(CellML,constantModelIndex,"cell/k_o",Err)
  CALL cmfe_CellML_VariableSetAsWanted(CellML,constantModelIndex,"cell/h",Err)

  !Finish the CellML environment
  CALL cmfe_CellML_CreateFinish(CellML,Err)

  !Start the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateStart(CellML,Err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,stwoField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"cell/s_2",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,soneField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"cell/s_1",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,kzeroField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"cell/k_o",CMFE_FIELD_VALUES_SET_TYPE,Err)


  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,NPF_ActiveField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"cell/A_npf",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"cell/A_npf",CMFE_FIELD_VALUES_SET_TYPE, &
    & NPF_ActiveField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,NPF_InActiveField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"cell/I_npf",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"cell/I_npf",CMFE_FIELD_VALUES_SET_TYPE, &
    & NPF_InActiveField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  CALL cmfe_CellML_CreateFieldToCellMLMap(CellML,FActinField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE, &
    & constantModelIndex,"cell/F_actin",CMFE_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_CellML_CreateCellMLToFieldMap(CellML,constantModelIndex,"cell/F_actin",CMFE_FIELD_VALUES_SET_TYPE, &
    & FActinField,CMFE_FIELD_U_VARIABLE_TYPE,1,CMFE_FIELD_VALUES_SET_TYPE,Err)

  !Finish the creation of CellML <--> OpenCMISS field maps
  CALL cmfe_CellML_FieldMapsCreateFinish(CellML,Err)

  !Start the creation of the CellML models field. This field is an integer field that stores which nodes have which cellml model
  CALL cmfe_Field_Initialise(CellMLModelsField,Err)
  CALL cmfe_CellML_ModelsFieldCreateStart(CellML, CellMLModelsFieldUserNumber, &
    & CellMLModelsField,Err)
  !Finish the creation of the CellML models field
  CALL cmfe_CellML_ModelsFieldCreateFinish(CellML,Err)
  !The CellMLModelsField is an integer field that stores which model is being used by which node.
  !By default all field parameters have default model value of 1, i.e. the first model. But, this command below is for example purposes
  CALL cmfe_Field_ComponentValuesInitialise(CellMLModelsField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,1_CMISSIntg,Err)

  !Start the creation of the CellML state field
  CALL cmfe_Field_Initialise(CellMLStateField,Err)
  CALL cmfe_CellML_StateFieldCreateStart(CellML, & 
    & CellMLStateFieldUserNumber,CellMLStateField,Err)
  !Finish the creation of the CellML state field
  CALL cmfe_CellML_StateFieldCreateFinish(CellML,Err)

  !Start the creation of the CellML intermediate field
  CALL cmfe_Field_Initialise(CellMLIntermediateField,Err)
  CALL cmfe_CellML_IntermediateFieldCreateStart(CellML, &
    & CellMLIntermediateFieldUserNumber,CellMLIntermediateField,Err)
  !Finish the creation of the CellML intermediate field
  CALL cmfe_CellML_IntermediateFieldCreateFinish(CellML,Err)

  !Start the creation of CellML parameters field
  CALL cmfe_Field_Initialise(CellMLParametersField,Err)
  CALL cmfe_CellML_ParametersFieldCreateStart(CellML, &
    & CellMLParametersFieldUserNumber,CellMLParametersField,Err)
  !Finish the creation of CellML parameters
  CALL cmfe_CellML_ParametersFieldCreateFinish(CellML,Err)
!______________________________________________________________________________________________________________  
  !Create the equations set equations for NPF_Active
  CALL cmfe_Equations_Initialise(NPF_ActiveEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(NPF_ActiveEquationsSet,NPF_ActiveEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(NPF_ActiveEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(NPF_ActiveEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(NPF_ActiveEquationsSet,Err)

  CALL cmfe_Equations_Initialise(NPF_InActiveEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(NPF_InActiveEquationsSet,NPF_InActiveEquations,Err)
  !Set the equations matrices sparsity type
  CALL cmfe_Equations_SparsityTypeSet(NPF_InActiveEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  !Set the equations set output
  CALL cmfe_Equations_OutputTypeSet(NPF_InActiveEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsTimingOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsMatrixOutput,Err)
  !CALL CMISSEquationsOutputTypeSet(Equations,CMISSEquationsElementMatrixOutput,Err)
  !Finish the equations set equations
  CALL cmfe_EquationsSet_EquationsCreateFinish(NPF_InActiveEquationsSet,Err)
!______________________________________________________________________________________________________________  
  !Create the problem
  WRITE(*,*) 'Create the problem'

  CALL cmfe_Problem_Initialise(Problem,Err)
  CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_CLASSICAL_FIELD_CLASS, &
    & CMFE_PROBLEM_REACTION_DIFFUSION_EQUATION_TYPE,CMFE_PROBLEM_CELLML_REAC_INTEG_REAC_DIFF_STRANG_SPLIT_SUBTYPE],Problem,Err)
  !Finish the creation of a problem.
  CALL cmfe_Problem_CreateFinish(Problem,Err)

  !Create the problem control
  CALL cmfe_Problem_ControlLoopCreateStart(Problem,Err)
  !Get the control loop
  CALL cmfe_ControlLoop_Initialise(ControlLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(Problem,CMFE_CONTROL_LOOP_NODE,ControlLoop,Err)
  !Set the times
  CALL cmfe_ControlLoop_TimesSet(ControlLoop,startT,endT,Tstep,Err)
  CALL cmfe_ControlLoop_TimeOutputSet(ControlLoop,outputfreq,Err)
  CALL cmfe_ControlLoop_OutputTypeSet(ControlLoop,CMFE_CONTROL_LOOP_PROGRESS_OUTPUT,Err)

  !Finish creating the problem control loop
  CALL cmfe_Problem_ControlLoopCreateFinish(Problem,Err)
!______________________________________________________________________________________________________________
  !CALL cmfe_DiagnosticsSetOn(cmfe_AllDiagType,[1,2,3,4,5],"Diagnostics",["cmfe_SolverTypeInitialise"],Err)
  WRITE(*,*) 'Set up the problem solvers for Strang splitting'
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)
  !First solver is a DAE solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)

  !Second solver is the dynamic solver for solving the parabolic equation
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Solver_Initialise(LinearSolver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  !set theta - backward vs forward time step parameter
  CALL cmfe_Solver_DynamicThetaSet(Solver,1.0_CMISSRP,Err)
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)
  !get the dynamic linear solver from the solver
  CALL cmfe_Solver_DynamicLinearSolverGet(Solver,LinearSolver,Err)
  !set linear solver to be direct solver. Note, I found this stuff in fluidmechanics/darcy/dynamic/src example
  !CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  !CALL CMISSSolverLibraryTypeSet(LinearSolver,CMISSSolverCMISSLibrary,Err)
  !CALL CMISSSolverLinearTypeSet(LinearSolver,CMISSSolverLinearDirectSolveType,Err)
  !CALL CMISSSolverLibraryTypeSet(LinearSolver,CMISSSolverMUMPSLibrary,Err)
  CALL cmfe_Solver_LinearIterativeMaximumIterationsSet(LinearSolver,1000,Err)


  !Third solver is another DAE solver
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL cmfe_Solver_DAETimeStepSet(Solver,ODE_TIME_STEP,Err) !set the third solver's integration time step
  CALL cmfe_Solver_OutputTypeSet(Solver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverTimingOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISSSolverSolverOutput,Err)
  !CALL CMISSSolverOutputTypeSet(Solver,CMISS_SOLVER_PROGRESS_OUTPUT,Err)

  !Finish the creation of the problem solver
  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)

!_______________________________________________________________________________________________________________
  !Start the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)
  !Get the first solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,1,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Get the third solver  
  !Get the CellML equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,3,Solver,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  CALL cmfe_Solver_CellMLEquationsGet(Solver,CellMLEquations,Err)
  !Add in the CellML environement
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellML,CellMLIndex,Err)

  !Finish the creation of the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

!_______________________________________________________________________________________________________________
  !Start the creation of the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)
  !Get the second solver  
  !Get the solver equations
  CALL cmfe_Solver_Initialise(Solver,Err)
  CALL cmfe_Problem_SolverGet(Problem,CMFE_CONTROL_LOOP_NODE,2,Solver,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquations,Err)
  CALL cmfe_Solver_SolverEquationsGet(Solver,SolverEquations,Err)
  !Set the solver equations sparsity
  !CALL CMISSSolverEquationsSparsityTypeSet(SolverEquations,CMISSSolverEquationsSparseMatrices,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquations,CMFE_SOLVER_SPARSE_MATRICES,Err)  
  !Add in the equations set for NPF_Active, NPF_InActive
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,NPF_ActiveEquationsSet,EquationsSetIndex,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquations,NPF_InActiveEquationsSet,EquationsSetIndex,Err)

  !Finish the creation of the problem solver equations
  CALL Cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

!_________________________________________________________________________________________________________
  !Create the solver equations set boundary conditions
  WRITE(*,*) 'Set up boundary conditions' 
  CALL cmfe_BoundaryConditions_Initialise(BoundaryConditions,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(SolverEquations,BoundaryConditions,Err)
 !Set no flux on cell boundary
  DO node=1,NUMBER_OF_NODES
    IF(NodeNums(node,2).EQ.1) THEN
      NODE_NUMBER = NodeNums(node,1)
      CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
      IF(NodeDomain==ComputationalNodeNumber) THEN
        CONDITION = CMFE_BOUNDARY_CONDITION_FIXED
        VALUE=0.0_CMISSDP
        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,NPF_ActiveField, &
          & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)

        CALL cmfe_BoundaryConditions_SetNode(BoundaryConditions,NPF_InActiveField, &
          & CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,CMFE_NO_GLOBAL_DERIV, &
          & NODE_NUMBER,1,CONDITION,VALUE,Err) !(neumann boundary condition - no flux)
      ENDIF
    ENDIF
  ENDDO

  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(SolverEquations,Err)
!__________________________________________________________________________________________________________
  !Solve the problem
  CALL cmfe_Problem_Solve(Problem,Err)
!__________________________________________________________________________________________________________
  EXPORT_FIELD=.TRUE.
  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"ActinWaves","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"ActinWaves","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)

  ENDIF
  

  WRITE(*,'(A)') "Program successfully completed."
  

  STOP
  
END PROGRAM ActinWaves
