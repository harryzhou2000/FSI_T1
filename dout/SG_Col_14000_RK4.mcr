#!MC 1410
$!VarSet |NumLoop| =221
$!Loop |NumLoop|      
$!Varset |num| = ( |Loop| *100 -100 + 14000)
$!ReadDataSet  '"E:\Documents\FSI_T1\Unstructured\UGSolver_Submit\dout\SG_Col_T1_SHOCK90_RKTestdata_0_AT|num|.dat" '
  ReadDataOption = New
  ResetStyle = No
  VarLoadMode = ByName
  AssignStrandIDs = Yes
  VarNameList = '"X" "Y" "rho" "U" "V" "E"'
$!ExportSetup ExportFName = 'E:\Documents\FSI_T1\Unstructured\UGSolver_Submit\dout\SG_Col_SHOCK_VOR_RK_|num|.png'
$!Export 
  ExportRegion = AllFrames
$!EndLoop