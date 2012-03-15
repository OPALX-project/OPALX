#!/bin/bash

## Diagnostics
rm -f Diagnostics/OPAL/*.T7
rm -f Diagnostics/OPAL/TwoStep.h5
ln -s ${PWD}/Diagnostics/TwoStep.h5 Diagnostics/OPAL/TwoStep.h5
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB02-MSLAC.T7 Diagnostics/OPAL/FINLB02-MSLAC.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FIND1-MB01.T7 Diagnostics/OPAL/FIND1-MB01.T7

rm -f Diagnostics/ImpactT/*.T7
rm -f Diagnostics/ImpactT/rfdata*
rm -f Diagnostics/ImpactT/TwoStep.h5
ln -s ${PWD}/Diagnostics/TwoStep.h5 Diagnostics/ImpactT/TwoStep.h5
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-MSLAC.T7 Diagnostics/ImpactT/1T6.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FIND1-MB01 Diagnostics/ImpactT/rfdata6

## Drift
rm -f DriftTest/OPAL/*.T7
rm -f DriftTest/OPAL/TwoStep.h5
ln -s ${PWD}/DriftTest/TwoStep.h5 DriftTest/OPAL/TwoStep.h5
ln -s ${PWD}/FIELDMAPS/OPAL/FINEG-MSL10.T7 DriftTest/OPAL/FINEG-MSL10.T7

rm -f DriftTest/ImpactT/*.T7
rm -f DriftTest/ImpactT/TwoStep.h5
ln -s ${PWD}/DriftTest/TwoStep.h5 DriftTest/ImpactT/TwoStep.h5
ln -s ${PWD}/FIELDMAPS/ImpactT/FINEG-MSL10.T7 DriftTest/ImpactT/1T2.T7

## ExternalFieldTest/SpaceCharge
rm -f ExternalFieldTest/SpaceCharge/OPAL/TwoStep.h5
rm -f ExternalFieldTest/SpaceCharge/OPAL/*.T7
ln -s ${PWD}/ExternalFieldTest/TwoStep.h5 ExternalFieldTest/SpaceCharge/OPAL/TwoStep.h5
ln -s ${PWD}/FIELDMAPS/OPAL/FINEG-MSL10.T7 ExternalFieldTest/SpaceCharge/OPAL/FINEG-MSL10.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB01-MSL.T7 ExternalFieldTest/SpaceCharge/OPAL/FINLB01-MSL.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB01-RACF.T7 ExternalFieldTest/SpaceCharge/OPAL/FINLB01-RACF.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB01-RACH.T7 ExternalFieldTest/SpaceCharge/OPAL/FINLB01-RACH.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB02-RAC.T7 ExternalFieldTest/SpaceCharge/OPAL/FINLB02-RAC.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB02-MSLAC.T7 ExternalFieldTest/SpaceCharge/OPAL/FINLB02-MSLAC.T7

rm -f ExternalFieldTest/SpaceCharge/ImpactT/TwoStep.h5
rm -f ExternalFieldTest/SpaceCharge/ImpactT/*.T7
rm -f ExternalFieldTest/SpaceCharge/ImpactT/rfdata*
ln -s ${PWD}/ExternalFieldTest/TwoStep.h5 ExternalFieldTest/SpaceCharge/ImpactT/TwoStep.h5
ln -s ${PWD}/FIELDMAPS/ImpactT/FINEG-MSL10.T7 ExternalFieldTest/SpaceCharge/ImpactT/1T2.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB01-MSL.T7 ExternalFieldTest/SpaceCharge/ImpactT/1T5.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB01-RACF.T7 ExternalFieldTest/SpaceCharge/ImpactT/1T3.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB01-RACH.T7 ExternalFieldTest/SpaceCharge/ImpactT/1T4.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-RAC_1 ExternalFieldTest/SpaceCharge/ImpactT/rfdata2
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-RAC_2 ExternalFieldTest/SpaceCharge/ImpactT/rfdata3
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-RAC_3 ExternalFieldTest/SpaceCharge/ImpactT/rfdata4
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-RAC_4 ExternalFieldTest/SpaceCharge/ImpactT/rfdata5
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-MSLAC.T7 ExternalFieldTest/SpaceCharge/ImpactT/1T6.T7

## ExternalFieldTest/NoSpaceCharge
rm -f ExternalFieldTest/NoSpaceCharge/OPAL/TwoStep.h5
rm -f ExternalFieldTest/NoSpaceCharge/OPAL/*.T7
ln -s ${PWD}/ExternalFieldTest/TwoStep.h5 ExternalFieldTest/NoSpaceCharge/OPAL/TwoStep.h5
ln -s ${PWD}/FIELDMAPS/OPAL/FINEG-MSL10.T7 ExternalFieldTest/NoSpaceCharge/OPAL/FINEG-MSL10.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB01-MSL.T7 ExternalFieldTest/NoSpaceCharge/OPAL/FINLB01-MSL.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB01-RACF.T7 ExternalFieldTest/NoSpaceCharge/OPAL/FINLB01-RACF.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB01-RACH.T7 ExternalFieldTest/NoSpaceCharge/OPAL/FINLB01-RACH.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB02-RAC.T7 ExternalFieldTest/NoSpaceCharge/OPAL/FINLB02-RAC.T7
ln -s ${PWD}/FIELDMAPS/OPAL/FINLB02-MSLAC.T7 ExternalFieldTest/NoSpaceCharge/OPAL/FINLB02-MSLAC.T7

rm -f ExternalFieldTest/NoSpaceCharge/ImpactT/TwoStep.h5
rm -f ExternalFieldTest/NoSpaceCharge/ImpactT/*.T7
rm -f ExternalFieldTest/NoSpaceCharge/ImpactT/rfdata*
ln -s ${PWD}/ExternalFieldTest/TwoStep.h5 ExternalFieldTest/NoSpaceCharge/ImpactT/TwoStep.h5
ln -s ${PWD}/FIELDMAPS/ImpactT/FINEG-MSL10.T7 ExternalFieldTest/NoSpaceCharge/ImpactT/1T2.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB01-MSL.T7 ExternalFieldTest/NoSpaceCharge/ImpactT/1T5.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB01-RACF.T7 ExternalFieldTest/NoSpaceCharge/ImpactT/1T3.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB01-RACH.T7 ExternalFieldTest/NoSpaceCharge/ImpactT/1T4.T7
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-RAC_1 ExternalFieldTest/NoSpaceCharge/ImpactT/rfdata2
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-RAC_2 ExternalFieldTest/NoSpaceCharge/ImpactT/rfdata3
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-RAC_3 ExternalFieldTest/NoSpaceCharge/ImpactT/rfdata4
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-RAC_4 ExternalFieldTest/NoSpaceCharge/ImpactT/rfdata5
ln -s ${PWD}/FIELDMAPS/ImpactT/FINLB02-MSLAC.T7 ExternalFieldTest/NoSpaceCharge/ImpactT/1T6.T7



