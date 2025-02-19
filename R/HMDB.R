#' Human Metabolome Database 
#'
#' This dataset contains fields from the HMDB database 
#' The variables are as follows:
#'
#' \enumerate{
#'  \item WebAddress web address for database (when combined with the 
#'  Unique_DB_ID entry this will provide a direct link to the online database
#'  entry e.g. Web address : "www.hmdb.ca/metabolites/" +
#'  Unique_DB_ID : "HMDB00001" = \url{http://www.hmdb.ca/metabolites/HMDB00001})
#'  \item Unique_DB_ID. HMDB unique reference number (HMDB00001 -- HMDB61388)
#'  \item name. HMDB entry name
#'  \item monoisotopic_weight. Monoisotopic weight of HMDB entry 
#'  (1.007825 -- 1499.87)
#'  \item SMILES. Canonical SMILES code of HMDB entry
#'  \item molecular_formula. Molecular formula of HMDB entry 
#'  \item cas_number. CAS registry number of HMDB entry
#'  \item compound_class. compound class information utilized by annoCompoundClass function.
#'  \item - 295. (n=286) all molecular descriptors extracted from the chemistry development kit using the rcdk package and a 1024bit chemical fingerprint obtained using the ChemmineR package in the form of an index rather than full length 1024 binary representation. see table below for details:
#'  }
#'  
#' \tabular{lll}{
#' \strong{category} \tab \strong{className} \tab \strong{molDescNames}
#' \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nSmallRings     \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nAromRings      \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nRingBlocks     \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nAromBlocks     \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nRings3         \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nRings4         \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nRings5         \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nRings6         \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nRings7         \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nRings8         \cr
#' topological    \tab SmallRingDescriptor                     \tab MD_nRings9         \cr
#' topological    \tab FractionalPSADescriptor                 \tab MD_tpsaEfficiency  \cr
#' topological    \tab ZagrebIndexDescriptor                   \tab MD_Zagreb          \cr
#' constitutional \tab XLogPDescriptor                         \tab MD_XLogP           \cr
#' topological    \tab WienerNumbersDescriptor                 \tab MD_WPATH           \cr
#' topological    \tab WienerNumbersDescriptor                 \tab MD_WPOL            \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Wlambda1.unity  \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Wlambda2.unity  \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Wlambda3.unity  \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Wnu1.unity      \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Wnu2.unity      \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Wgamma1.unity   \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Wgamma2.unity   \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Wgamma3.unity   \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Weta1.unity     \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Weta2.unity     \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_Weta3.unity     \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_WT.unity        \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_WA.unity        \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_WV.unity        \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_WK.unity        \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_WG.unity        \cr
#' hybrid         \tab WHIMDescriptor                          \tab MD_WD.unity        \cr
#' topological    \tab WeightedPathDescriptor                  \tab MD_WTPT.1          \cr
#' topological    \tab WeightedPathDescriptor                  \tab MD_WTPT.2          \cr
#' topological    \tab WeightedPathDescriptor                  \tab MD_WTPT.3          \cr
#' topological    \tab WeightedPathDescriptor                  \tab MD_WTPT.4          \cr
#' topological    \tab WeightedPathDescriptor                  \tab MD_WTPT.5          \cr
#' constitutional \tab WeightDescriptor                        \tab MD_MW              \cr
#' topological    \tab VAdjMaDescriptor                        \tab MD_VAdjMat         \cr
#' topological    \tab VABCDescriptor                          \tab MD_VABC            \cr
#' topological    \tab TPSADescriptor                          \tab MD_TopoPSA         \cr
#' constitutional \tab RuleOfFiveDescriptor                    \tab MD_LipinskiFailures\cr
#' constitutional \tab RotatableBondsCountDescriptor           \tab MD_nRotB           \cr
#' topological    \tab PetitjeanShapeIndexDescriptor           \tab MD_topoShape       \cr
#' topological    \tab PetitjeanShapeIndexDescriptor           \tab MD_geomShape       \cr
#' topological    \tab PetitjeanNumberDescriptor               \tab MD_PetitjeanNumber \cr
#' geometrical    \tab MomentOfInertiaDescriptor               \tab MD_MOMI.X          \cr
#' geometrical    \tab MomentOfInertiaDescriptor               \tab MD_MOMI.Y          \cr
#' geometrical    \tab MomentOfInertiaDescriptor               \tab MD_MOMI.Z          \cr
#' geometrical    \tab MomentOfInertiaDescriptor               \tab MD_MOMI.XY         \cr
#' geometrical    \tab MomentOfInertiaDescriptor               \tab MD_MOMI.XZ         \cr
#' geometrical    \tab MomentOfInertiaDescriptor               \tab MD_MOMI.YZ         \cr
#' geometrical    \tab MomentOfInertiaDescriptor               \tab MD_MOMI.R          \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.11         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.12         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.13         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.14         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.22         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.23         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.24         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.33         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.34         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEC.44         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEO.11         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEO.12         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEO.22         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEN.11         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEN.12         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEN.13         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEN.22         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEN.23         \cr
#' topological    \tab MDEDescriptor                           \tab MD_MDEN.33         \cr
#' constitutional \tab MannholdLogPDescriptor                  \tab MD_MLogP           \cr
#' constitutional \tab LongestAliphaticChainDescriptor         \tab MD_nAtomLAC        \cr
#' geometrical    \tab LengthOverBreadthDescriptor             \tab MD_LOBMAX          \cr
#' geometrical    \tab LengthOverBreadthDescriptor             \tab MD_LOBMIN          \cr
#' constitutional \tab LargestPiSystemDescriptor               \tab MD_nAtomP          \cr
#' constitutional \tab LargestChainDescriptor                  \tab MD_nAtomLC         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sLi         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssBe        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssssBe      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssBH        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssB        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssssB       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sCH3        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dCH2        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssCH2       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.tCH         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dsCH        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aaCH        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssCH       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ddC         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.tsC         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dssC        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aasC        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aaaC        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssssC       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sNH3        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sNH2        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssNH2       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dNH         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssNH        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aaNH        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.tN          \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssNH       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dsN         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aaN         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssN        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ddsN        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aasN        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssssN       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sOH         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dO          \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssO         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aaO         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sF          \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sSiH3       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssSiH2      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssSiH      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssssSi      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sPH2        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssPH        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssP        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dsssP       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssssP      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sSH         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dS          \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssS         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aaS         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dssS        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ddssS       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sCl         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sGeH3       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssGeH2      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssGeH      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssssGe      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sAsH2       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssAsH       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssAs       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssdAs      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssssAs     \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sSeH        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dSe         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssSe        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.aaSe        \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.dssSe       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ddssSe      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sBr         \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sSnH3       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssSnH2      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssSnH      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssssSn      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sI          \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sPbH3       \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssPbH2      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.sssPbH      \cr
#' topological    \tab KierHallSmartsDescriptor                \tab MD_khs.ssssPb      \cr
#' topological    \tab KappaShapeIndicesDescriptor             \tab MD_Kier1           \cr
#' topological    \tab KappaShapeIndicesDescriptor             \tab MD_Kier2           \cr
#' topological    \tab KappaShapeIndicesDescriptor             \tab MD_Kier3           \cr
#' topological    \tab HybridizationRatioDescriptor            \tab MD_HybRatio        \cr
#' electronic     \tab HBondDonorCountDescriptor               \tab MD_nHBDon          \cr
#' electronic     \tab HBondAcceptorCountDescriptor            \tab MD_nHBAcc          \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAV.1          \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAV.2          \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAV.3          \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAVH.1         \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAVH.2         \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAVH.3         \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAV.4          \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAV.5          \cr
#' geometrical    \tab GravitationalIndexDescriptor            \tab MD_GRAV.6          \cr
#' topological    \tab FragmentComplexityDescriptor            \tab MD_fragC           \cr
#' topological    \tab FMFDescriptor                           \tab MD_FMF             \cr
#' topological    \tab EccentricConnectivityIndexDescriptor    \tab MD_ECCEN           \cr
#' electronic     \tab CPSADescriptor                          \tab MD_PPSA.1          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_PPSA.2          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_PPSA.3          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_PNSA.1          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_PNSA.2          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_PNSA.3          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_DPSA.1          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_DPSA.2          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_DPSA.3          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_FPSA.1          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_FPSA.2          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_FPSA.3          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_FNSA.1          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_FNSA.2          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_FNSA.3          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_WPSA.1          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_WPSA.2          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_WPSA.3          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_WNSA.1          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_WNSA.2          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_WNSA.3          \cr
#' electronic     \tab CPSADescriptor                          \tab MD_RPCG            \cr
#' electronic     \tab CPSADescriptor                          \tab MD_RNCG            \cr
#' electronic     \tab CPSADescriptor                          \tab MD_RPCS            \cr
#' electronic     \tab CPSADescriptor                          \tab MD_RNCS            \cr
#' electronic     \tab CPSADescriptor                          \tab MD_THSA            \cr
#' electronic     \tab CPSADescriptor                          \tab MD_TPSA            \cr
#' electronic     \tab CPSADescriptor                          \tab MD_RHSA            \cr
#' electronic     \tab CPSADescriptor                          \tab MD_RPSA            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_SP.0            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_SP.1            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_SP.2            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_SP.3            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_SP.4            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_SP.5            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_SP.6            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_SP.7            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_VP.0            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_VP.1            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_VP.2            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_VP.3            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_VP.4            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_VP.5            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_VP.6            \cr
#' topological    \tab ChiPathDescriptor                       \tab MD_VP.7            \cr
#' topological    \tab ChiPathClusterDescriptor                \tab MD_SPC.4           \cr
#' topological    \tab ChiPathClusterDescriptor                \tab MD_SPC.5           \cr
#' topological    \tab ChiPathClusterDescriptor                \tab MD_SPC.6           \cr
#' topological    \tab ChiPathClusterDescriptor                \tab MD_VPC.4           \cr
#' topological    \tab ChiPathClusterDescriptor                \tab MD_VPC.5           \cr
#' topological    \tab ChiPathClusterDescriptor                \tab MD_VPC.6           \cr
#' topological    \tab ChiClusterDescriptor                    \tab MD_SC.3            \cr
#' topological    \tab ChiClusterDescriptor                    \tab MD_SC.4            \cr
#' topological    \tab ChiClusterDescriptor                    \tab MD_SC.5            \cr
#' topological    \tab ChiClusterDescriptor                    \tab MD_SC.6            \cr
#' topological    \tab ChiClusterDescriptor                    \tab MD_VC.3            \cr
#' topological    \tab ChiClusterDescriptor                    \tab MD_VC.4            \cr
#' topological    \tab ChiClusterDescriptor                    \tab MD_VC.5            \cr
#' topological    \tab ChiClusterDescriptor                    \tab MD_VC.6            \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_SCH.3           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_SCH.4           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_SCH.5           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_SCH.6           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_SCH.7           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_VCH.3           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_VCH.4           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_VCH.5           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_VCH.6           \cr
#' topological    \tab ChiChainDescriptor                      \tab MD_VCH.7           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C1SP1           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C2SP1           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C1SP2           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C2SP2           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C3SP2           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C1SP3           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C2SP3           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C3SP3           \cr
#' topological    \tab CarbonTypesDescriptor                   \tab MD_C4SP3           \cr
#' electronic     \tab BPolDescriptor                          \tab MD_bpol            \cr
#' constitutional \tab BondCountDescriptor                     \tab MD_nB              \cr
#' hybrid         \tab BCUTDescriptor                          \tab MD_BCUTw.1l        \cr
#' hybrid         \tab BCUTDescriptor                          \tab MD_BCUTw.1h        \cr
#' hybrid         \tab BCUTDescriptor                          \tab MD_BCUTc.1l        \cr
#' hybrid         \tab BCUTDescriptor                          \tab MD_BCUTc.1h        \cr
#' hybrid         \tab BCUTDescriptor                          \tab MD_BCUTp.1l        \cr
#' hybrid         \tab BCUTDescriptor                          \tab MD_BCUTp.1h        \cr
#' constitutional \tab BasicGroupCountDescriptor               \tab MD_nBase           \cr
#' topological    \tab AutocorrelationDescriptorPolarizability \tab MD_ATSp1           \cr
#' topological    \tab AutocorrelationDescriptorPolarizability \tab MD_ATSp2           \cr
#' topological    \tab AutocorrelationDescriptorPolarizability \tab MD_ATSp3           \cr
#' topological    \tab AutocorrelationDescriptorPolarizability \tab MD_ATSp4           \cr
#' topological    \tab AutocorrelationDescriptorPolarizability \tab MD_ATSp5           \cr
#' topological    \tab AutocorrelationDescriptorMass           \tab MD_ATSm1           \cr
#' topological    \tab AutocorrelationDescriptorMass           \tab MD_ATSm2           \cr
#' topological    \tab AutocorrelationDescriptorMass           \tab MD_ATSm3           \cr
#' topological    \tab AutocorrelationDescriptorMass           \tab MD_ATSm4           \cr
#' topological    \tab AutocorrelationDescriptorMass           \tab MD_ATSm5           \cr
#' topological    \tab AutocorrelationDescriptorCharge         \tab MD_ATSc1           \cr
#' topological    \tab AutocorrelationDescriptorCharge         \tab MD_ATSc2           \cr
#' topological    \tab AutocorrelationDescriptorCharge         \tab MD_ATSc3           \cr
#' topological    \tab AutocorrelationDescriptorCharge         \tab MD_ATSc4           \cr
#' topological    \tab AutocorrelationDescriptorCharge         \tab MD_ATSc5           \cr
#' constitutional \tab AtomCountDescriptor                     \tab MD_nAtom           \cr
#' constitutional \tab AromaticBondsCountDescriptor            \tab MD_nAromBond       \cr
#' constitutional \tab AromaticAtomsCountDescriptor            \tab MD_naAromAtom      \cr
#' electronic     \tab APolDescriptor                          \tab MD_apol            \cr
#' constitutional \tab ALOGPDescriptor                         \tab MD_ALogP           \cr
#' constitutional \tab ALOGPDescriptor                         \tab MD_ALogp2          \cr
#' constitutional \tab ALOGPDescriptor                         \tab MD_AMR             \cr
#' constitutional \tab AcidicGroupCountDescriptor              \tab MD_nAcid           
#' } 
#' @docType data
#' @keywords datasets
#' @name HMDB
#' @usage data(HMDB)
#' @source \url{http://www.hmdb.ca/downloads}
#' @format A data frame with 39192 rows and 294 variables
NULL
