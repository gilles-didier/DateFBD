ComputeDiagnosticMCMC.o: ComputeDiagnosticMCMC.c Utils.h MyR.h Tree.h \
 Fossil.h FossilInt.h Model.h Uncertainty.h Distribution.h \
 MCMCImportanceSampling.h
ComputeExtinctionComparison.o: ComputeExtinctionComparison.c Utils.h \
 MyR.h Tree.h Fossil.h FossilInt.h Model.h Uncertainty.h Distribution.h \
 MCMCImportanceSampling.h
ComputeExtinctionDistribution.o: ComputeExtinctionDistribution.c Utils.h \
 MyR.h Tree.h Fossil.h FossilInt.h Model.h Uncertainty.h Distribution.h \
 MCMCImportanceSampling.h
ComputeExtinctionDistributionSingleTree.o: \
 ComputeExtinctionDistributionSingleTree.c Utils.h MyR.h Tree.h Fossil.h \
 FossilInt.h Model.h Uncertainty.h Distribution.h MinimizeNLOpt.h \
 DrawTreeCairo.h DrawTreeGeneric.h DrawTreePSTricks.h DrawTreeTikz.h \
 DrawFossilInt.h DrawDensity.h DrawTimePosition.h \
 MCMCImportanceSampling.h
ComputeExtinctionUpperBound.o: ComputeExtinctionUpperBound.c Utils.h \
 MyR.h Tree.h Fossil.h FossilInt.h Model.h Uncertainty.h Distribution.h \
 MCMCImportanceSampling.h
ComputeSpeciationDistribution.o: ComputeSpeciationDistribution.c Utils.h \
 MyR.h Tree.h Fossil.h FossilInt.h Model.h Uncertainty.h Distribution.h \
 MCMCImportanceSampling.h
ComputeSpeciationDistributionSingleTree.o: \
 ComputeSpeciationDistributionSingleTree.c Utils.h MyR.h Tree.h Fossil.h \
 FossilInt.h Model.h Uncertainty.h Distribution.h MinimizeNLOpt.h \
 DrawTreeCairo.h DrawTreeGeneric.h DrawTreePSTricks.h DrawTreeTikz.h \
 DrawFossilInt.h DrawDensity.h MCMCImportanceSampling.h
Distribution.o: Distribution.c Utils.h MyR.h Distribution.h
DrawDensity.o: DrawDensity.c DrawDensity.h DrawTreeGeneric.h Tree.h \
 Utils.h MyR.h Distribution.h
DrawFossil.o: DrawFossil.c DrawFossil.h DrawTreeGeneric.h Tree.h Utils.h \
 MyR.h Fossil.h
DrawFossilInt.o: DrawFossilInt.c DrawFossilInt.h DrawTreeGeneric.h Tree.h \
 Utils.h MyR.h FossilInt.h Fossil.h
DrawFossilTree.o: DrawFossilTree.c Utils.h MyR.h Tree.h Fossil.h \
 FossilInt.h DrawTreeCairo.h DrawTreeGeneric.h DrawTreePSTricks.h \
 DrawTreeTikz.h DrawFossilInt.h DrawDensity.h Distribution.h \
 DrawTimePosition.h
DrawTimePosition.o: DrawTimePosition.c DrawTimePosition.h \
 DrawTreeGeneric.h Tree.h Utils.h MyR.h Distribution.h
DrawTreeCairo.o: DrawTreeCairo.c /usr/include/cairo/cairo.h \
 /usr/include/cairo/cairo-version.h /usr/include/cairo/cairo-features.h \
 /usr/include/cairo/cairo-deprecated.h /usr/include/cairo/cairo-ps.h \
 /usr/include/cairo/cairo.h /usr/include/cairo/cairo-pdf.h \
 /usr/include/cairo/cairo-svg.h Utils.h MyR.h DrawTreeCairo.h Tree.h \
 DrawTreeGeneric.h
DrawTreeFossil.o: DrawTreeFossil.c Utils.h MyR.h Tree.h Fossil.h \
 DrawTreeCairo.h DrawTreeGeneric.h DrawTreePSTricks.h DrawTreeTikz.h \
 DrawFossilInt.h FossilInt.h
DrawTreeGeneric.o: DrawTreeGeneric.c Utils.h MyR.h DrawTreeGeneric.h \
 Tree.h
DrawTreePSTricks.o: DrawTreePSTricks.c Utils.h MyR.h DrawTreePSTricks.h \
 Tree.h DrawTreeGeneric.h
DrawTreeTikz.o: DrawTreeTikz.c Utils.h MyR.h DrawTreeTikz.h Tree.h \
 DrawTreeGeneric.h
EstimateRates.o: EstimateRates.c Utils.h MyR.h Tree.h Fossil.h \
 FossilInt.h Model.h Uncertainty.h Distribution.h MinimizeNLOpt.h
ExtinctionCI.o: ExtinctionCI.c Utils.h MyR.h Tree.h ExtinctionCI.h \
 Fossil.h
Fossil.o: Fossil.c Utils.h MyR.h Fossil.h Tree.h
FossilInt.o: FossilInt.c Utils.h MyR.h FossilInt.h Tree.h Fossil.h
MCMCImportanceSampling.o: MCMCImportanceSampling.c Utils.h MyR.h \
 Uncertainty.h Tree.h Fossil.h Model.h Distribution.h \
 MCMCImportanceSampling.h FossilInt.h
MinimizeNLOpt.o: MinimizeNLOpt.c Utils.h MyR.h MinimizeNLOpt.h Model.h \
 Tree.h Fossil.h FossilInt.h
Model.o: Model.c Utils.h MyR.h TreeExtras.h Tree.h Model.h Fossil.h
MyR.o: MyR.c MyR.h
SampleTreeFossil.o: SampleTreeFossil.c Utils.h MyR.h Tree.h SimulTree.h \
 Fossil.h SimulFossil.h DrawTreeCairo.h DrawTreeGeneric.h \
 DrawTreePSTricks.h DrawTreeTikz.h DrawFossil.h DrawDensity.h \
 Distribution.h
SimulFossil.o: SimulFossil.c SimulFossil.h Fossil.h Utils.h MyR.h Tree.h
SimulTree.o: SimulTree.c SimulTree.h Tree.h Utils.h MyR.h
Tree.o: Tree.c Tree.h Utils.h MyR.h
TreeExtras.o: TreeExtras.c Utils.h MyR.h TreeExtras.h Tree.h
Uncertainty.o: Uncertainty.c Utils.h MyR.h TreeExtras.h Tree.h \
 Uncertainty.h Fossil.h Model.h Distribution.h
Utils.o: Utils.c Utils.h MyR.h
