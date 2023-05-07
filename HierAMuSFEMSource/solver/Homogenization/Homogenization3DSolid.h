#pragma once

#include "HomogeniztaionBase.h"
#include "MatrixTypes.h"

namespace HierAMuS {

class Homogenization3DSolid : public HomogenizationBase {
public:
  Homogenization3DSolid();
  Homogenization3DSolid(const Homogenization3DSolid &other);
  ~Homogenization3DSolid(){};

  void init(PointerCollection &pointers, ParameterList &parameters) override;
  void computeAMatrix(PointerCollection &pointers) override;
  auto getNumberOfStrains() -> indexType override;

  auto getAMatrix() -> Types::MatrixXX<prec> & override;
  auto getDv() -> prec override;

  auto getType() -> indexType override { return 3; };

  auto getDisplacementIncrement(Types::VectorX<prec> &strainIncrement)
      -> Types::VectorX<prec> override;

  void setPeriodicDisplacements(PointerCollection &pointers,
                                Types::VectorX<prec> &strains,
                                Types::VectorX<prec> &strainIncrement) override {};
  
  void toFile(PointerCollection &pointers, std::ofstream &out) override;
  void fromFile(PointerCollection &pointers, std::ifstream &in) override;

private:
  void setDisplacementBoundaryConditions(PointerCollection &pointers);


  Types::MatrixXX<prec> homogenizationMatrix;

  indexType meshIdDisp;
  indexType dispOrder;
  indexType bctype;

  std::vector<indexType> leftFaces, rightFaces, topFaces, bottomFaces, frontFaces, backFaces ;

  Types::Vector3<prec> xmin, xmax;
};
} // namespace FEMProject