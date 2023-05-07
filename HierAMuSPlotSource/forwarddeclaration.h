// Copyright 2016 - 2023 Simon Klarmann <simon.klarmann@gmail.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
//#ifdef BUILDCPPFEMFEM


namespace HierAMuS {

//-----------------control-------------------------------------------------
class FEMProgram;
struct InfoData;
class OutputHandler;
class stringCommandHandler;
template <typename prec> class Timer;
//-----------------ControlProgram------------------------------------------
class controlProgram;
    //-----------------commands------------------------------------------------
namespace Commands {
class InputParser;
class GenericCommand;
} // namespace Commands

//-----------------elementFormulations-----------------------------------
namespace Elementformulations {
class ElementFormulation201;
class GenericElementFormulation;
} // namespace Elementformulations

//-----------------equations-----------------------------------------------
class DegreeOfFreedom;
class EquationHandler;
class GenericNodes;
class NodeSet;

//-----------------finiteElements------------------------------------------
namespace FiniteElement {
class ElementList;
class GenericFiniteElement;
class beamInterfaceElement2D;
class LinearTriangleElement;
class LinearPrism;
class QuadrilateralNodal;
} // namespace FiniteElement
//-----------------geometry------------------------------------------------
namespace Geometry {
class GeometryData;
class Vertex;
class Edges;
class LinearEdge;
class Faces;
class LinearTriangle;
class LinearQuadrilateral;
class Volumes;
class LinearBrick;
class Special;
class BeamInterface2D;
} // namespace Geometry
  //-----------------loads---------------------------------------------------
class DofLoad;
class GenericLoad;
class loadList;
class PropfunctionHandler;
class Timefunction;

//-----------------materials-----------------------------------------------
namespace Materials {
class GenericMaterialFormulation;
class Material;
class MaterialList;
class MaterialFormulationList;
class ElementFormulationList;
} // namespace Materials

//-----------------math----------------------------------------------------
class Plane;
class Userconstants;

//-----------------plot----------------------------------------------------
class vtkPlotInterface;

class vtkPlotInterfaceBase;

//-----------------pointercollection---------------------------------------
class PointerCollection;

//-----------------shapefunctions------------------------------------------
class IntegrationPointsManagement;
class IntegrationPointsBase;
class GaussPoints1D;
class GaussPoints2D;
class GaussPoints3D;

//-----------------solver--------------------------------------------------
class EigenPardisoLDLT;
class EigenPardisoLLT;
class EigenPardisoLU;
class EigenSimplicialLDLT;
class EigenSimplicialLLT;
class EigenSparseLU;
class EigenSparseQR;
class GenericSolutionState;
class GenericSolver;
class StaticSolutionState;

} // namespace HierAMuS
//#endif // blub
