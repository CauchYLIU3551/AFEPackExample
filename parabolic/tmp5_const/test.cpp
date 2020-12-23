/**
 * @file   test.cpp
 * @author Robert Lie
 * @date   Wed Feb 28 11:45:08 2007
 *
 * @brief  这是个抛物型方程的例子，非常简单，中间就没有解释太多了
 *
 */

#include <iostream>
#include <vector>
#include <math.h>
#include <iterator>
#include <algorithm>

#include <base/exceptions.h>
#include <lac/full_matrix.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>


#include <AFEPack/EasyMesh.h>
#include <AFEPack/TemplateElement.h>
#include <AFEPack/FEMSpace.h>
#include <AFEPack/BilinearOperator.h>
#include <AFEPack/Operator.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/Geometry.h>
#include <AFEPack/Miscellaneous.h>
#include <AFEPack/Functional.h>

#define DIM 2
#define PI (4.0*atan(1.0)) 

/// 初值和边值的表达式


// Try to use the global variable to get the equation with diff time t;
double t=0;
double max=0.0;
double _u_(const double * p)
{
  return sin(PI*(p[0]+t))*sin(2*PI*p[1]);
  //return sin(PI*p[0]) * sin(2*PI*p[1]);
  //return p[0]*exp(p[1]);
}

/// 右端项
double _f_(const double * p)
{
  return 5*PI*PI*_u_(p)+PI*cos(PI*(p[0]))*sin(2*PI*p[1]);
  //return 10+5*PI*PI*_u_(p);
  //return p[0]*p[1] + sin(p[1]);
}


/// 1/dt - \Delta 离散出来的矩阵
class Matrix : public L2InnerProduct<DIM,double>
{
private:
  double _dt;
public:
  Matrix(FEMSpace<double,DIM>& sp, double dt) :
    L2InnerProduct<DIM,double>(sp, sp), _dt(dt) {}
/// Get the element matrix;
  virtual void getElementMatrix(const Element<double,DIM>& e0,
                                const Element<double,DIM>& e1,
                                const ActiveElementPairIterator< DIM >::State s)
  {
    double vol = e0.templateElement().volume();
    u_int acc = algebricAccuracy();
    const QuadratureInfo<DIM>& qi = e0.findQuadratureInfo(acc);
    u_int n_q_pnt = qi.n_quadraturePoint();
    std::vector<double> jac = e0.local_to_global_jacobian(qi.quadraturePoint());
    //AFEPack::Point<DIM> test1;
    //std::vector<AFEPack::Point<DIM>> test2;
    // 
    // Here always get an error! vector<Point<DIM> >;
    // Reason: Because there is a same name class in deal.ii: Point class. While
    // using Point, it might use the dealii.Point defaultly!!!! So you just add 
    // the namespace to control the class will solve the problem;
    std::vector<AFEPack::Point<DIM>> q_pnt = e0.local_to_global(qi.quadraturePoint());
    //std::cout<<"ATTENTION!This is a test flag!!!"<<std::endl;
    std::vector<std::vector<double> > bas_val = e0.basis_function_value(q_pnt);
    std::vector<std::vector<std::vector<double> > > bas_grad = e0.basis_function_gradient(q_pnt);
    u_int n_ele_dof = e0.dof().size();
    for (u_int l = 0;l < n_q_pnt;++ l) {
      double Jxw = vol*qi.weight(l)*jac[l];
      for (u_int i = 0;i < n_ele_dof;++ i) {
        for (u_int j = 0;j < n_ele_dof;++ j) {
          elementMatrix(i,j) += Jxw*(bas_val[i][l]*bas_val[j][l]/_dt +
                                     innerProduct(bas_grad[i][l], bas_grad[j][l]));
        }
      }
    }
  }
};

int main(int argc, char * argv[])
{
  /// 准备网格
  EasyMesh mesh;
  mesh.readData(argv[1]);

  /// 准备参考单元
  TemplateGeometry<DIM> tmp_geo;
  tmp_geo.readData("triangle.tmp_geo");
  CoordTransform<DIM,DIM> crd_trs;
  crd_trs.readData("triangle.crd_trs");
  TemplateDOF<DIM> tmp_dof(tmp_geo);
  tmp_dof.readData("triangle.1.tmp_dof");
  BasisFunctionAdmin<double,DIM,DIM> bas_fun(tmp_dof);
  bas_fun.readData("triangle.1.bas_fun");

  std::vector<TemplateElement<double,DIM> > tmp_ele(1);
  tmp_ele[0].reinit(tmp_geo, tmp_dof, crd_trs, bas_fun);

  /// 定制有限元空间
  FEMSpace<double,DIM> fem_space(mesh, tmp_ele);
  u_int n_ele = mesh.n_geometry(DIM);
  fem_space.element().resize(n_ele);
  for (u_int i = 0;i < n_ele;++ i) {
    fem_space.element(i).reinit(fem_space, i, 0);
  }
  fem_space.buildElement();
  fem_space.buildDof();
  fem_space.buildDofBoundaryMark();

  /// 准备初值
  FEMFunction<double,DIM> u_h(fem_space);
  FEMFunction<double,DIM> u_h2(fem_space);
  //Operator::L2Interpolate(&_u_, u_h);

  //FEMSpace<double,DIM>& fem_space = f1.femSpace();
  FEMSpace<double,DIM>::ElementIterator the_element = fem_space.beginElement();
  FEMSpace<double,DIM>::ElementIterator end_element = fem_space.endElement();

  //double t=0;//

  for (;the_element != end_element;the_element ++) {
    const std::vector<int>& element_dof = the_element->dof();
    unsigned int n_element_dof = element_dof.size();
    for (unsigned int i = 0;i < n_element_dof;i ++) {
      unsigned int j = element_dof[i];
      const AFEPack::Point<DIM>& interp_point = fem_space.dofInfo(j).interp_point;
      u_h(j) = _u_(interp_point);
    }
  }

  /// 准备边界条件
  BoundaryFunction<double,DIM> boundary(BoundaryConditionInfo::DIRICHLET,
                                        1,
                                        &_u_);
  BoundaryConditionAdmin<double,DIM> boundary_admin(fem_space);
  boundary_admin.add(boundary);

  //double t;//

  do {
    double dt = 0.01; /// 简单起见，随手取个时间步长算了

    int num=t/dt+1;// The number of iterations;

    /// 准备线性系统的矩阵
    Matrix mat(fem_space, dt);
    mat.algebricAccuracy() = 4;
    mat.build();



    /*
    // Build the original function to help to compare;
    for (;the_element != end_element;the_element ++) {
    const std::vector<int>& element_dof2 = the_element->dof();
    unsigned int n_element_dof2 = element_dof2.size();
    for (unsigned int i = 0;i < n_element_dof2;i ++) {
      unsigned int j = element_dof2[i];
      const AFEPack::Point<DIM>& interp_point2 = fem_space.dofInfo(j).interp_point;
      u_h2(j) = _u_(interp_point2);
      }
    }
    std::stringstream result2;
    result2.setf(std::ios::fixed);
    result2.precision(0);
    result2 << "u_original_" << int(t/dt) << ".dx";
    u_h2.writeOpenDXData(result2.str());
*/
    /// 准备右端项
    Vector<double> rhs(fem_space.n_dof());
    FEMSpace<double,DIM>::ElementIterator the_ele = fem_space.beginElement();
    FEMSpace<double,DIM>::ElementIterator end_ele = fem_space.endElement();
    for (;the_ele != end_ele;++ the_ele) {
      double vol = the_ele->templateElement().volume();
      const QuadratureInfo<DIM>& qi = the_ele->findQuadratureInfo(4);
      u_int n_q_pnt = qi.n_quadraturePoint();
      std::vector<double> jac = the_ele->local_to_global_jacobian(qi.quadraturePoint());
      std::vector<AFEPack::Point<DIM>> q_pnt = the_ele->local_to_global(qi.quadraturePoint());
      std::vector<std::vector<double> > bas_val = the_ele->basis_function_value(q_pnt);

      // Give the u_h value in fact it is the value of the u^{n};
      //
      std::vector<double> u_h_val = u_h.value(bas_val, *the_ele);
      std::vector<std::vector<double> > u_h_grad = u_h.gradient(q_pnt, *the_ele);
      const std::vector<int>& ele_dof = the_ele->dof();
      u_int n_ele_dof = ele_dof.size();
      for (u_int l = 0;l < n_q_pnt;++ l) {
        double Jxw = vol*qi.weight(l)*jac[l];
        double f_val = _f_(q_pnt[l]);
        for (u_int i = 0;i < n_ele_dof;++ i) {
          // Build the rhs vector;
          rhs(ele_dof[i]) += Jxw*bas_val[i][l]*(u_h_val[l]/dt + f_val);
        }
      }

      // Get the original functions;
      const std::vector<int>& element_dof2 = the_ele ->dof();
      unsigned int n_element_dof2 = element_dof2.size();
      for (unsigned int i = 0;i < n_element_dof2;i ++) {
	unsigned int j = element_dof2[i];
	const AFEPack::Point<DIM>& interp_point2 = fem_space.dofInfo(j).interp_point;
	u_h2(j) = _u_(interp_point2);
      }
    }

    // Save the result of the original function at time = t s;
    std::stringstream result2;
    result2.setf(std::ios::fixed);
    result2.precision(0);
    result2 << "u_original_" << t/dt << ".dx";
    u_h2.writeOpenDXData(result2.str());

    /// 应用边界条件
    boundary_admin.apply(mat, u_h, rhs);

    /// 求解线性系统
    AMGSolver solver(mat);
    //solver.lazyReinit(mat);
    solver.solve(u_h, rhs, 1.0e-08, 200);
   
    // calculate the maximum absolute value of the u_h;
    double max1=0.0;
    for (int j=0;j<u_h.size();j++)
    {
      if(abs(u_h[j])>max1)
      {
        max1=abs(u_h[j]);
      }
    }
    if(abs(max1-1)>max)
    {
      max=abs(max1-1);
    }

    std::stringstream result;
    result.setf(std::ios::fixed);
    result.precision(0);
    result << "u_h_" << t/dt+1 << ".dx";
    u_h.writeOpenDXData(result.str());

    t += dt; /// 更新时间

    //std::cout << "Press ENTER to continue or CTRL+C to stop ..." << std::flush;
    //getchar();

    double error = Functional::L2Error(u_h, FunctionFunction<double>(&_u_), 3);
    std::cerr << "\nL2 error = " << error << std::endl;


    std::cout << "\n\tt = " <<  t << std::endl;
  } while (t<=2);
 
// Get the max error between 1 and the top of the u_h;
//
  std::cout<<"The max error of the calculation is : "<<max<<std::endl; 
  return 0;
}

/**
 * end of file
 *
 */
