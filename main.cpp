#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/jet.h>
#include <igl/doublearea.h>
#include <igl/principal_curvature.h>
#include <igl/lscm.h>
#include <igl/boundary_loop.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle/triangulate.h>
#include <igl/sort.h>
#include <igl/sort_vectors_ccw.h>
#include <igl/remove_unreferenced.h>

#include <map>
#include <vector>
#include <set>
#include <algorithm>
using namespace std;
using namespace Eigen;
Eigen::MatrixXd V;
Eigen::MatrixXi F; 
vector<bool> marked;
vector<bool> deletedTriangles;

#define REMOVE_VERTEX_NULL 0

typedef std::set<std::pair<double,int> > PriorityQueue;

void removeVertex(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXi &F2, int vi)
{
  vector<vector<int> > VF;
  vector<vector<int> > VFi;
  igl::vertex_triangle_adjacency(V, F, VF, VFi);
  
  int N = VF[vi].size();
  map<int, int> conn;
  
  int first;
  for (int i = 0; i < VF[vi].size(); i++){
    vector<int> a;
    int loc_vi;
    for (int j = 0; j < 3; j++) {
      if (F(VF[vi][i], j) != vi) a.push_back(F(VF[vi][i], j));
      else loc_vi = j;
    }
    if (loc_vi == 1) conn[a[1]] = a[0];
    else conn[a[0]] = a[1];
    first = a[1];
  }

  int last = conn[first];
  vector<int> circle;
  circle.push_back(last);
  while (last != first) {
    last = conn[last];
    circle.push_back(last);
  }


  for(auto it : circle) {
    marked[it] = true;
  }
  marked[vi] = true;
  
  int K = circle.size();
  int j_last = circle[K-1];
  VectorXd vk_last = V.row(j_last) - V.row(vi);
  double thetak = 0;
  vector<double> thetacum;
  vector<double> rs;
  for (int i = 0; i < K; i++) {
    VectorXd vk = V.row(circle[i]) - V.row(vi);
    rs.push_back(vk.norm());
    double theta = acos(vk.dot(vk_last)/(vk.norm()*vk_last.norm()));
    thetak += theta;
    thetacum.push_back(thetak);
    j_last = circle[i];
    vk_last = vk;
  }

  MatrixXd V2d(K,2);
  MatrixXi E2d(K, 2);
  MatrixXd H2d(0,0);

  double a = 2.* igl::PI / thetak;

  for (int i = 0; i < K; i++) {
    double x = pow(rs[i], a) * cos(thetacum[i]*a);
    double y = pow(rs[i], a) * sin(thetacum[i]*a);

    V2d(i, 0) = x;
    V2d(i, 1) = y; 
    E2d(i, 0) = i;
    if (i == K-1) E2d(i, 1) = 0;
    else E2d(i, 1) = i+1;
  }
  //H2d << 0 , 0;
  MatrixXd Vt;
  MatrixXi Ft;

  igl::triangle::triangulate(V2d,E2d,H2d,"pQ",Vt,Ft);
  MatrixXi Ftsorted;
  MatrixXi I;
  
  F2.resize(F.rows(),3);

  int m = 0;
  for(int f = 0;f<F.rows();f++)
  {
    if(
      F(f,0) == vi || 
      F(f,1) == vi || 
      F(f,2) == vi)
    {
      continue;
    } else {
      F2.row(m) = F.row(f);
      m++;
    }
  }

  for (int i=0; i < Ft.rows();i++)
  {
    F2(m, 0) = circle[Ft(i, 0)];
    F2(m, 1) = circle[Ft(i, 1)];
    F2(m, 2) = circle[Ft(i, 2)];
    //cout << "new face: " << F2.row(m) << endl;
    m++;
  }
  F2.conservativeResize(m,F2.cols());

}

void resetMaker() {
  for (int i = 0; i < V.rows(); i++) {
    marked[i] = false;
  }
}

auto StarArea = [](Eigen::MatrixXd V, Eigen::MatrixXi F) 
{ 
  Eigen::VectorXd triA;
  Eigen::VectorXd starA = Eigen::VectorXd::Zero(V.rows());
  igl::doublearea(V, F, triA);
  // Loop over faces
  for(int i = 0;i<F.rows();i++)
  {
    // loop around face
    for(int j = 0;j<F.cols();j++)
    {
      starA(F(i,j)) += triA(i);
    }
  }
  starA.array() /= starA.maxCoeff();

  return starA;
};

auto curvatureScore = [](Eigen::MatrixXd V, Eigen::MatrixXi F) 
{ 
  Eigen::MatrixXd PD1,PD2;
  Eigen::VectorXd PV1,PV2;
  igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
  Eigen::VectorXd H = (PV1.array().abs()+PV2.array().abs());
  H.array() /= H.maxCoeff();
  return H;
};

void calcPriority(Eigen::MatrixXd V, Eigen::MatrixXi F, PriorityQueue& Q, std::vector<PriorityQueue::iterator >& Qit) 
{ 
  Eigen::VectorXd area(V.rows());
  area = StarArea(V, F);
  Eigen::VectorXd H = curvatureScore(V, F);
  Eigen::VectorXd score = 0.5*(H+area);

  Qit.resize(V.rows());
  for (int i = 0; i < V.rows(); i++) {
    Qit[i] = Q.insert(std::pair<double, int>(score(i), i)).first;
  }
};


int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("bunny.off", V, F);
  MatrixXd V2;
  MatrixXi F2;
  VectorXi I2;
  int L = floor(log2(V.rows()));

  cout << V.rows() << endl;

  PriorityQueue Q;
  std::vector<PriorityQueue::iterator > Qit;
  for (int i = 0; i < V.rows(); i++) marked.push_back(false);
  calcPriority(V,F, Q, Qit);

  const auto &pre_draw = [&](igl::viewer::Viewer & viewer)->bool
  {
    if(viewer.core.is_animating)
    {
      if (!Q.empty()) {
        std::pair<double,int> p = *(Q.begin());
        Q.erase(Q.begin());
        int vi = p.second;
        //cout << Q.size() << endl;
        if (!marked[vi]) {
          removeVertex(V, F, F2, vi);
          F = F2;
          //cout << F.rows();
          viewer.data.clear();
          viewer.data.set_mesh(V,F);  
          viewer.data.set_face_based(true);
        }
      } else if (L > 0) {
        cout << "L = " << L << endl;
        igl::remove_unreferenced(V, F, V2, F2, I2);
        V = V2;
        F = F2;
        calcPriority(V,F, Q, Qit);
        resetMaker();
        L--;
      } else {
        viewer.core.is_animating ^= 1;
      }
    }
    return false;
  };

  const auto &key_down =
    [&](igl::viewer::Viewer &viewer,unsigned char key,int mod)->bool
  {
    switch(key)
    {
      case ' ':
        viewer.core.is_animating ^= 1;
        break;
      
      default:
        return false;
    }
    return true;
  };

  /*
  vector<vector<int> > VF;
  vector<vector<int> > VFi;
  igl::vertex_triangle_adjacency(V, F, VF, VFi);
  for (int i = 0; i < VF[p.second].size(); i++) 
  {
    int tri = VF[p.second][i];
    cout << tri << endl;
  }
  */
  // remove all IGL_COLLAPSE_EDGE_NULL faces
 
    // Plot the mesh
  igl::viewer::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.core.is_animating = false;
  viewer.callback_key_down = key_down;
  viewer.callback_pre_draw = pre_draw;

  viewer.launch();
}
