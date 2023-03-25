#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <cmath>
#include <set>
#include <utility>
#include <Eigen/Dense>
#include <ostream>
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/per_face_normals.h>


using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V,V_new; // matrix storing vertex coordinates of the input mesh (n rows, 3 columns)
MatrixXi F,F_new; // incidence relations between faces and edges (f columns)
double threshold;

std::vector<std::vector<int>> vertex_face_adjacency;
std::vector<std::vector<int>> VFi;
Eigen::MatrixXd face_normals;

void draw_bounding_box(igl::opengl::glfw::Viewer &viewer, const MatrixXd &V)
{
    // compute the corners of the bounding box
    Vector3d m = V.colwise().minCoeff();
    Vector3d M = V.colwise().maxCoeff();
    
    MatrixXd V_box(8, 3);  // Corners of the bounding box
    MatrixXi E_box(12, 2); // edges of the bounding box
    
    V_box << m(0), m(1), m(2),
    M(0), m(1), m(2),
    M(0), M(1), m(2),
    m(0), M(1), m(2),
    m(0), m(1), M(2),
    M(0), m(1), M(2),
    M(0), M(1), M(2),
    m(0), M(1), M(2);
    
    E_box << 0, 1,
    1, 2,
    2, 3,
    3, 0,
    4, 5,
    5, 6,
    6, 7,
    7, 4,
    0, 4,
    1, 5,
    2, 6,
    7, 3;
    
    viewer.append_mesh();
    viewer.data(1).add_points(V_box, Eigen::RowVector3d(1, 0, 0));
    
    for (unsigned i = 0; i < E_box.rows(); ++i) // Plot the edges of the bounding box
        viewer.data().add_edges(
                                V_box.row(E_box(i, 0)),
                                V_box.row(E_box(i, 1)),
                                Eigen::RowVector3d(1, 0, 0));
}

bool vector_contains (vector<pair<int,int>> veclist, pair<int,int> a){
    return count(veclist.begin(), veclist.end(), a);
}
bool vector_contains (vector<int> veclist, int a){
    return count(veclist.begin(), veclist.end(), a);
}



RowVector4d face_vector(RowVector3i f1, int v) {  //calculates face vector
    //std::cout<<"started here"<<std::endl;
    RowVector3d edge_ab = V.row(f1(1))-V.row(f1(0));//(B - A) ; Edge from A to B
    //std::cout<<"1st edge "<< edge_ab<<std::endl;
    //std::cout<<"2nd edge "<<f1<< " "<< V.rows() <<" "<<V.row(f1(2))-V.row(f1(0))<<std::endl;
    RowVector3d edge_ac = V.row(f1(2))-V.row(f1(0));//(C - A) ; Edge from A to C
    
    //std::cout<<"comes here"<<std::endl;
    RowVector3d normal = edge_ab.cross(edge_ac);
    
    normal = normal.normalized();
    //    cout<<"\nNORMAL "<<normal<<endl;
    double d = -RowVector3d(V.row(v))*normal.transpose();
//        std::cout<<"ax+by+cz+d=> "<<normal(0)* V(f1(0),0)<<" + " <<normal(1)*V(f1(0),1)<<" + " <<normal(2)*V(f1(0),2)<<" + "<<d<<" = "<<normal(0)* V(f1(0),0)+ normal(1)*V(f1(0),1)+ normal(2)*V(f1(0),2)+d<<endl;
//    std::cout<<"a2+b2+c2= "<<normal(0)*normal(0) +normal(1)*normal(1) +normal(2)*normal(2) <<endl;
    //    cout<<"d= "<<d<<endl;
    return RowVector4d((normal(0)), (normal(1)), (normal(2)), (d));
//    return RowVector4d(fabs(normal(0)), fabs(normal(1)), fabs(normal(2)), (d));
}
double distance(Eigen::Vector3d a , Eigen::Vector3d b){
    return (b-a).norm();
}
double faceArea(int faceIndex) {
    // Get the vertex indices of the face
    int v1 = F(faceIndex, 0);
    int v2 = F(faceIndex, 1);
    int v3 = F(faceIndex, 2);
    
    // Get the vertex coordinates
    Eigen::Vector3d p1 = V.row(v1);
    Eigen::Vector3d p2 = V.row(v2);
    Eigen::Vector3d p3 = V.row(v3);

    // Calculate the cross product of two edges to find the area
    Eigen::Vector3d edge1 = p2 - p1;
    Eigen::Vector3d edge2 = p3 - p1;
    double area = 0.5 * edge1.cross(edge2).norm();
    
    return area;
}

//; Calculates the quadric error metric for the given vertex
MatrixXd calculate_quadric(int v) { // index of the vertex
    
    MatrixXd Q = MatrixXd::Zero(4,4);
//    vector<int> otherV;//, v2; //other vertices
//    vector<vector<int>> faces_of_v;
    //    std::cout<<
    vector<int> neigh_face;
    for(int i=0; i<F.rows(); i++)
    {
        if(F(i,0)==v)
        {
            neigh_face.push_back(i);
//            if(!vector_contains(otherV,F(i,1)))
//                otherV.push_back(F(i,1));
//            if(!vector_contains(otherV,F(i,2)))
//                otherV.push_back(F(i,2));
//            faces_of_v.push_back(otherV);
        }
        else if(F(i,1)==v)
        {
            neigh_face.push_back(i);
//            if(!vector_contains(otherV,F(i,0)))
//                otherV.push_back(F(i,0));
//            if(!vector_contains(otherV,F(i,2)))
//                otherV.push_back(F(i,2));
//            faces_of_v.push_back(otherV);
        }
        else if(F(i,2)==v)
        {
            neigh_face.push_back(i);
//            if(!vector_contains(otherV,F(i,0)))
//                otherV.push_back(F(i,0));
//            if(!vector_contains(otherV,F(i,1)))
//                otherV.push_back(F(i,1));
//            faces_of_v.push_back(otherV);
        }
//        otherV.clear();
    }
    
//    for(int i=0; i < faces_of_v.size(); i++ ){
    for(int i=0; i < neigh_face.size(); i++ ){
        MatrixXd Kp = MatrixXd::Identity(4,4);
//        int v1,v2;
//
//        v1= faces_of_v[i][0];
//        v2= faces_of_v[i][1];
        RowVector4d face = face_vector(F.row(neigh_face[i]), v);
        Kp = faceArea(neigh_face[i])*face.transpose()*face;
        Q += Kp;
    }
 
//    faces_of_v.clear();
    return Q;
}

struct Edge {
    int a , b;
    Edge(  int c ,  int d ) : a( std::min< int>(c,d) ) , b( std::max< int>(c,d) ) {}
    bool operator < ( Edge const & o ) const {   return a < o.a  ||  (a == o.a && b < o.b);  }
    bool operator == ( Edge const & o ) const {   return a == o.a  &&  b == o.b;  }
};

bool edge_contains (vector<Edge> veclist, Edge a){
    return count(veclist.begin(), veclist.end(), a);
}

std::map< Edge , std::vector< int > > to_find_boundary_edges;
vector< bool > vertex_is_collapsable;
//vector<MatrixXd> errors(V.rows());
std::vector<Edge> edge_map_;



void create_edgemap(){
    to_find_boundary_edges.clear();
    edge_map_.clear();
    for(int i=0; i<F.rows(); i++)
    {
        Edge edge1= Edge(F(i,0) , F(i,1));
        Edge edge2= Edge(F(i,1) , F(i,2));
        Edge edge3= Edge(F(i,2) , F(i,0));
        to_find_boundary_edges[edge1].push_back(i);
        to_find_boundary_edges[edge2].push_back(i);
        to_find_boundary_edges[edge3].push_back(i);
        if(!edge_contains(edge_map_, edge1))
            edge_map_.push_back(edge1);
        if(!edge_contains(edge_map_, edge2))
            edge_map_.push_back(edge2);
        if(!edge_contains(edge_map_, edge3))
            edge_map_.push_back(edge3);
    }
//    cout<<"Edge Map created!! -- No. of edges (Non - duplicate) are "<<edge_map_.size()<<endl;
//    cout<<"...."<<flush;
}

void vertex_is_collapsable_check(){
    vertex_is_collapsable.clear();
    vertex_is_collapsable.resize(F.rows());
    for(int i=0; i< F.rows(); i++)
    {
        if(to_find_boundary_edges[Edge(F(i,0),F(i,1))].size()==2)
        {
            vertex_is_collapsable[F(i,0)]=true;
            vertex_is_collapsable[F(i,1)]=true;
        }
        else
        {
            vertex_is_collapsable[F(i,0)]=false;
            vertex_is_collapsable[F(i,1)]=false;
        }
        if(to_find_boundary_edges[Edge(F(i,1),F(i,2))].size()==2)
        {
            vertex_is_collapsable[F(i,2)]=true;
            vertex_is_collapsable[F(i,1)]=true;
        }
        else
        {
            vertex_is_collapsable[F(i,2)]=false;
            vertex_is_collapsable[F(i,1)]=false;
        }
        
        if(to_find_boundary_edges[Edge(F(i,2),F(i,0))].size()==2)
        {
            vertex_is_collapsable[F(i,0)]=true;
            vertex_is_collapsable[F(i,2)]=true;
        }
        else
        {
            vertex_is_collapsable[F(i,0)]=false;
            vertex_is_collapsable[F(i,2)]=false;
        }
    }
//    cout<<"...."<<endl;
//    cout<<"Collapsable check!! -- done"<<endl;
}

//void initialize_error_quadrics(){
//    errors.clear();
//    errors.resize(V.rows());
//    for(int i=0; i<V.rows(); i++)
//    {
//        errors.at(i) = calculate_quadric(i);
//        cout<<"until "<<i<<endl;
//    }
//    cout<<"initializing error quadrics!! -- done"<<endl;
//
//}
struct Pair {
    int a , b;
    Pair(  int c ,  int d ) : a( std::min< int>(c,d) ) , b( std::max< int>(c,d) ) {}
    bool operator < ( Pair const & o ) const {   return a < o.a  ||  (a == o.a && b < o.b);  }
    bool operator == ( Pair const & o ) const {   return a == o.a  &&  b == o.b;  }
};

struct VertexPairs {
    Pair p;
    double cost;
    RowVector3d target;
//    std::vector<int> adj_faces;
    bool operator==(const VertexPairs& other) const {
            return (p == other.p);
        }
};
bool comp(const VertexPairs& a, const VertexPairs& b) {
    if( a.cost < b.cost)
        return true;
    else
        return false;
}
bool Pair_contains (vector<Pair> veclist, Pair a){
    return count(veclist.begin(), veclist.end(), a);
}
bool VertexPairs_contains (vector<VertexPairs> veclist, VertexPairs a){
    return std::find(veclist.begin(), veclist.end(), a) != veclist.end();
}
vector<int> removed_indices;
vector<VertexPairs> vp;

void initialize_vertex_pairs(){
    //todo create all the initial vertex pairs
//    vertex_pairs.clear();
    for(int i=0; i<V.rows(); i++)
    {
        for(int j=0; j<V.rows(); j++)
        {
            if(i!=j)
            {
                if(edge_contains(edge_map_, Edge(i,j)))
                {
                    VertexPairs temp = { Pair(i,j), 0.0, RowVector3d(0.0,0.,0.)};
                    if(!VertexPairs_contains(vp,temp))
                    {
                        vp.push_back(temp);
                    }
                }
                else if(distance(V.row(i),V.row(j)) < threshold)
                {
                    VertexPairs temp = { Pair(i,j), 0.0, RowVector3d(0.0,0.,0.)};
                    if(!VertexPairs_contains(vp, temp))
                    {
                        vp.push_back(temp);
                    }
                }
            }
        }
    }
//    cout<<"... "<<endl;
//    cout<<"initializin vertexpairs -- done "<<endl;
}
//vector<double> cost;
//vector<RowVector3d> target;
//vector<std::pair<double,Pair>> costof_pair;
double computecost(int v1, int v2, RowVector3d &target)
{
    //        costof_pair
    //compute cost of each pair;
    //and also store the new position
    //think how to store the new position;
    //( Q11  Q12  Q13  Q14 )        ( 0 )
    //( Q21  Q22  Q23  Q24 ) * vt = ( 0 )
    //( Q13  Q14  Q33  Q34 )        ( 0 )
    //(  0    0    0    1  )        ( 1 )
    
    MatrixXd Q = MatrixXd::Zero(4,4);
    MatrixXd tQ = MatrixXd::Zero(4,4);
    //    int v1 = vertex_pairs[i].a, v2=vertex_pairs[i].b;
//    std::cout<<"starts here "<<v1<<" , "<<v2<<"  "<<calculate_quadric(v1)<< " "<< calculate_quadric(v2) <<endl;
    tQ = calculate_quadric(v1) + calculate_quadric(v2);
//    std::cout<<"ends here "<<v1<<" , "<<v2<<endl;
    Q.row(0) = RowVector4d(tQ(0,0),tQ(0,1),tQ(0,2),tQ(0,3));
    Q.row(1) = RowVector4d(tQ(0,1),tQ(1,1),tQ(1,2),tQ(1,3));
    Q.row(2) = RowVector4d(tQ(0,2),tQ(1,2),tQ(2,2),tQ(2,3));
    Q.row(3) = RowVector4d(0.,0.,0.,1.0);
    //        RowVector3d target;
    
    double cost;
    RowVector3d vavg = (V.row(v1) + V.row(v2))/2.0;
//    std::cout<<"coming here "<<endl;
   
    if (Q.determinant() != 0 )
    {
//        std::cout<<"starts here 1"<<std::endl;
//        std::cout<<"q.determinant is zero for "<<v1<<" , "<<v2<<endl;
        Vector4d vt = Q.inverse() * Vector4d(0., 0., 0., 1.);
//        cout<<"\n"<<V.row(v1)<<" \n"<< V.row(v2)<<"  \n\n"<<vt.transpose()<<endl;
        target = RowVector3d(vt(0),vt(1),vt(2));
        cost = (double) (vt.transpose() * Q * vt);
//        cout<<"from here1"<<endl;
//        std::cout<<"ends here1 "<<std::endl;
    }
    else
    {
        
//        std::cout<<"q.determinant is not zero "<<v1<<" , "<<v2<<endl;
//        std::cout<<"starts here2 "<<endl;
        double v1_cost = (double) (RowVector4d(V(v1,0),V(v1,1),V(v1,2),1.) * Q * RowVector4d(V(v1,0),V(v1,1),V(v1,2),1.).transpose());
        double v2_cost = (double) (RowVector4d(V(v2,0),V(v2,1),V(v2,2),1.) * Q * RowVector4d(V(v2,0),V(v2,1),V(v2,2),1.).transpose());
        double vavg_cost = (double) (RowVector4d(vavg(0),vavg(1),vavg(2),1.) * Q * RowVector4d(vavg(0),vavg(1),vavg(2),1.).transpose());
        if (v1_cost <= v2_cost && v1_cost<=vavg_cost)
        {
//            cout<<"from here2"<<endl;
            cost = v1_cost;
            target = RowVector3d(V(v1,0),V(v1,1),V(v1,2));
        } else if(v2_cost <= v1_cost && v2_cost<=vavg_cost)
        {
//            cout<<"from here3"<<endl;
            cost = v2_cost;
            target = RowVector3d(V(v2,0),V(v2,1),V(v2,2));
        } else
        {
//            cout<<"from here4"<<endl;
            cost = vavg_cost;
            target = vavg;
        }
//        std::cout<<"ends here 2"<<endl;
    }
        if(vector_contains(removed_indices,v1) )
        {
//            std::cout<<"before "<<cost<<" and after ";
            cost+= 1e9;//10.0+100.0;
//            std::cout<<"cost "<<cost<<endl;;
        }
    else if(vector_contains(removed_indices,v2) )
        {
//            std::cout<<"before "<<cost<<" and after ";
            cost+= 1e9;//10.0+100.0;
//            std::cout<<"cost "<<cost<<endl;;
        }
//    if((distance(target,V.row(v1))+distance(target,V.row(v2))) > 1.5* (distance(vavg,V.row(v1))+distance(vavg,V.row(v2))))
//    {
//        cout<<"going inside: "<<target<<endl;
//
//        target = target.normalized();
//    }
    return cost;
}
void compute_costof_each_pair()
{
    //    costof_pair.resize(vertex_pairs.size());
//    cost.clear();
//    cost.resize(vertex_pairs.size());
//    target.clear();
//    target.resize(vertex_pairs.size());
    for(int i=0; i< vp.size(); i++)
    {
        vp[i].cost = computecost(vp[i].p.a, vp[i].p.b, vp[i].target);
    }
}
void updateF(int i, MatrixXi &Face ){
    MatrixXi new_F(Face.rows()-1,Face.cols());
    int x = 0;
    for(int k=0; k<Face.rows(); k++) {
        if(k != i)
        {
            new_F.row(x++) = Face.row(k);
        }
    }
//    std::cout << "\nbefore " << Face.rows() <<"  ";
    Face = new_F;
//    std::cout<<"index "<<i <<" after   "<<Face.rows()<< std::endl;
}
void updatetempF(int i, MatrixXi &Face ){
    MatrixXi new_F(Face.rows()-1,Face.cols());
    int x = 0;
    for(int k=0; k<Face.rows(); k++) {
        if(k != i)
        {
            new_F.row(x++) = Face.row(k);
        }
    }
//    std::cout << "\nbefore " << Face.rows() <<"  ";
    Face = new_F;
//    std::cout<<"index "<<i <<" after   "<<Face.rows()<< std::endl;
}
int updateF(Pair ab)
{
    MatrixXd tempV = V;
    MatrixXi tempF = F;
    int survived = ab.a;
    int deleted  = ab.b;
    for(int i=0 ; i< tempF.rows(); i++)
    {
        for(int j=0 ; j< tempF.cols(); j++)
        {
            if(tempF(i,j)==deleted)
            {
                tempF(i,j)=survived;
                //                std::cout<<"going 1st ("<<i<<","<<j<<") "<<F.row(i)<<" And "<<index<<std::endl;
            }
        }
    }
    for(int i=0 ; i< tempF.rows(); i++)
    {
        for(int j=0 ; j< tempF.cols(); j++)
        {
//            std::cout<<"face dosms "<<F.row(i)<<std::endl;
            if(tempF(i,j)> deleted)
            {
                
//                std::cout<<"going 2nd ("<<i<<","<<j<<") "<<F.row(i)<<std::endl;
                tempF(i,j)= tempF(i,j)-1;
//                std::cout<<"changed to "<<F.row(i)<<std::endl;
            }
            
        }
    }
    int flag =0;
    for(int i=0 ; i< tempF.rows(); i++)
    {
        if(tempF(i,0)==tempF(i,1)||tempF(i,2)==tempF(i,1)||tempF(i,0)==tempF(i,2))
        {
//            std::cout<<"\n\nWRONG!!!"<< i<<"\n\n"<<std::endl;
//            std::cout<<tempF.row(i)<<endl;
//            if(flag<3)
                updateF(i,tempF);
            flag++;
        }
    }
    if(flag < 3)
    {
        F=tempF;
    }
    return flag;
}

int deleteRow(Pair ab, MatrixXd &VV){
    int index = ab.b;
    MatrixXd tempV(VV.rows()-1,VV.cols());
    int j = 0;
    for(int i=0 ; i< VV.rows(); i++)
    {
        if(i!= index)
        {
            tempV.row(j++)= VV.row(i);
        }
        //        else
        //            std::cout<<"hii "<<index;
    }
    
    int flag = updateF(ab);
    if(flag<3)
        VV= tempV;
   
//    std::cout<<"\n\n faces updateddddd\n\n"<<std::endl;
//    std::cout<<"\n updated V is\n"<<V<<"\n\n before V is\n"<<V_new<<std::endl;
//    std::cout<<"\n updated F is\n"<<F<<"\n\n before F is\n"<<F_new<<std::endl;
    return flag;
}
void update_vpairs(int index)
                       {
    int survived = vp[index].p.a;
    int deleted = vp[index].p.b;
    vp.erase(vp.begin() + index);
//    target.erase(target.begin() + index);
//    cost.erase(cost.begin() + index);
    //std::cout<<" size is  "<<vertex_pairs[15].a<<" "<<vertex_pairs[15].b<<endl;
    for(int i=0; i<vp.size(); i++)
    {
        //        std::cout<<" started "<<i<<endl;
        if(vp[i].p.a==deleted)
        {
//            cout<<"comes here"<<endl;
            //            std::cout<<" started going 1 "<<i<<endl;
            vp[i].p.a=survived;
            //            if(vp[i].p.a>deleted)
            //            {
            //                //                std::cout<<" started going 1.1 "<<i<<endl;
            //                cost[i]=computecost(vertex_pairs[i].a-1, vertex_pairs[i].b-1, target[i]);
            //                //                std::cout<<" ending 1.1 "<<i<<endl;
            //            }
        }
         if(vp[i].p.b==deleted)
        {
//            cout<<"comes here3"<<endl;
            vp[i].p.b=survived;
            if(vp[i].p.a==vp[i].p.b)
            {
                
                vp.erase(vp.begin() + i);
            }
            else
            {
//                cout<<"comes here18"<<endl;
                vp[i].cost = computecost(vp[i].p.a,vp[i].p.b,vp[i].target);
//                cout<<"comes here22"<<endl;
            }
        }
    }
    for(int i=0; i<vp.size(); i++)
    {
        if(vp[i].p.a>deleted)
        {
//            cout<<"comes here6"<<endl;
            vp[i].p.a=vp[i].p.a-1;
        }
        if(vp[i].p.b>deleted)
        {
//            cout<<"comes here7"<<endl;
            vp[i].p.b=vp[i].p.b-1;
        }
    }
    for(int i=0; i<vp.size(); i++)
    {
        if(vp[i].p.a==vp[i].p.b)
        {
//            cout<<"comes here8"<<endl;
            vp.erase(vp.begin() + i);
        }
        else{
            if(vp[i].p.a==survived || vp[i].p.b==survived){
                vp[i].cost= computecost(vp[i].p.a,vp[i].p.b,vp[i].target);
//                cout<<"comes here9"<<endl;
            }
        }
    }
//        else if(vertex_pairs[i].b>=deleted)
//        {
//            //                std::cout<<" started going 1.2 "<<i<<endl;
//            //                std::cout<<" values "<<vertex_pairs[i].a<<" "<<vertex_pairs[i].b-1<<endl;
//            cost[i]=computecost(vertex_pairs[i].a, vertex_pairs[i].b-1, target[i]);
//            //                std::cout<<" ending 1.2 "<<i<<endl;
//        }
//        else
//        {
//            //                std::cout<<" started going 1.3 "<<i<<endl;
//            cost[i]= computecost(vertex_pairs[i].a, vertex_pairs[i].b, target[i]);
//            //                std::cout<<" ending 1.3 "<<i<<endl;
//        }
//        //            std::cout<<" target1 "<< target[i]<<endl;
//        //            std::cout<<" ending 1"<<endl;
//
//
//        vertex_pairs[i].b=survived;
//        //            std::cout<<" started going 2 "<<i<<endl;
//        if(vertex_pairs[i].a>=deleted && vertex_pairs[i].b>=deleted)
//        {
//            //                std::cout<<" started going 2.1 "<<i<<endl;
//            cost[i]=computecost(vertex_pairs[i].a-1, vertex_pairs[i].b-1, target[i]);
//            //                std::cout<<" ending 2.1 "<<i<<endl;
//        }
//        else if(vertex_pairs[i].b>=deleted)
//        {
//            //                std::cout<<" started going 2.2 "<<i<<endl;
//            cost[i]=computecost(vertex_pairs[i].a, vertex_pairs[i].b-1, target[i]);
//            //                std::cout<<" ending 2.2 "<<i<<endl;
//        }
//        else if(vertex_pairs[i].a>=deleted)
//        {
//            //                std::cout<<" started going 2.2 "<<i<<endl;
//            cost[i]=computecost(vertex_pairs[i].a-1, vertex_pairs[i].b, target[i]);
//            //                std::cout<<" ending 2.2 "<<i<<endl;
//        }
//        else
//        {
//            //                std::cout<<" started going 2.3 "<<i<<endl;
//            cost[i]= computecost(vertex_pairs[i].a, vertex_pairs[i].b, target[i]);
//            //                std::cout<<" ending 2.3 "<<i<<endl;
//        }
        //            std::cout<<" ending 2"<<endl;
        //            std::cout<<" target2 "<< target[i]<<endl;
        
        
//            vp[i].p.b = vp[i].p.b -1;
            //            std::cout<<" going 3 "<<endl;
       
            //            std::cout<<" going 31 "<<endl;
        
        //        std::cout<<" coming till end "<<i<<endl;
    
//    std::cout<<" coming till end 1"<<endl;
}

bool isTriangleOrientationValid(Pair vp,int index) {
    

    
    if(to_find_boundary_edges[Edge(vp.a,vp.b)].size()<2)
    {
        std::cout << "isTri comes here" << std::endl;
        return false;
    }
    else if(to_find_boundary_edges[Edge(vp.a,vp.b)].size()==2 )
    {
        
            int f1= to_find_boundary_edges[Edge(vp.a,vp.b)][0];
            int f2= to_find_boundary_edges[Edge(vp.a,vp.b)][1];
            RowVector3d edge_ab1 = V.row(F(f1,1))-V.row(F(f1,0));
            RowVector3d edge_ac1 = V.row(F(f1,2))-V.row(F(f1,0));
            RowVector3d n1 = edge_ab1.cross(edge_ac1).normalized();
            
            RowVector3d edge_ab2 = V.row(F(f2,1))-V.row(F(f2,0));
            RowVector3d edge_ac2 = V.row(F(f2,2))-V.row(F(f2,0));
            RowVector3d n2 = edge_ab2.cross(edge_ac2).normalized();
            double angle = acos(n1.dot(n2) / (n1.norm() * n2.norm()));
//        cout<<"angle = "<<angle<<" "<<M_PI/2.0<<endl;
        if(angle > M_PI/2.0)
        {
            std::cout << "istri 2 comes here" << std::endl;
            return false;
        }
        int v1= vp.a;
        int v2= vp.b;
        std::vector<int> connected_faces;
        for (int i = 0; i < F.rows(); ++i) {
            if (F(i, 0) == v1 || F(i, 1) == v1 || F(i, 2) == v1 ||
                F(i, 0) == v2 || F(i, 1) == v2 || F(i, 2) == v2) {
                connected_faces.push_back(i);
            }
        }

        // Average the normals of the connected triangles
    Vector3d n3=Vector3d(0.,0.,0.);
        for (int i : connected_faces) {
            Eigen::Vector3d v0 = V.row(F(i, 0)).transpose();
            Eigen::Vector3d v1 = V.row(F(i, 1)).transpose();
            Eigen::Vector3d v2 = V.row(F(i, 2)).transpose();
            Eigen::Vector3d n = (v1 - v0).cross(v2 - v0).normalized();
            n3 += n;
        }
        n3.normalize();
    double threshold = 0.9; // Set threshold to 0.9

    // Calculate the dot product of the two normals before and after the collapse
    double dot_before = n1.dot(n2);
    double dot_after = n3.dot(n2);
//        cout<<"dot "<<dot_before<<" - "<< dot_after<< endl;// > threshold
        if (dot_before - dot_after > threshold) {
            cout<<"dot "<<dot_before<<" - "<< dot_after<< endl;// > threshold
            return false;
        }
    }
    
    return true;
//    Vector3d v1=V.row(vp.a) - target[index];
//    Vector3d v2= target[index] - V.row(vp.b);
//    v1=v1.normalized();
//    v2=v2.normalized();
//    double dot = v1.dot(v2);
//    double angle = std::acos(dot / (v1.norm() * v2.norm())) ;
//    std::cout<<" v's = "<< V.row(vp.a)<<" " << V.row(vp.b)<<" "<<target[index]<<" " <<"angle "<<angle<<endl;
//    if(distance(V.row(vp.a),target[index]) + distance(V.row(vp.a),target[index]) > 1.9* distance(V.row(vp.a),V.row(vp.b)))
//        return false;
//    return true;
        
//    igl::vertex_triangle_adjacency(V, F, vertex_face_adjacency, VFi);
//    igl::per_face_normals(V, F, face_normals);
//
//    // Get the set of vertices adjacent to the vertex pair
//       vector<int> adjacent_vertices;
//       for (auto f : vertex_face_adjacency[vp.a])
//           for (int v=0 ; v<3; v++)//; F.col)
//               if (F(f,v)!= vp.a && F(f,v)!=  vp.b)
//                   adjacent_vertices.push_back(F(f,v));
//       for (auto f : vertex_face_adjacency[vp.b])
//           for (int v=0 ; v<3; v++)
//               if (F(f,v)!= vp.a && F(f,v)!=  vp.b)
//                   adjacent_vertices.push_back(F(f,v));
//
//       // Get the normals of the adjacent faces before the collapse
//       vector<Vector3d> adjacent_normals_before;
//       for (int v : adjacent_vertices)
//           for (auto f : vertex_face_adjacency[v])
//               if (F(f,0)==(vp.a) || F(f,1)==(vp.a) || F(f,2)==(vp.a) || F(f,0)==(vp.b) || F(f,1)==(vp.b) || F(f,2)==(vp.b) )
//                   adjacent_normals_before.push_back((face_normals.row(f)));
//
//       // Collapse the vertex pair
////       VertexData new_vertex = collapseVertex(vp);
//
//       // Get the normals of the adjacent faces after the collapse
//       vector<Vector3d> adjacent_normals_after;
//       for (int v : adjacent_vertices)
//           for (auto f : vertex_face_adjacency[v])
//               if (F(f,0)==(vp.a) || F(f,1)==(vp.a) || F(f,2)==(vp.a))
//                   adjacent_normals_after.push_back((face_normals.row(f)));
//
//       // Check that the dot product of each pair of normals is positive
//       for (auto n1 : adjacent_normals_before)
//           for (auto n2 : adjacent_normals_after)
//               if (n1.dot(n2) < 0)
//                   return false;
//
//       // Collapse is valid
//       return true;

//    Eigen::RowVector3d n_new = RowVector3d((V.row(v2) - V.row(v1)).cross(RowVector3d(V.row(v3) - V.row(v1)))).normalized();
//    Eigen::RowVector3d n_avg = Eigen::RowVector3d(0.,0.,0.);
//    for (int i = 0; i < F.rows(); i++) {
//        if (F(i, 0) == v1 || F(i, 1) == v1 || F(i, 2) == v1 ||
//            F(i, 0) == v2 || F(i, 1) == v2 || F(i, 2) == v2 ||
//            F(i, 0) == v3 || F(i, 1) == v3 || F(i, 2) == v3) {
//            Eigen::RowVector3d e1 = V.row(F(i, 1)) - V.row(F(i, 0));
//            Eigen::RowVector3d e2 = V.row(F(i, 2)) - V.row(F(i, 0));
//            n_avg += e1.cross(e2);
//        }
//    }
//    n_avg.normalize();
//    return n_avg.dot(n_new) >= 0;
}
bool hasSelfIntersections(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    // Iterate over all pairs of faces
    for (int i = 0; i < F.rows(); i++) {
        for (int j = i+1; j < F.rows(); j++) {
            
            // If the faces share a vertex, they cannot intersect
            bool shareVertex = false;
            for (int k = 0; k < 3; k++) {
                if (F(i,k) == F(j,0) || F(i,k) == F(j,1) || F(i,k) == F(j,2)) {
                    shareVertex = true;
                    break;
                }
            }
            if (shareVertex) {
                continue;
            }
            
            // Calculate the intersection point between the two faces
            Eigen::Vector3d v1 = V.row(F(i,0));
            Eigen::Vector3d v2 = V.row(F(i,1));
            Eigen::Vector3d v3 = V.row(F(i,2));
            Eigen::Vector3d w1 = V.row(F(j,0));
            Eigen::Vector3d w2 = V.row(F(j,1));
            Eigen::Vector3d w3 = V.row(F(j,2));
            Eigen::Vector3d n1 = (v2 - v1).cross(v3 - v1);
            Eigen::Vector3d n2 = (w2 - w1).cross(w3 - w1);
            Eigen::Vector3d r = v1 - w1;
            double den = n1.dot(n2.cross(-n1));
            if (den == 0) {
                continue;
            }
            double t = n2.dot(r.cross(n1)) / den;
            if (t < 0 || t > 1) {
                continue;
            }
            double u = -n1.dot(r.cross(n2)) / den;
            if (u < 0 || u > 1) {
                continue;
            }
            double v = 1 - t - u;
            if (v < 0 || v > 1) {
                continue;
            }
            std::cout << "intersections part" << std::endl;
            // The faces intersect
            return true;
        }
    }
    
    // The mesh has no self-intersections
    return false;
}

bool isthisvalid(Pair ab, int index){
    std::vector<Edge> temp_edge_map_;
    std::map< Edge , std::vector< int > > temp_to_find_boundary_edges;
    
        MatrixXd tempV = V;
        MatrixXi tempF = F;
        vector<VertexPairs> tempVP = vp;
    
    int survived = ab.a;
    int deleted  = ab.b;
    for(int i=0 ; i< tempF.rows(); i++)
    {
        for(int j=0 ; j< tempF.cols(); j++)
        {
            if(tempF(i,j)==deleted)
            {
                tempF(i,j)=survived;
            }
        }
    }
    for(int i=0 ; i< tempF.rows(); i++)
    {
        for(int j=0 ; j< tempF.cols(); j++)
        {

            if(tempF(i,j)> deleted)
            {
                tempF(i,j)= tempF(i,j)-1;
            }
            
        }
    }
    for(int i=0 ; i< tempF.rows(); i++)
    {
        if(tempF(i,0)==tempF(i,1)||tempF(i,2)==tempF(i,1)||tempF(i,0)==tempF(i,2))
        {
            std::cout<<"going here"<<endl;
            updatetempF(i,tempF);
        }
    }
    if(hasSelfIntersections(tempV,tempF))
    {
        std::cout << "intersec" << std::endl;
        return false;
    }
    
    for(int i=0; i<tempF.rows(); i++)
    {
        Edge edge1= Edge(tempF(i,0) , tempF(i,1));
        Edge edge2= Edge(tempF(i,1) , tempF(i,2));
        Edge edge3= Edge(tempF(i,2) , tempF(i,0));
        temp_to_find_boundary_edges[edge1].push_back(i);
        temp_to_find_boundary_edges[edge2].push_back(i);
        temp_to_find_boundary_edges[edge3].push_back(i);
//        if(!edge_contains(temp_edge_map_, edge1))
//            temp_edge_map_.push_back(edge1);
//        if(!edge_contains(temp_edge_map_, edge2))
//            temp_edge_map_.push_back(edge2);
//        if(!edge_contains(temp_edge_map_, edge3))
//            temp_edge_map_.push_back(edge3);
    }
    
    for(int i=0; i< tempF.rows(); i++)
    {
        if( temp_to_find_boundary_edges[Edge(tempF(i,0),tempF(i,1))].size()<2 || temp_to_find_boundary_edges[Edge(tempF(i,1),tempF(i,2))].size()<2 || temp_to_find_boundary_edges[Edge(tempF(i,2),tempF(i,0))].size()<2)
        {
            std::cout<<" isthisvalid comes here 1 "<<endl;
            return false;
        }
        else if (temp_to_find_boundary_edges[Edge(tempF(i,0),tempF(i,1))].size()>2 ||
                 temp_to_find_boundary_edges[Edge(tempF(i,1),tempF(i,2))].size()>2 ||
                 temp_to_find_boundary_edges[Edge(tempF(i,2),tempF(i,0))].size()>2)
        {
            std::cout<<" isthisvalid comes here 2 "<<endl;
            return false;
        }
//        else if (temp_to_find_boundary_edges[Edge(tempF(i,0),tempF(i,1))].size()==2 ||
//                 temp_to_find_boundary_edges[Edge(tempF(i,1),tempF(i,2))].size()==2 ||
//                 temp_to_find_boundary_edges[Edge(tempF(i,2),tempF(i,0))].size()==2)
//        {
////            std::cout<<" yessyessyess "<<endl;
//        }
    }
    
    
    return true;
}

//bool anothervalidcheck(Pair vp, int index){
//    int v1 = vp.a;
//    int v2 = vp.b;
//    // Find all triangles that are connected to v1 or v2
//        std::vector<int> connected_faces;
//        for (int i = 0; i < F.rows(); ++i) {
//            if (F(i, 0) == v1 || F(i, 1) == v1 || F(i, 2) == v1 ||
//                F(i, 0) == v2 || F(i, 1) == v2 || F(i, 2) == v2) {
//                connected_faces.push_back(i);
//            }
//        }
//
//        // Average the normals of the connected triangles
//    Vector3d n3=Vector3d(0.,0.,0.);
//        for (int i : connected_faces) {
//            Eigen::Vector3d v0 = V.row(F(i, 0)).transpose();
//            Eigen::Vector3d v1 = V.row(F(i, 1)).transpose();
//            Eigen::Vector3d v2 = V.row(F(i, 2)).transpose();
//            Eigen::Vector3d n = (v1 - v0).cross(v2 - v0).normalized();
//            n3 += n;
//        }
//        n3.normalize();
//    double threshold = 0.9; // Set threshold to 0.9
//
//    // Calculate normals of the two neighboring faces before the collapse
//    Eigen::Vector3d n1 = face_normals.row(f1).normalized();
//    Eigen::Vector3d n2 = face_normals.row(f2).normalized();
//    // Calculate the dot product of the two normals before and after the collapse
//    double dot_before = n1.dot(n2);
//    double dot_after = n3.dot(n2);
//
//    // Check if the change in the dot product is above the threshold
//    if (dot_before - dot_after > threshold) {
//        // Collapse should not be performed
//    }
//    else {
//        // Collapse can be performed
//    }
//}

void find_mincost(){
    auto result = std::min_element(vp.begin(), vp.end(), comp);
    
    int index = std::distance(vp.begin(), result);
    //std::cout<<vertex_pairs[index].a<<" "<<vertex_pairs[index].b<<"<--index "<<target[index]<<" <- "<<V.row(vertex_pairs[index].a)<<" and "<<V.row(vertex_pairs[index].b)<<" cost is "<<cost[index]<<endl;
    
//        removed_indices.push_back(vp[index].p.a);
    
        if(vertex_is_collapsable[vp[index].p.b]&&vertex_is_collapsable[vp[index].p.a])
        {
            
            if(isTriangleOrientationValid(vp[index].p, index) && isthisvalid(vp[index].p,index))
            {
               
                int flag= deleteRow(vp[index].p, V);
                if(flag <3)
                {
//                    cout<<"comes here4"<<endl;
                    create_edgemap();
                    
                    vertex_is_collapsable_check();
                    
                    V.row(vp[index].p.a) = vp[index].target;
                    removed_indices.push_back(vp[index].p.a);
                    //        std::cout<<"going"<<endl;
                    update_vpairs(index);
//                    cout<<"comes here5"<<endl;
                }
//                else
//                    vp[index].cost += flag*(10.0);
//                std::cout<<"end"<<std::endl;
            }
            else
            {
//                std::cout<<"FALSE"<<endl;
                vp[index].cost += 1e9;
            }
        }
        else
        {
                        std::cout<<"split edge"<<std::endl;
            vp[index].cost +=1e9;
//            vp[index].target = (V.row(vp[index].p.a) + V.row(vp[index].p.b) ) *0.5;
//            int flag= deleteRow(vp[index].p, V);
//            if(flag <3)
//            {
//                create_edgemap();
//                vertex_is_collapsable_check();
//                V.row(vp[index].p.a) = vp[index].target;
//                removed_indices.push_back(vp[index].p.a);
//                //        std::cout<<"going"<<endl;
//                update_vpairs(index);
//            }
//            else
//                vp[index].cost += flag*(10.0);
//            std::cout<<"end"<<std::endl
        }
            
}

// function to find the smallest error
//in which update pairs will be called and update costof each pair also will be called
//And remove vertex function should be created along with update faces function;

/* find the smallest error and new position
 update the first vertex's position and remove the second vertex and update faces
 //for now may be initialize pairs again or else update the pairs and cost of the pairs and remove duplicate pairs update edge map with new connections.
 */
/*void update_edge_map(int v1, int v2){
 //v2 is removed, v1 has v2's connections
 }
 
 void remove_duplicate_pairs(){
 //    std::unordered_set<Pair> seen;
 //        int i = 0;
 //        for (const auto &p : vertex_pairs) {
 //            if (seen.count(p) == 0) {
 //                seen.insert(p);
 //                vertex_pairs[i++] = p;
 //            }
 //        }
 //        vertex_pairs.resize(i);
 }*/

/*void update_vertexpairs(int v1, int v2){
 
 todo:replace v2 with v1 in all the pairs
 for (auto &p : vertex_pairs) {
 if (p.a == v2) {
 p.a = v1;
 }
 if (p.b == v2) {
 p.b = v1;
 }
 }
 remove_duplicate_pairs();
 }*/

void init(){
    
    create_edgemap();
    vertex_is_collapsable_check();
    initialize_vertex_pairs();
    compute_costof_each_pair();
}
int  abc;
int target_vertex_count;
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier){
    switch(key){
        case '1':
        {
            find_mincost();
            viewer.data(0).clear();
            viewer.data(0).set_mesh(V, F);
//            std::cout << "old Vertices: \n" << V_new.rows() << std::endl;
            std::cout << "new Vertices: \n" << V.rows() << std::endl;
        }
            return true;
        case '2':
        {
            abc=V.rows();
            while(V.rows() >abc-1 )
                find_mincost();
            viewer.data(0).clear();
            viewer.data(0).set_mesh(V, F);
            //            std::cout << "old Vertices: \n" << V_new.rows() << std::endl;
            std::cout << "new Vertices: \n" << V.rows() << std::endl;
            std::cout<<"hi"<<endl;
            return true;
        }
        case '3':
        {
            while (V.rows()  > V_new.rows()-target_vertex_count) {
                find_mincost();
            }
            viewer.data(0).clear();
            viewer.data(0).set_mesh(V, F);
//            std::cout << "old Vertices: \n" << V_new.rows() << std::endl;
            std::cout << "new Vertices: \n" << V.rows() << std::endl;
        }
        default: break;
    }
    return false;
}







// ------------ main program ----------------
int main(int argc, char *argv[])
{
    string f = "../data/bunny.off";
    threshold = 0.0;
    target_vertex_count=100;
    if (argc>=2) {
      string w = argv[1];
      f = "../data/" + w + ".off";
    }
    if (argc>=3) {
        target_vertex_count = atoi(argv[2]);
    }
    if (argc>=4) {
        threshold = atoi(argv[3]);
    }
    
            igl::readOFF(f, V_new, F_new); // Load an input mesh in OFF format
    //    V.resize(7,3);
    //    V << 0.0, 0.0, 0.0,
    //    1.0, 0.0, 0.0,
    //    2.0, 0.0, 0.0,
    //    2.0, 1.0, 0.0,
    //    1.0, 1.0, 0.0,
    //    0.0, 1.0, 0.0,
    //    1.0, 2.0, 0.0;
    //    F.resize(6,3);
    //    F << 0, 1, 5,
    //    1,2,4,
    //    1,4,5,
    //    2,3,4,
    //    4,3,6,
    //    4,6,5;

//    V_new.resize(8,3);
//    V_new << 0.0, 0.0, 0.0,
//    1.0, 0.0, 0.0,
//    1.0, 1.0, 0.0,
//    0.0, 1.0, 0.0,
//    0.0, 0.0, 1.0,
//    1.0, 0.0, 1.0,
//    1.0, 1.0, 1.0,
//    0.0, 1.0, 1.0;
//
//    F_new.resize(12,3);
//    F_new << 0, 1, 2,
//    0, 2, 3,
//    1, 5, 6,
//    1, 6, 2,
//    5, 4, 7,
//    5, 7, 6,
//    4, 0, 3,
//    4, 3, 7,
//    3, 2, 6,
//    3, 6, 7,
//    0, 4, 5,
//    0, 5, 1;

    
    //    0, 1, 2,
    //         2, 3, 0,
    //         4, 5, 6,
    //         6, 7, 4,
    //         0, 3, 7,
    //         7, 4, 0,
    //         1, 5, 6,
    //         6, 2, 1,
    //         0, 5, 1,
    //         1, 2, 0,
    //         3, 6, 7,
    //         7, 2, 3;
    
    //    Eigen::MatrixXd dV;
    //    Eigen::MatrixXi dF;
    //    Eigen::VectorXi SVI,SVJ,SF;
    V= V_new;

    F=F_new;
    
    //  print the number of mesh elements
    //    std::cout << "Vertices: " << V.rows() <<"\n coming here "<<calculate_quadric(4)<<std::endl;
    //    create_edgemap();
    
    
    std::cout<<"Target vertices count = "<<target_vertex_count<<std::endl;
    int count=0;
    init();
//    std::cout<<"..."<<std::endl;
//    while (V.rows()  > target_vertex_count) {
//        if(count>1)
//        {
//
//        }
//        std::cout<<"...";
//        find_mincost();
//        count++;
        
//    }
    
    

//    std::ofstream file("vertices.txt");
//        if (file.is_open())
//        {
//            file << V << std::endl;
//            file.close();
//        }
//        else
//        {
//            std::cerr << "Unable to open file!" << std::endl;
//        }
    std::cout << "old Vertices: \n" << V_new.rows() << std::endl;
    std::cout << "new Vertices: \n" << V.rows() << std::endl;
    //    std::cout << "old faces: \n" << F_new << std::endl;
    //    std::cout << "new faces: \n" << F << std::endl;
    //    std::cout << "new Faces:    " << F << std::endl;
    igl::opengl::glfw::Viewer viewer; // create the 3d viewer
    viewer.callback_key_down = &key_down;
    viewer.data(0).set_mesh(V, F); // load a face-based representation of the input 3d shape
    //    draw_bounding_box(viewer, V); // draw the boundaing box (red edges and vertices)
    viewer.launch(); // run the editor
}
