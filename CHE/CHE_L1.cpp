/**
* @file    CHE_L1.cpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmidia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 1)
*/


#include <ctime>

#include <stack>

#include "CHE_L1.hpp"



using namespace std;

//--------------------------------------------------//
vector<Vid> CHE_L1::R_00(const Vid v)
//--------------------------------------------------//
/** Computes the vertices in the star of a given vertex.*/
{
  vector<Vid> star;

  if( !v_valid(v) ) { 
    cout << "CHE_L1::Vertex Star ERROR: invalid vertex id" << endl;
    star.push_back(INV); 
    return star; 
  }

  HEid h = 0, h0 = 0;
  
  //Gets the first incident half--edge
  for(HEid i=0; i<3*ntrig(); ++i)
    if( V(i) == v ) { h0 = h = i; break; }
  
  //Walks through the first wise.  
  do 
    star.push_back( V(next(h)) );
  while( (O(h) != -1) && ( (h = next(O(h))) != h0) );

  //Walks through the second wise.
  if( h != h0 ){
   h = h0;
   do
     star.push_back( V(prev(h)) );
   while( (h = O(prev(h))) != -1 );
  }
  return star;
}
//--------------------------------------------------//
vector<TRid> CHE_L1::R_02(const Vid v)
//--------------------------------------------------//
/** Computes the triangles in the star of a given vertex.*/
{
  vector<Vid> star;
  
  if( !v_valid(v) ) { 
    cout << "CHE_L1::Vertex Star ERROR: invalid vertex id" << endl;
    star.push_back(INV); 
    return star; 
  }

  HEid h = 0, h0 = 0;
  
  //Gets the first incident half--edge
  for(HEid i=0; i<3*ntrig(); ++i)
    if( V(i) == v ) { h0 = h = i; break; }
  
  //Walks through the first wise.  
  do 
    star.push_back( trig(h) );
  while( (O(h) != -1) && ( (h = next(O(h))) != h0) );

  //Walks through the second wise.
  if( h != h0 ){
   h = h0;
   while( (h = O(prev(h))) != -1 );
      star.push_back( trig(h) );
  }
  return star;
}
//--------------------------------------------------//
vector<Vid> CHE_L1::R_10( const HEid h )
//--------------------------------------------------//
/** Computes the vertices in the star of a given edge.*/
{
  vector<Vid> star;

  if( !he_valid(h) || !he_valid(next(h)) ) 
  { 
    cout << "CHE_L1::Edge Star ERROR: invalid edge id" << endl;
    star.push_back(INV); 
    return star; 
  }

  star.push_back( V(prev(   h  )) );
  star.push_back( V(prev( O(h) )) );

  return star;
}
//--------------------------------------------------//
vector<TRid> CHE_L1::R_12( const HEid h )
//--------------------------------------------------//
/** Computes the triangles in the star of a given edge.*/
{
  vector<Vid> star;

  if( !he_valid(h) || !he_valid(next(h)) ) 
  { 
    cout << "CHE_L1::Edge Star ERROR: invalid edge id" << endl;
    star.push_back(INV); 
    return star; 
  }

  star.push_back( trig( h  ) );
  star.push_back( trig(O(h)) );

  return star;
}
//--------------------------------------------------//
vector<TRid> CHE_L1::R_22(const TRid t)
//--------------------------------------------------//
/** Computes the triangle incidents of a given triangle.*/
{
  vector<TRid> star;

  if( !tr_valid(t) ) 
  { 
    cout << "CHE_L1::Triangle Star ERROR: invalid triangle id" << endl;
    star.push_back(INV); 
    return star; 
  }

  star.push_back( trig(O(3*t)) );
  star.push_back( trig(O(3*t+1)) );
  star.push_back( trig(O(3*t+2)) );

  return star;
}
//--------------------------------------------------//
void CHE_L1::compute_opposites()
//--------------------------------------------------//
/** Computes the opposite of each half-edge.*/
{
  if( _V.size() == 0 && _G.size() == 0 ) return;	

  cout << "Pet_CHE::compute_opposites...  " ;
  
  Vid a, b;

  map<pair<Vid,Vid>,HEid> adjacency ;
  map<pair<Vid,Vid>,HEid>::iterator pos ;

  _O.clear();
  _O.resize( 3*ntrig(), -1 );

  for ( HEid c = 0 ; c < 3*ntrig() ; ++c )
  {
    if( !he_valid(c) ) continue;

    a = V(c) ;  b = V(next(c)) ;

    if( b < a) { Vid tmp = a ; a = b ; b = tmp ; }

    pos = adjacency.find( pair<Vid,Vid>(a,b) ) ;

    if( pos != adjacency.end() )
    {
      set_O(c, pos->second) ;
      set_O(O(c),c) ;

      adjacency.erase( pos ) ;
    }
    else
    {
      adjacency[pair<Vid,Vid>(a,b)] = c ;
    }
  }
  adjacency.clear() ;
  
  cout << " done." << endl;
}
//--------------------------------------------------//
void CHE_L1::compute_connected()
//--------------------------------------------------//
/** Computes the compound of each vertex.*/
{
  if( _V.size() == 0 && _G.size() == 0 ) return;	

  cout << "Pet_CHE::compute_bounds...  " ;

  _C.clear () ;
  _C.resize( nvert(), -1 ) ;

  for( Vid v=0; v< nvert(); ++v)
    set_C( v, v);

  for( TRid j=0; j< ntrig(); ++j)
  {
    if( !tr_valid(j) ) continue;

    Cid b0= get_component( V(3*j) );
    Cid b1= get_component( V(3*j+1) );
    Cid b2= get_component( V(3*j+2) );

    if( b0< 0 || b1<0 || b2<0 ) continue;
    
    if( b0 < b1 )
    {
      if( b0 < b2 )
      { /** b0 min*/
        set_C( b1, b0 ) ;
        set_C( b2, b0 ) ;
      }
      else
      { /** b2 min*/
        set_C( b0, b2 ) ;
        set_C( b1, b2 ) ;
      }
    }
    else
    {
      if( b1 < b2 )
      { /** b1 min*/
        set_C( b0, b1 ) ;
        set_C( b2, b1 ) ;
      }
      else
      { /** b2 min*/
        set_C( b0, b2 ) ;
        set_C( b1, b2 ) ;
      }
    }
  }

  set_nbound(0);

  map<Cid,Cid> m;
  map<Cid,Cid>::iterator it;
  for( Vid v=0; v< nvert(); ++v)
  {
    if( !v_valid(v) ) continue;

    Cid b= get_component( v );
    if(b<0) continue;

    it = m.find( b ) ;
    if( it== m.end() )
      m.insert(it, make_pair(b, _ncomp++) );
  }

  for(Vid v=0; v< nvert(); ++v)
  {
   if( !v_valid(v) ) continue;

    Cid b= get_component( v );
    if(b<0) continue;

    set_C( v, m[b]);        
  }
  m.clear();

  cout <<" "<< ncomp() << " connected compound(s) found." << endl;
}
//--------------------------------------------------//
void CHE_L1::check()
//--------------------------------------------------//
/** Checks the mesh.*/
{
  CHE_L0::check();

  if(3*ntrig() != static_cast<int>(_O.size()) ){ cout << "CHE_L1:: Erro 3*ntrig()!= O.size" << endl; return;}
  if(  nvert() !=  static_cast<int>(_C.size()) ){ cout << "CHE_L1:: Erro   nvert()!= C.size" << endl; return;}

  for(int i=0; i<3*ntrig(); ++i)
  {
    if( O(i)>=0 && O(O(i)) != i )
    {
      cout << "CHE_L1:: Erro O(O(" << i <<")) != "<< i << endl;
      return;
    }
    if( O(i)>=0 && !orient_check(i, O(i)) )
    {
      cout << "CHE_L1:: triangle "<< trig(i) << " bad oriented." << endl;
      return;
    }
  }

  for(int i=0; i<nvert(); ++i)
  {
    if(C(i) > ncomp() )
    {
      cout << "CHE_L1:: C("<< i << ") > nbound." << endl;
      return;
    } 
    if(C(i) < 0 )
    {
      cout << "CHE_L1:: C("<< i << ") < 0." << endl;
      return;
    }
  }
}
//--------------------------------------------------//
void CHE_L1::orient()
//--------------------------------------------------//
/** Orients the mesh.*/
{
  if( _V.size() == 0 && _G.size() == 0 ) return;	

  cout << "Pet_CHE::orient... " ;

  stack<HEid> s;
  vector<bool> visited; 
  visited.resize( ntrig(), false );

  s.push( 0 ); s.push( 1 ); s.push( 2 );
  visited[0] = true;

  while( !s.empty() )
  {
    HEid h = s.top();
    s.pop();
    
    /** Avoid null edges*/
    HEid o = O(h); 
    if( o == -1 ) continue;

    /** Avoid Loops*/
    TRid t = trig(o);
    if( visited[t] ) continue;

    /** Repairs orientation*/
    if( !orient_check(h,o) ) orient_change( t );

    /** Marks as visited*/
    visited[t] = true;

    /** Push half-edges of t*/
    s.push( 3*t ); s.push( 3*t+1 ); s.push( 3*t+2 );
  }
  cout << endl << "Pet_CHE::orient done." << endl;
}
//--------------------------------------------------//
void CHE_L1::orient_change(const TRid t)
//--------------------------------------------------//
/** Changes the orientation of a triangle.*/
{
  cout << endl << "Pet_CHE::orient_change... " << t ;	
  HEid h1= 3*t;
  HEid h2= 3*t+1;
  HEid h3= 3*t+2;

  /** exchange the geometrical vertices */
  Vid v = V(h1);
  set_V( h1, V(h2));
  set_V( h2,   v  );

  /** exchange mates*/
  HEid o2 = O(h2);
  HEid o3 = O(h3);
  set_O( h2, o3 );
  set_O( h3, o2 );
  set_O( o2, h3 );
  set_O( o3, h2 );
}
//--------------------------------------------------//
const bool CHE_L1::orient_check(const HEid h, const HEid o)
//--------------------------------------------------//
/** Checks the orientation between two half-edges.*/
{
  Vid v1, v2, v3, v4;
  
  v1 = V(h);
  v2 = V(next(h));
  
  v3 = V(o);
  v4 = V(next(o));

  return( v1 == v4 && v2 == v3 );
}
//--------------------------------------------------//
Cid CHE_L1::get_component( Cid i )
//--------------------------------------------------//
/** Gets the connected compound. Sets _C*/
{
  Cid b = C(i) ;
  if( b < 0 || b == i ) return b ;

  b = get_component(b) ;
  set_C( i,b );

  return b;
}
//--------------------------------------------------//
void CHE_L1::read_ply( const char* file )
//--------------------------------------------------//
/** Gets a boundary compound. Sets _C*/
{
  // Stores Start && L1 time;
  clock_t start_time = 0,
             L1_time = 0;

  start_time = clock();


	CHE_L0::read_ply( file );
	compute_opposites();
	orient();
	compute_connected();
	compute_normals();


  L1_time = clock();

  //cout << "L1 load time:" << static_cast<double>(L1_time-start_time)/static_cast<double>(CLOCKS_PER_SEC) << endl;

}
//------------------------------------//
