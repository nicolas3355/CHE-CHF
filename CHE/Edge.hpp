/**
* @file    Edge.hpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief  (Edge Class).
*/
//--------------------------------------------------//



#ifndef _Edge_HPP_

#define _Edge_HPP_



#include <iostream>



using namespace std;

//--------------------------------------------------//

/** Edge class for CHF data structure 

  * \brief Edge class*/

class Edge

//--------------------------------------------------//

{

//-- Edge protected data --//

protected:

  /** Edge Vertices

    * \brief first vertice of the edge*/

	int  _a;    

  /** Edge Vertices

    * \brief first vertice of the edge*/

	int  _b;    

  /** Edge Half-Edge

    * \brief first vertice of the edge*/

	int _h;    



//-- Edge constructors.--// 

public:

  /** \brief Default constructor.*/

	Edge():_a(-1), _b(-1), _h(-1){}



  /** \brief First constructor

    * \param int  a  Edge _a

    * \param int  b  Edge _b

    * \param int  h  Edge _h */

	Edge( int  a, int  b, int h ):_a(a), _b(b), _h(h){}

  

  /** \brief Copy constructor

    * \param Edge& _e*/

	Edge( const Edge& e ):_a(e.a()), _b(e.b()), _h(e.h()){}



  /** \brief Destructor.*/

  ~Edge(){}



//-- Edge methods.--//

public:

	/** \brief Access to the vertex a of the edge*/

  inline const int   a() const{ return _a ; }



	/** \brief Access to the vertex b of the edge*/

  inline const int   b() const{ return _b ; }



	/** \brief Access to the a half-edge of the edge*/

  inline const int   h() const{ return _h ; }



public:

	/** \brief Sets the vertex a of the edge

	  * \param const int a */

	inline const void  set_a (const int a ) { _a=a; }

	

	/** \brief Sets the vertex b of the edge

    * \param const int b */

  inline const void  set_b (const int b ) { _b=b; }

	

	/** \brief Sets the a half-edge of the edge

	  * \param const int h */

  inline const void  set_h (const int h ) { _h=h  ; }



//-- Edge i/o operators --//

public:

  /** \brief Read operator for Edge objects 
    * \param istream& s
    * \param Edge& e */
  friend istream& operator >> (istream& s,  Edge& e)

  {

    int a, b, h;



    s >> a >> b >> h  ;



    e.set_a(a); e.set_b(b); e.set_h(h);

	  

	return s;

  }

  /** \brief Write operator for Edge objects 
    * \param ostream& s
    * \param Edge& e */
  friend ostream& operator << (ostream& s, const Edge& e)

  {

    s << "  Verts: ( " << e.a() << " , " << e.b() <<  " ) " << endl ;

    s << "     EH:   " << e.h() << endl ;

	return s;

  }



//-- Edge bool operators --//

public: 
  /** \brief Equal operator for Edge objects 
    * \param Edge& e */
   bool operator == (const Edge& e){ return (a()==e.a() && b()==e.b() && h()==e.h()) ;}

  
  /** \brief Equal operator for Edge objects 
    * \param Edge& e */
   bool operator != (const Edge& e){ return (a()!=e.a() && b()!=e.b() && h()!=e.h()) ;}



};

#endif
