/**
* @file    Test_Pet_CHE.cpp
* @author  Marcos Lage         <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner      <thomas.lewiner@polytechnique.org>
* @author  Helio  Lopes        <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Test)
*/
//--------------------------------------------------//



#include <GL/glui.h>

#include <GL/glut.h>



#include <map>

#include "CHE_L3.hpp"



using namespace std;

//------GLUI---------//

GLUI  *glui;

GLUI_Rotation *view_rot;

GLUI_Translation *trans_xy;

GLUI_Translation *zoom;

GLUI_Checkbox *c;

GLUI_Spinner *simplexid;



float scale= 1.0;

float view_rotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };

float       obj_pos[] = { 0.0, 0.0, 0.0 };



int view_w=720,view_h=650, window_id;



//------Lights-------//

GLfloat light0_ambient [] =	{0.2, 0.2, 0.2, 1.0};

GLfloat light0_diffuse [] =	{.4, .8, 0.7, 1.0};

GLfloat light0_position[] = {10.0, 10.0, 10.0, 0.0};



GLfloat light1_ambient [] =	{0.2, 0.2, 0.2, 1.0};

GLfloat light1_diffuse [] =	{.6, .6, 0.3,1.0};

GLfloat light1_position[] = {-10.0, -10.0, 10.0, 0.0};



//------Material-----//

GLfloat material_shininess[1] = { 50.0 } ;



//------Globals-----//

CHE_L0 ch0;

CHE_L1 ch1;

CHE_L2 ch2;

CHE_L3 ch3;

char file[1024]= ".ply";

vector<int> star; 

int smooth=true, wire=true, points=false, vstar=true, level = 0, vid = 0, dim = 0, nverts = 0;

//----------------------------------------------------------------//

vector<int> test_vstar(int vid, int dim)

//----------------------------------------------------------------//

{

  vector<int> st;



  if( dim == 0 ){

    if(level == 0) st = ch0.R_00(vid);

    if(level == 1) st = ch1.R_00(vid);

    if(level == 2) st = ch2.R_00(vid);

    if(level == 3) st = ch3.R_00(vid);

  }

  else{

    if(level == 0) st = ch0.R_02(vid);

    if(level == 1) st = ch1.R_02(vid);

    if(level == 2) st = ch2.R_02(vid);

    if(level == 3) st = ch3.R_02(vid);

  }



  return st;

}

//----------------------------------------------------------------//

void draw_vstar(int dim)

//----------------------------------------------------------------//

{

  for(int i=0; i< static_cast<int>(star.size()); ++i)

 {

    //--- Vertex ---//

    if( dim == 0 ){

      Vertex v1;



		  if(level == 0) v1 = ch0.G(star[i]);

		  if(level == 1) v1 = ch1.G(star[i]);

		  if(level == 2) v1 = ch2.G(star[i]);

		  if(level == 3) v1 = ch3.G(star[i]);

			

		  glColor3f(1.0,0.3,0.1);

		  glPointSize(2.0);

		  glBegin(GL_POINTS);

			  glVertex3d(v1.x() ,  v1.y() ,  v1.z() );

		  glEnd();

    }



    //--- Faces ---//

    else 

	{

		  Vertex v1, v2, v3;

	

		  if(level == 0) {

			v1 = ch0.G( ch0.V(3*star[i])   );

			v2 = ch0.G( ch0.V(3*star[i]+1) );

			v3 = ch0.G( ch0.V(3*star[i]+2) );

		  }

		  if(level == 1) {

			v1 = ch1.G( ch1.V(3*star[i])   );

			v2 = ch1.G( ch1.V(3*star[i]+1) );

			v3 = ch1.G( ch1.V(3*star[i]+2) );

		  }

		  if(level == 2){

			v1 = ch2.G( ch2.V(3*star[i])   );

			v2 = ch2.G( ch2.V(3*star[i]+1) );

			v3 = ch2.G( ch2.V(3*star[i]+2) );

		  }

		  if(level == 3){

			v1 = ch3.G( ch3.V(3*star[i])   );

			v2 = ch3.G( ch3.V(3*star[i]+1) );

			v3 = ch3.G( ch3.V(3*star[i]+2) );

		  }

	

		  glColor3f(1.0,0.0,0.0);

			  glBegin(GL_TRIANGLES);

				  glNormal3d(v1.nx(),  v1.ny(),  v1.nz());

				  glVertex3d(v1.x()  , v1.y(),   v1.z() );

									

				  glNormal3d(v2.nx(),  v2.ny(),  v2.nz());

				  glVertex3d(v2.x()  , v2.y(),   v2.z() );

								

				  glNormal3d(v3.nx(),  v3.ny(),  v3.nz());

				  glVertex3d(v3.x()  , v3.y(),   v3.z() );

			  glEnd();

	}

  }

}

//----------------------------------------------------------------//

void start_config()

//----------------------------------------------------------------//

{

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

	glOrtho(-1, 1, -1, 1, -1, 1);



	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

	glTranslatef (obj_pos[0], obj_pos[1], obj_pos[2] );

	glMultMatrixf( view_rotate );

	glScalef(scale,scale,scale);



	glEnable (GL_DEPTH_TEST);

	glDepthFunc (GL_LESS);

	glClearDepth(100.0);

	glClear(GL_DEPTH_BUFFER_BIT);

}

//----------------------------------------------------------------//

void display()

//----------------------------------------------------------------//

{

	glClearColor(1.0,1.0,1.0,0.0);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	start_config();

		

  glShadeModel(GL_SMOOTH);

	glEnable(GL_LIGHTING);

	glEnable(GL_NORMALIZE);



	glEnable(GL_LIGHT0);

	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);

	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);

	glLightfv(GL_LIGHT0, GL_POSITION,light0_position);



	glEnable(GL_LIGHT1);

	glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);

	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);

	glLightfv(GL_LIGHT1, GL_POSITION,light1_position);



	glEnable(GL_POLYGON_OFFSET_FILL);

	glPolygonOffset(1,3);



	glEnable(GL_COLOR_MATERIAL);

	glMaterialfv(GL_FRONT, GL_SHININESS, material_shininess );



	

  if(vstar)

  {

    draw_vstar(dim);

  }

 

  if(smooth)

	{

    if(level == 0 ) ch0.draw_smooth();

    if(level == 1 ) ch1.draw_smooth();

    if(level == 2 ) ch2.draw_smooth();

    if(level == 3 ) ch3.draw_smooth();

	}

	if(points)

	{

    if(level == 0 ) ch0.draw_verts();

    if(level == 1 ) ch1.draw_verts();

    if(level == 2 ) ch2.draw_verts();

    if(level == 3 ) ch3.draw_verts();

	}

	if(wire)

	{

   if(level == 0 ) ch0.draw_wire();

   if(level == 1 ) ch1.draw_wire();

   if(level == 2 ) ch2.draw_wire();

   if(level == 3 ) ch3.draw_wire();

	}



  glDisable(GL_LIGHTING);

	glDisable(GL_COLOR_MATERIAL);

	glDisable(GL_POLYGON_OFFSET_FILL);

	glDisable(GL_LIGHT0);

	glDisable(GL_LIGHT1);

	glutSwapBuffers();

}

//---------------------------------------------------------------//

void myGlutIdle( void )

//---------------------------------------------------------------//

{

  if ( glutGetWindow() != window_id)

  glutSetWindow(window_id);

  glutPostRedisplay();

}

//----------------------------------------------------------------//

void Reshape  (int w, int h)

//----------------------------------------------------------------//

{

	view_w = w;

	view_h = h;

	glViewport (0,0,w,h);

}

//----------------------------------------------------------------//

void read_key(int key)

//----------------------------------------------------------------//

{

	switch(key)

	{

	case 1:

    if(level == 0) {

      ch0.read_ply( file );

      ch0.check();

      nverts = ch0.nvert();      

    }

    if(level == 1) { 

      ch1.read_ply( file );

      ch1.check();

      nverts = ch1.nvert();

    }

    if(level == 2) {

      ch2.read_ply( file );

      ch2.check();

      nverts = ch2.nvert();

   }

    if(level == 3) {

      ch3.read_ply( file );

      ch3.check();

      nverts = ch3.nvert();

    }

    simplexid->set_int_limits(0, nverts);

    read_key(5);

  break;



	case 2:

		if(level == 0) ch0.write_ply( file );

		if(level == 1) ch1.write_ply( file );

		if(level == 2) ch2.write_ply( file );

		if(level == 3) ch3.write_ply( file );

	break;



	case 3:

	break;



	case 4:

		view_rot->reset();
		trans_xy->set_x(0.0);
		trans_xy->set_y(0.10);
		zoom->set_z(1.0);
	break;



	case 5:

    star.clear();

    star = test_vstar(vid, dim);

	break;

	

	case 6:

	break;



	case 10:

		exit(1);

	break;

	}

}

//----------------------------------------------------------------//

int main(int argc, char **argv)

//----------------------------------------------------------------//

{

	glutInit(&argc,argv);

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);

	glutInitWindowSize(view_w,view_h);

	glutInitWindowPosition(50,50);

	window_id= glutCreateWindow("Pet_CHE_Test");



	glutDisplayFunc(display);



	GLUI_Master.set_glutReshapeFunc ( Reshape );

    GLUI_Master.set_glutIdleFunc( myGlutIdle );



    glui=  GLUI_Master.create_glui_subwindow( window_id, GLUI_SUBWINDOW_BOTTOM );

	glui-> set_main_gfx_window( window_id );



	//----------------------------------------------------------//

	/*Main Window	*/

	//----------------------------------------------------------//



	//--------PANEL_Main

	GLUI_Panel *panelm = glui->add_panel( "Main", GLUI_PANEL_NONE );

	panelm->set_alignment(GLUI_ALIGN_CENTER);



	//---PANEL1

	GLUI_Panel *panel1 = glui->add_panel_to_panel(panelm, "CHE i/o", GLUI_PANEL_EMBOSSED );

	panel1->set_alignment(GLUI_ALIGN_LEFT);



	//TEXTTEXT

	GLUI_EditText *edittext =

    glui->add_edittext_to_panel(panel1 , "File ", GLUI_EDITTEXT_TEXT, file );

	edittext->set_alignment(GLUI_ALIGN_RIGHT);

	

	//BUTTON

	GLUI_Button *b =

	glui->add_button_to_panel( panel1, "Read Ply" , 1, read_key );

	b->set_alignment(GLUI_ALIGN_RIGHT);



	//BUTTON

	GLUI_Button *b1 =

	glui->add_button_to_panel( panel1, "Write Ply" , 2, read_key );

	b1->set_alignment(GLUI_ALIGN_RIGHT);



	glui->add_column_to_panel(panelm, false);

	

	//---PANEL2

	GLUI_Panel *panel2 = glui->add_panel_to_panel(panelm, "Controls", GLUI_PANEL_EMBOSSED );

	panel2->set_alignment(GLUI_ALIGN_LEFT);



	//ROTATION

	view_rot = glui->add_rotation_to_panel( panel2, "Rotation", view_rotate );

    view_rot->set_spin( 1.0f );



	glui->add_column_to_panel(panel2, false);



	//ZOOM

	zoom = glui->add_translation_to_panel(panel2, "Zoom", GLUI_TRANSLATION_Z, &scale );

	zoom->set_speed( .01f );



	glui->add_column_to_panel(panel2, false);



	//TRANSLATION

	trans_xy = glui->add_translation_to_panel(panel2, "Translation", GLUI_TRANSLATION_XY, obj_pos );

	trans_xy->set_speed( .01f );



	glui->add_column_to_panel(panelm, false);



	//---PANEL3

	GLUI_Panel *panel3 = glui->add_panel_to_panel(panelm, "CHE Level", GLUI_PANEL_EMBOSSED );

	panel3->set_alignment(GLUI_ALIGN_LEFT);



	glui->add_column_to_panel(panel3, false);



	//CHECK RADIO

	GLUI_RadioGroup *group1 =

	glui->add_radiogroup_to_panel(panel3, &level, 3, read_key);

		glui->add_radiobutton_to_group( group1, "L0" );

		glui->add_radiobutton_to_group( group1, "L1" );

		glui->add_radiobutton_to_group( group1, "L2" );

		glui->add_radiobutton_to_group( group1, "L3" );



	glui->add_column_to_panel(panel3, false);



	//CHECK BOX

	glui->add_checkbox_to_panel(panel3, "Smooth", &smooth);

	

	//CHECK BOX

	c=	glui->add_checkbox_to_panel(panel3, "Wire", &wire);

	//c->disable();

		

	//CHECK BOX

	glui->add_checkbox_to_panel(panel3, "Points", &points);



	glui->add_column_to_panel(panelm, false);



 	//---PANEL4

	GLUI_Panel *panel4 = glui->add_panel_to_panel(panelm, "Stars", GLUI_PANEL_EMBOSSED );

	panel4->set_alignment(GLUI_ALIGN_LEFT);



  glui->add_checkbox_to_panel(panel4, "vStars", &vstar);



	GLUI_RadioGroup *group2 =

	glui->add_radiogroup_to_panel(panel4, &dim, 5, read_key);

		glui->add_radiobutton_to_group( group2, "R00" );

		glui->add_radiobutton_to_group( group2, "R02" );



  simplexid = glui->add_spinner_to_panel(panel4, "Vertex id", GLUI_SPINNER_INT, &vid, 5, read_key );

  simplexid->set_int_limits(0, nverts, GLUI_LIMIT_WRAP);

  simplexid->set_speed(0.1);



	glutMainLoop();

}

//----------------------------------------------------------------//
