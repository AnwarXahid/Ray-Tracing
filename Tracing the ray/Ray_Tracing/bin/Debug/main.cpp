#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <windows.h>
#include <glut.h>
#include "FILE2.h"
#include <math.h>
#include <iostream>
#include <typeinfo>


#define pi (2*acos(0.0))
#define Window_width 500
#define Window_height 500
#define fovy 80

using namespace std;


double cameraHeight;
double cameraAngle;
int drawaxes;


typedef struct point
{
	double x,y,z;
}Point;

Point pos, u, r, l;

vector<point> lights;
double lightThickness = 1.0;
double image_width;

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 0, 0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);
		}glEnd();
		glColor3f(0, 1.0, 0);
		glBegin(GL_LINES);{
			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);
        }glEnd();
        glColor3f(0, 0, 1.0);
        glBegin(GL_LINES);{
			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}

point Sub(point a,point b) {
	point ret;
	ret.x=a.x-b.x;
	ret.y=a.y-b.y;
	ret.z=a.z-b.z;
	return ret;
}

point Add(point a,point b) {
	point ret;
	ret.x=a.x+b.x;
	ret.y=a.y+b.y;
	ret.z=a.z+b.z;
	return ret;
}


point scalarProduct(double a,point b) {
	point ret;
	ret.x=a*b.x;
	ret.y=a*b.y;
	ret.z=a*b.z;
	return ret;
}


void drawPoint(point pt) {

    glColor3f(1.0, 1.0, 1.0);

    glBegin(GL_QUADS);
    {
        glVertex3f(pt.x+lightThickness, pt.y, pt.z+lightThickness);
        glVertex3f(pt.x+lightThickness, pt.y, pt.z-lightThickness);
        glVertex3f(pt.x-lightThickness, pt.y, pt.z-lightThickness);
        glVertex3f(pt.x-lightThickness, pt.y, pt.z+lightThickness);
    }
    glEnd();



}

Point dir_of_perp_vector(Point vect, Point perp, int dir)
{
    double c = cos(pi/180);
    double s = dir * sin(pi/180);
    Point point;
    point.x = c * vect.x + s * perp.x;
    point.y = c * vect.y + s * perp.y;
    point.z = c * vect.z + s * perp.z;
    c = sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
    point.x /= c;
    point.y /= c;
    point.z /= c;
    return point;
}

void capture(){
    point** frameBuffer;
    point colour;
    colour.x = 0; colour.y = 0; colour.z = 0;
    frameBuffer = new point* [int(image_width)];
    for(int i = 0; i < image_width; ++i){
        frameBuffer[i] = new point[int(image_width)];
        for(int j = 0; j < image_width; ++j){
            frameBuffer[i][j] = colour;
        }
    }

    double plane_distance = (Window_height/2) / tan(fovy*pi/360);
    point topleft = Add(Sub(Add(pos,scalarProduct(plane_distance, l)),scalarProduct((Window_width*1.0)/2, r)),scalarProduct((Window_height*1.0)/2, u));

    double du=(Window_width*1.0)/image_width;
    double dv=(Window_height*1.0)/image_width;
    //cout<<"debug"<<image_width;
    for(int i =0; i<image_width;++i){
        //cout<<"debug";
        for(int j =0; j < image_width; ++j){
            point corner = Sub(Add(topleft,scalarProduct(j*du,r)),scalarProduct(i*dv,u));
            vector3 point_to_vect, vect_pos;
            Point temp_p = Sub(corner, pos);
            point_to_vect.x = temp_p.x; point_to_vect.y = temp_p.y; point_to_vect.z = temp_p.z;
            vect_pos.x = pos.x; vect_pos.y = pos.y; vect_pos.z = pos.z;
            Ray ray(vect_pos, point_to_vect);
            int nearest=-1;
            double color_0,color_1,color_2;
            double t_min = 99999;
            //cout<<"debug";
            for (int k=0; k < objects.size(); k++) {
                //cout<<k;
                double t = objects[k]->intersect(&ray, &color_0, &color_1, &color_2, 0);
                if(t<=0)
                    continue;
                else if(t<t_min){
                    t_min = t;
                    nearest = k;

                }
            }

            if(nearest!=-1){
                double t = objects[nearest]->intersect(&ray, &color_0, &color_1, &color_2, 1);

                frameBuffer[i][j].x = color_0;
                frameBuffer[i][j].y = color_1;
                frameBuffer[i][j].z = color_2;
                //cout<<objects[nearest]->name<<endl;
                //cout << color_0 << " "<< color_1 << " "<< color_2 <<endl;

            }
        }
    }


    bitmap_image image(image_width, image_width);

    for (int i=0; i<image_width; i++) {
        for (int j=0; j<image_width; j++) {
            image.set_pixel(j, i, frameBuffer[i][j].x*255, frameBuffer[i][j].y*255, frameBuffer[i][j].z*255);
            //cout<<frameBuffer[i][j].x*255;
        }
    }

    image.save_image("capture.bmp");
    cout<<"image saved successfully"<<endl;
}


void keyboardListener(unsigned char key, int x,int y){
	switch(key){

	case '0': {
            capture();
            break;
		}
		case '1': {
            Point l1 = dir_of_perp_vector(l, r, -1);
            r = dir_of_perp_vector(r, l, 1);
            l = l1;
            break;
		}
		case '2': {
            Point l1 = dir_of_perp_vector(l, r, 1);
            r = dir_of_perp_vector(r, l, -1);
            l = l1;
            break;
		}

		case '3': {
            Point u1 = dir_of_perp_vector(u, l, -1);
            l = dir_of_perp_vector(l, u, 1);
            u = u1;
            break;
        }
        case '4': {
            Point u1 = dir_of_perp_vector(u, l, 1);
            l = dir_of_perp_vector(l, u, -1);
            u = u1;
            break;
        }
        case '5': {
            Point r1 = dir_of_perp_vector(r, u, -1);
            u = dir_of_perp_vector(u, r, 1);
            r = r1;
            break;
        }
        case '6':{
            Point r1 = dir_of_perp_vector(r, u, 1);
            u = dir_of_perp_vector(u, r, -1);
            r = r1;
            break;
        }
		default:
			break;
	}
}

void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN: {
            pos.x -= l.x*3;
            pos.y -= l.y*3;
            pos.z -= l.z*3;
            break;
		}
		case GLUT_KEY_UP: {
            pos.x += l.x*3;
            pos.y += l.y*3;
            pos.z += l.z*3;
            break;
		}
		case GLUT_KEY_RIGHT: {
            pos.x += r.x*3;
            pos.y += r.y*3;
            pos.z += r.z*3;
            break;
		}
		case GLUT_KEY_LEFT: {
            pos.x -= r.x*3;
            pos.y -= r.y*3;
            pos.z -= r.z*3;
            break;
		}
		case GLUT_KEY_PAGE_UP: {
            pos.x += u.x*3;
            pos.y += u.y*3;
            pos.z += u.z*3;
            break;
		}
		case GLUT_KEY_PAGE_DOWN: {
            pos.x -= u.x*3;
            pos.y -= u.y*3;
            pos.z -= u.z*3;
            break;
		}
		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;
		default:
			break;
	}
}


void display(){

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

    gluLookAt(pos.x, pos.y, pos.z, pos.x+l.x, pos.y+l.y, pos.z+l.z, u.x, u.y, u.z);

	glMatrixMode(GL_MODELVIEW);

	drawAxes();
	//           infinite drawing loop            //
	for (int i=0; i < objects.size(); i++) {
        objects[i]->draw();
    }

    for (int i=0; i < lights.size(); i++) {
        drawPoint(lights[i]);
    }

	glutSwapBuffers();
}


void animate(){
	glutPostRedisplay();
}

void init(){
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;

	glClearColor(0,0,0,0);

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();
    gluPerspective(fovy,	1,	1,	1000.0);
}


void loadTestData(){
    image_width = 768;
    recursion_level = 3;
    Object *temp;
    vector3 Center = {0, 0, 10};
    double Radius = 10;
    temp=new Sphere(Center, Radius); // Center(0,0,10), Radius 10
    temp->setColor(1,0,0);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);

    point light1 = {-50,50,50};
    vector3 light1_vect = {-50,50,50};
    lights.push_back(light1);
    lights_vect.push_back(light1_vect);

    temp=new Floor(1000, 20);
    temp->setCoEfficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);
    objects.push_back(temp);


}


void loadActualData(){
    freopen("scene.txt", "r", stdin);

    cin>>recursion_level;
    cin>>image_width;

    int numOfObjects;
    cin>>numOfObjects;

    string str;
    Object *temp_object;

    for(int i = 0; i<numOfObjects; ++i){
        cin>>str;
        double x, y, z, radius;
        if(str=="sphere"){
            cin>>x>>y>>z;
            vector3 center = {x, y, z};

            cin>>radius;
            temp_object = new Sphere(center, radius);

            cin>>x>>y>>z;
            temp_object->setColor(x, y, z);

            cin>>x>>y>>z>>radius;
            temp_object->setCoEfficients(x, y, z, radius);

            cin>>x;
            temp_object->setShine(x);

            objects.push_back(temp_object);
        }
        else if(str=="triangle"){
            cin>>x>>y>>z;
            vector3 v1 = {x, y, z};
            cin>>x>>y>>z;
            vector3 v2 = {x, y, z};
            cin>>x>>y>>z;
            vector3 v3 = {x, y, z};

            temp_object = new Triangle(v1,v2,v3);

            cin>>x>>y>>z;
            temp_object->setColor(x, y, z);

            cin>>x>>y>>z>>radius;
            temp_object->setCoEfficients(x, y, z, radius);

            cin>>x;
            temp_object->setShine(x);

            objects.push_back(temp_object);
        }
        else if(str=="general"){
            double coefficient[10];
            for (int l=0; l<10; l++) {
                cin>>coefficient[l];
            }

            cin>>x>>y>>z;
            vector3 refference = {x, y, z};

            cin>>x>>y>>z;
            temp_object = new generalQuad(coefficient, refference, y, z, x);

            cin>>x>>y>>z;
            temp_object->setColor(x, y, z);

            cin>>x>>y>>z>>radius;
            temp_object->setCoEfficients(x, y, z, radius);

            cin>>x;
            temp_object->setShine(x);

            objects.push_back(temp_object);
        }
    }

    cin>>numOfObjects;
    for (int i=0; i<numOfObjects; i++) {
        double x, y, z;
        cin>>x>>y>>z;

        point light1 = {x, y, z};
        vector3 light1_vect = {x, y, z};
        lights.push_back(light1);
        lights_vect.push_back(light1_vect);
    }


    temp_object = new Floor(1000, 20);
    temp_object->setCoEfficients(0.4,0.2,0.2,0.2);
    temp_object->setShine(1);
    objects.push_back(temp_object);

    return ;
}


int main(int argc, char **argv){
    //loadTestData();      //for customizing different shape and configuration
    loadActualData();
    //cout<<"debug";




    double d = 1 / sqrt(2);
    pos = {100, 100, 50};
    u = {0, 0, 1};
    r = {-d, d, 0};
    l = {-d, -d, 0};
	glutInit(&argc,argv);
	glutInitWindowSize(Window_width, Window_height);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("Ray Tracing");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
