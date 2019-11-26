#ifndef FILE2_H_INCLUDED
#define FILE2_H_INCLUDED

#include "bitmap_image.hpp"
#include <iostream>
#include <vector>

#define pi (2*acos(0.0))


enum {AMBIENT, DIFFUSE, SPECULAR, REFLECTION};

int recursion_level;

typedef struct Vector3
{
	double x,y,z;

	Vector3 normalize() {
	    double temp = sqrt(x*x + y*y + z*z);
	     x /= temp;
	     y /= temp;
	     z /= temp;
        return  *this;
    }
}vector3;

std::vector<vector3> lights_vect;

vector3 Sub(vector3 a,vector3 b) {
	vector3 ret;
	ret.x=a.x-b.x;
	ret.y=a.y-b.y;
	ret.z=a.z-b.z;
	return ret;
}

vector3 Add(vector3 a,vector3 b) {
	vector3 ret;
	ret.x=a.x+b.x;
	ret.y=a.y+b.y;
	ret.z=a.z+b.z;
	return ret;
}

vector3 cross_product(vector3 a, vector3 b) {
	vector3 ret;
    ret.x = a.y*b.z-a.z*b.y;
    ret.y = a.z*b.x-a.x*b.z;
    ret.z = a.x*b.y-a.y*b.x;


	return ret;
}

double dot_product(Vector3 a,Vector3 b) {
    double ret;
    ret = a.x*b.x + a.y*b.y + a.z*b.z;
    return ret;
}

vector3 scalarProduct(double a,Vector3 b) {
	vector3 ret;
	ret.x=a*b.x;
	ret.y=a*b.y;
	ret.z=a*b.z;
	return ret;
}


class Ray{
    public:
    vector3 start;
    vector3 dir;

    Ray(vector3 start, vector3 dir){
        this->start = start;
        this->dir = dir.normalize();
    }
};


//creating base class
class Object{
    public:
    vector3 reference_point;
    double height, width, length, source_factor = 1.0;
    int Shine;
    double color[3];
    double co_efficients[4];
    std::string name;

    Object(){ }
    virtual void draw() = 0;
    virtual double getIntersectionT(Ray *r) = 0;
    virtual double intersect(Ray *r, double *color_0, double *color_1, double *color_2, int level) = 0;

    void setColor(double R, double G, double B) {color[0] = R; color[1] = G; color[2] = B;};
    void setShine(int s) {Shine = s;};
    void setCoEfficients(double cof_1, double cof_2, double cof_3, double cof_4) { co_efficients[AMBIENT] = cof_1; co_efficients[DIFFUSE]=cof_2; co_efficients[SPECULAR]=cof_3; co_efficients[REFLECTION]=cof_4;};
};

std::vector<Object*> objects;
//derived class
class Sphere: public Object{
    public:
    Sphere(vector3 Center, double Radius) {
        reference_point.x=Center.x;
        reference_point.y=Center.y;
        reference_point.z=Center.z;
        length=Radius;
        name = "sphere";
    }

    void draw()
    {
        double radius = length;
        int slices = 20; int stacks = 20;
        vector3 points[100][100];
        int i,j;
        double h,r;
        glColor3f(color[0],color[1],color[2]);
        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0;i<stacks;i++)
        {
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);{
                    //upper hemisphere
                    glVertex3f(points[i][j].x + reference_point.x,points[i][j].y + reference_point.y,points[i][j].z + reference_point.z);
                    glVertex3f(points[i][j+1].x + reference_point.x,points[i][j+1].y + reference_point.y,points[i][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j+1].x + reference_point.x,points[i+1][j+1].y + reference_point.y,points[i+1][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j].x + reference_point.x,points[i+1][j].y + reference_point.y,points[i+1][j].z + reference_point.z);
                    //lower hemisphere
                    glVertex3f(points[i][j].x + reference_point.x,points[i][j].y + reference_point.y,-points[i][j].z + reference_point.z);
                    glVertex3f(points[i][j+1].x + reference_point.x,points[i][j+1].y + reference_point.y,-points[i][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j+1].x + reference_point.x,points[i+1][j+1].y + reference_point.y,-points[i+1][j+1].z + reference_point.z);
                    glVertex3f(points[i+1][j].x + reference_point.x,points[i+1][j].y + reference_point.y,-points[i+1][j].z + reference_point.z);
                }glEnd();
            }
        }
    }

    double getIntersectionT(Ray* r) {
        vector3 temp_start = Sub(r->start, reference_point);

        double A = dot_product(r->dir, r->dir);
        double B = 2 * dot_product(r->dir, temp_start);
        double C = dot_product(temp_start, temp_start) - length * length;
        double D = B*B - 4*A*C;


        if(D<0)
            return -1;
        double t1 = (-(B - sqrt(D)))/(2*A);
        double t2 = (-(B + sqrt(D)))/(2*A);

        if(t1<t2)
            return t1;
        else
            return t2;
    }


    vector3 getNormal(vector3 intersecPoint){
        vector3 normal = Sub(intersecPoint, reference_point);
        return normal.normalize();
    }

    vector3 getReflection(Ray *r, vector3 normal){
        vector3 temp = Sub(r->dir,scalarProduct(2*dot_product(r->dir,normal),normal));
        return temp.normalize();
    }


    //intersect function
    double intersect(Ray *r, double *color_0, double *color_1, double *color_2, int level) {
        double t = getIntersectionT(r);
        //std::cout<<t<<std::endl;
        if (t <= 0) {
            return -1;
        }

        if (level == 0 ) {
            return t;
        }
        //std::cout<<t<<std::endl;
        *color_0 = color[0] * co_efficients[AMBIENT];
        *color_1 = color[1] * co_efficients[AMBIENT];
        *color_2 = color[2] * co_efficients[AMBIENT];
        if(*color_0 < 0)
            *color_0 = 0;
        else if(*color_0 > 1)
            *color_0 = 1;
        if(*color_0 < 0)
            *color_1 = 0;
        else if(*color_0 > 1)
            *color_1 = 1;
            if(*color_0 < 0)
            *color_2 = 0;
        else if(*color_0 > 1)
            *color_2 = 1;

        vector3 intersectionPoint = Add(r->start, scalarProduct(t,r->dir));
        vector3 normal=getNormal(intersectionPoint);
        vector3 reflection=getReflection(r, normal);

        for(int i = 0; i<lights_vect.size(); ++i){
            vector3 dir = Sub(lights_vect[i], intersectionPoint);
            double distance = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
            dir = dir.normalize();

            vector3 start = Add(intersectionPoint,scalarProduct(1,dir));
            Ray L(start, dir);

            bool vivid = false;

            for (int j=0; j < objects.size(); j++) {

                double temp_t = objects[j]->getIntersectionT(&L);

                if(temp_t > 0 || abs(temp_t) > distance) {
                    continue;
                }

                vivid = true;
                break;
            }

            if (vivid){

                double lambert = dot_product(L.dir, normal);
                double phong = pow(dot_product(reflection, r->dir), Shine);

                lambert = lambert > 0 ? lambert : 0;
                phong = phong > 0 ? phong : 0;


                *color_0 += source_factor * lambert * co_efficients[DIFFUSE] * color[0];
                *color_0 += source_factor * phong * co_efficients[SPECULAR] * color[0];
                *color_1 += source_factor * lambert * co_efficients[DIFFUSE] * color[1];
                *color_1 += source_factor * phong * co_efficients[SPECULAR] * color[1];
                *color_2 += source_factor * lambert * co_efficients[DIFFUSE] * color[2];
                *color_2 += source_factor * phong * co_efficients[SPECULAR] * color[2];
                if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;

            }

            if (level < recursion_level) {
                vector3 start = Add(intersectionPoint, scalarProduct(1,reflection));
                Ray reflectionRay(start, reflection);

                int nearest=-1;
                double min_relection = 99999;
                double color_reflection_0, color_reflection_1, color_reflection_2;

                for (int k=0; k < objects.size(); k++) {

                    double temp_reflection = objects[k]->getIntersectionT(&reflectionRay);

                    if(temp_reflection <= 0)
                    {
                        continue;
                    }
                    else if(temp_reflection < min_relection)
                    {
                        min_relection = temp_reflection;
                        nearest = k;
                    }
                }

                if(nearest!=-1) {

                    double debug_grabage = objects[nearest]->intersect(&reflectionRay, &color_reflection_0, &color_reflection_1, &color_reflection_2, level+1);
                    if(debug_grabage == -1)
                        return -1;
                    *color_0 += color_reflection_0 * co_efficients[REFLECTION];
                    *color_1 += color_reflection_1 * co_efficients[REFLECTION];
                    *color_2 += color_reflection_2 * co_efficients[REFLECTION];
                    if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;

                }

            }



        }



        return t;
    }

};


class Floor: public Object {

public:

    //bitmap_image bitImage;
    //double biHeight, biWidth;

    Floor(double floorWidth, double tileWidth) {
        reference_point.x = -floorWidth/2;
        reference_point.y = -floorWidth/2;
        reference_point.z =  0;
        length = tileWidth;
        name = "Floor";
        //bitImage = bitmap_image("capture.bmp");
        //biHeight = bitImage.height()/1000.0;
        //biWidth = bitImage.width()/1000.0;
    }

    void draw() {

        int numOfTiles = abs(reference_point.x*2/length);

        for (int i=0; i<numOfTiles; i++) {
            for (int j=0; j<numOfTiles; j++) {

                if ((i+j)%2) {
                    glColor3f(0, 0, 0);
                } else {
                    glColor3f(1, 1, 1);
                }

                glBegin(GL_QUADS);
                {
                    glVertex3f(reference_point.x+length*i, reference_point.y+length*j, reference_point.z);
                    glVertex3f(reference_point.x+length*(i+1), reference_point.y+length*j, reference_point.z);
                    glVertex3f(reference_point.x+length*(i+1), reference_point.y+length*(j+1), reference_point.z);
                    glVertex3f(reference_point.x+length*i, reference_point.y+length*(j+1), reference_point.z);
                }
                glEnd();
            }
        }

    }

    vector3 getNormal(){
        vector3 v = {0, 0, 1};
        return v;
    }

    double getIntersectionT(Ray* r) {

        vector3 normal = getNormal();

        double temp = dot_product(normal, r->start) * (-1) / dot_product(normal, r->dir);
        //std::cout<<"bal"<<std::endl;
        return temp;
    }

    vector3 getReflection(Ray *r, vector3 normal){
        vector3 temp = Sub(r->dir,scalarProduct(2*dot_product(r->dir,normal),normal));
        return temp.normalize();
    }

    double intersect(Ray* r, double *color_0, double *color_1, double *color_2, int level) {

        double temp = getIntersectionT(r);
        //std::cout<<temp<<std::endl;

        vector3 intersectionPoint = Add(r->start, scalarProduct(temp, r->dir));

        double x_Min = reference_point.x;
        double x_Max = x_Min * (-1);

        double y_Min = reference_point.y;
        double y_Max = y_Min * (-1);



        if (x_Min > intersectionPoint.x || intersectionPoint.x > x_Max ||
                y_Min > intersectionPoint.y || intersectionPoint.y > y_Max) {
            return -1;
        }

        int temp_x = intersectionPoint.x / length;
        int temp_y = intersectionPoint.y / length;

        if ((temp_x+temp_y)%2) {
            color[0] = 0;
            color[1] = 0;
            color[2] = 0;
        } else {
            color[0] = 1;
            color[1] = 1;
            color[2] = 1;

        }

        //unsigned char red, green, blue;
        //int x = (intersectionPoint.x + abs(reference_point.x)) * biWidth;
        //int y = (intersectionPoint.y + abs(reference_point.y)) * biHeight;

        //bitImage.get_pixel(x, y, red, green, blue);
        *color_0 = color[0] * co_efficients[AMBIENT]; //* red / 255.0;
        *color_1 = color[1] * co_efficients[AMBIENT]; //* green / 255.0;
        *color_2 = color[2] * co_efficients[AMBIENT]; //* blue / 255.0;
        if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;
        //std::cout << *color_0 << " "<< *color_1 << " "<< *color_2 <<std::endl;
        //return temp;
        vector3 normal=getNormal();
        vector3 reflection=getReflection(r, normal);

        for(int i = 0; i<lights_vect.size(); ++i){
            vector3 dir = Sub(lights_vect[i], intersectionPoint);
            double distance = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
            dir = dir.normalize();

            vector3 start = Add(intersectionPoint,scalarProduct(1,dir));
            Ray L(start, dir);

            bool vivid = false;

            for (int j=0; j < objects.size(); j++) {

                double temp_t = objects[j]->getIntersectionT(&L);

                if(temp_t > 0 || abs(temp_t) > distance) {
                    continue;
                }

                vivid = true;
                break;
            }
            if (vivid){

                double lambert = dot_product(L.dir, normal);
                double phong = pow(dot_product(reflection, r->dir), Shine);

                lambert = lambert > 0 ? lambert : 0;
                phong = phong > 0 ? phong : 0;


                *color_0 += source_factor * lambert * co_efficients[DIFFUSE] * color[0];
                *color_0 += source_factor * phong * co_efficients[SPECULAR] * color[0];
                *color_1 += source_factor * lambert * co_efficients[DIFFUSE] * color[1];
                *color_1 += source_factor * phong * co_efficients[SPECULAR] * color[1];
                *color_2 += source_factor * lambert * co_efficients[DIFFUSE] * color[2];
                *color_2 += source_factor * phong * co_efficients[SPECULAR] * color[2];
                if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;

            }

            if (level < recursion_level) {
                vector3 start = Add(intersectionPoint, scalarProduct(1,reflection));
                Ray reflectionRay(start, reflection);

                int nearest=-1;
                double min_relection = 99999;
                double color_reflection_0, color_reflection_1, color_reflection_2;

                for (int k=0; k < objects.size(); k++) {

                    double temp_reflection = objects[k]->getIntersectionT(&reflectionRay);

                    if(temp_reflection <= 0)
                    {
                        continue;
                    }
                    else if(temp_reflection < min_relection)
                    {
                        min_relection = temp_reflection;
                        nearest = k;
                    }
                }

                if(nearest!=-1) {
                    double debug_grabage = objects[nearest]->intersect(&reflectionRay, &color_reflection_0, &color_reflection_1, &color_reflection_2, level+1);
                    if(debug_grabage==-1) return -1;
                    *color_0 += color_reflection_0 * co_efficients[REFLECTION];
                    *color_1 += color_reflection_1 * co_efficients[REFLECTION];
                    *color_2 += color_reflection_2 * co_efficients[REFLECTION];
                    if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;

                }

            }



        }
        //std::cout<<temp<<std::endl;
        return temp;
    }



};

class Triangle: public Object {
    public:
    vector3 a, b, c;

    Triangle(vector3 a, vector3 b, vector3 c) {
        this->a.x = a.x;
        this->a.y = a.y;
        this->a.z = a.z;
        this->b.x = b.x;
        this->b.y = b.y;
        this->b.z = b.z;
        this->c.x = c.x;
        this->c.y = c.y;
        this->c.z = c.z;
        name = "Triangle";
    }

    void draw() {
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    vector3 getNormal(vector3 intersection) {

        vector3 temp_1 = Sub(b, a);
        vector3 temp_2 = Sub(c, a);

        vector3 normal = cross_product(temp_1, temp_2);

        return normal.normalize();
    }

    double getIntersectionT(Ray* r) {
        const float EPSILON = 0.0000001;
        vector3 e1 = Sub(b, a);
        vector3 e2 = Sub(c, a);

        vector3 c_product = cross_product(r->dir, e2);
        double determinant = dot_product(e1, c_product);

        if(determinant > -EPSILON && determinant < EPSILON) {
            return -1;
        }

        double invert_det = 1.0 / determinant;

        vector3 t_ray = Sub(r->start, a);

        double u = dot_product(t_ray, c_product) * invert_det;

        if(u < 0 || u > 1) {
            return -1;
        }

        vector3 temp_product = cross_product(t_ray, e1);

        double v = dot_product(r->dir, temp_product) * invert_det;

        if(v < 0 || u + v  > 1) {
            return -1;
        }

        double ret = dot_product(e2, temp_product) * invert_det;

        if(ret > EPSILON) { //ray intersection
            return ret;
        }

        return -1;
    }


    vector3 getReflection(Ray *r, vector3 normal){
        vector3 temp = Sub(r->dir,scalarProduct(2*dot_product(r->dir,normal),normal));
        return temp.normalize();
    }


double intersect(Ray *r, double *color_0, double *color_1, double *color_2, int level) {
        double t = getIntersectionT(r);
        if (t <= 0) {
            return -1;
        }

        if (level == 0 ) {
            return t;
        }

        *color_0 = color[0] * co_efficients[AMBIENT];*color_1 = color[1] * co_efficients[AMBIENT];*color_2 = color[2] * co_efficients[AMBIENT];
        if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;

        vector3 intersectionPoint = Add(r->start, scalarProduct(t,r->dir));
        vector3 normal=getNormal(intersectionPoint);
        vector3 reflection=getReflection(r, normal);

        for(int i = 0; i<lights_vect.size(); ++i){
            vector3 dir = Sub(lights_vect[i], intersectionPoint);
            double distance = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
            dir = dir.normalize();

            if (dot_product(dir, normal) > 0) {
                normal = scalarProduct(-1,normal);
            }

            vector3 start = Add(intersectionPoint,scalarProduct(1,dir));
            Ray L(start, dir);

            bool vivid = false;

            for (int j=0; j < objects.size(); j++) {

                double temp_t = objects[j]->getIntersectionT(&L);

                if(temp_t > 0 || abs(temp_t) > distance) {
                    continue;
                }

                vivid = true;
                break;
            }
            if (vivid){

                double lambert = dot_product(L.dir, normal);
                double phong = pow(dot_product(reflection, r->dir), Shine);

                lambert = lambert > 0 ? lambert : 0;
                phong = phong > 0 ? phong : 0;


                *color_0 += source_factor * lambert * co_efficients[DIFFUSE] * color[0];
                *color_0 += source_factor * phong * co_efficients[SPECULAR] * color[0];
                *color_1 += source_factor * lambert * co_efficients[DIFFUSE] * color[1];
                *color_1 += source_factor * phong * co_efficients[SPECULAR] * color[1];
                *color_2 += source_factor * lambert * co_efficients[DIFFUSE] * color[2];
                *color_2 += source_factor * phong * co_efficients[SPECULAR] * color[2];
                if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;

            }

            if (level < recursion_level) {
                vector3 start = Add(intersectionPoint, scalarProduct(1,reflection));
                Ray reflectionRay(start, reflection);

                int nearest=-1;
                double min_relection = 99999;
                double color_reflection_0, color_reflection_1, color_reflection_2;

                for (int k=0; k < objects.size(); k++) {

                    double temp_reflection = objects[k]->getIntersectionT(&reflectionRay);

                    if(temp_reflection <= 0)
                    {
                        continue;
                    }
                    else if(temp_reflection < min_relection)
                    {
                        min_relection = temp_reflection;
                        nearest = k;
                    }
                }

                if(nearest!=-1) {

                    double debug_grabage = objects[nearest]->intersect(&reflectionRay, &color_reflection_0, &color_reflection_1, &color_reflection_2, level+1);

                    if(debug_grabage==-1) return -1;*color_0 += color_reflection_0 * co_efficients[REFLECTION];
                    *color_1 += color_reflection_1 * co_efficients[REFLECTION];
                    *color_2 += color_reflection_2 * co_efficients[REFLECTION];
                    if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;

                }

            }



        }



        return t;
    }




};

class generalQuad: public Object{
    public:
    double A, B, C, D, E, F, G, H, I, J;

    generalQuad(double co_eff[10], vector3 reff,double width, double height,  double length) {
        A = co_eff[0];
        B = co_eff[1];
        C = co_eff[2];
        D = co_eff[3];
        E = co_eff[4];
        F = co_eff[5];
        G = co_eff[6];
        H = co_eff[7];
        I = co_eff[8];
        J = co_eff[9];
        reference_point.x = reff.x;
        reference_point.y = reff.y;
        reference_point.z = reff.z;
        this->width = width;
        this->height = height;
        this->length = length;
        name = "general";
    }

    void draw() {}

    vector3 getNormal(vector3 intersection) {

        double u = 2 * A * intersection.x + D * intersection.y + F * intersection.z  + G;
        double v = 2 * B * intersection.y + D * intersection.x + E * intersection.z  + H;
        double z = 2 * C * intersection.z + E * intersection.y + F * intersection.x  + I;

        vector3 normal = {u, v, z};

        return normal.normalize();
    }

    double getIntersectionT(Ray* r) {

        double a = A * r->dir.x * r->dir.x + B * r->dir.y * r->dir.y + C * r->dir.z * r->dir.z + D * r->dir.x * r->dir.y + E * r->dir.y * r->dir.z + F * r->dir.z * r->dir.x;
        double b = 2 * (A * r->start.x * r->dir.x + B * r->start.y * r->dir.y + C * r->start.z * r->dir.z) + D * (r->start.x * r->dir.y + r->dir.x * r->start.y) + E * (r->start.y * r->dir.z + r->dir.y * r->start.z) + F * (r->start.z * r->dir.x + r->dir.z * r->start.x) + G * r->dir.x + H * r->dir.y + I * r->dir.z;
        double c = A * r->start.x * r->start.x + B * r->start.y * r->start.y + C * r->start.z * r->start.z + D * r->start.x * r->start.y + E * r->start.y * r->start.z + F * r->start.z * r->start.x + G * r->start.x + H * r->start.y + I * r->start.z + J;
        double d = b*b - 4*a*c;

        if (d < 0) {
            return -1;
        }

        double t1 = (- b + sqrt(d)) / (2.0*a);
        double t2 = (- b - sqrt(d)) / (2.0*a);


        vector3 intersectionPoint_1 = Add(r->start, scalarProduct(t1, r->dir));
        vector3 intersectionPoint_2 = Add(r->start, scalarProduct(t2,r->dir));

        double x_min = reference_point.x;
        double x_max = x_min + length;

        double y_min = reference_point.y;
        double y_max = y_min + width;

        double z_min = reference_point.z;
        double z_max = z_min + height;


        bool flag1 = (length > 0 && (x_min > intersectionPoint_1.x || intersectionPoint_1.x > x_max) ||
                      width > 0 && (y_min > intersectionPoint_1.y || intersectionPoint_1.y > y_max) ||
                      height > 0 && (z_min > intersectionPoint_1.z || intersectionPoint_1.z > z_max));

        bool flag2 = (length > 0 && (x_min > intersectionPoint_2.x || intersectionPoint_2.x > x_max) ||
                      width > 0 && (y_min > intersectionPoint_2.y || intersectionPoint_2.y > y_max) ||
                      height > 0 && (z_min > intersectionPoint_2.z || intersectionPoint_2.z > z_max));


        if (flag1 && flag2) {
            return -1;
        }
        else if (flag1)
        {
            return t2;
        }
        else if (flag2){
            return t1;
        }
        else {
            if (t1 < t2) {
                return t1;
            }
            else{
                return t2;
            }

        }
    }


    vector3 getReflection(Ray *r, vector3 normal){
        vector3 temp = Sub(r->dir,scalarProduct(2*dot_product(r->dir,normal),normal));
        return temp.normalize();
    }


    double intersect(Ray *r, double *color_0, double *color_1, double *color_2, int level) {
        double t = getIntersectionT(r);
        if (t <= 0) {
            return -1;
        }

        if (level == 0 ) {
            return t;
        }

        *color_0 = color[0] * co_efficients[AMBIENT];*color_1 = color[1] * co_efficients[AMBIENT];*color_2 = color[2] * co_efficients[AMBIENT];
        if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;

        vector3 intersectionPoint = Add(r->start, scalarProduct(t,r->dir));
        vector3 normal=getNormal(intersectionPoint);
        vector3 reflection=getReflection(r, normal);

        for(int i = 0; i<lights_vect.size(); ++i){
            vector3 dir = Sub(lights_vect[i], intersectionPoint);
            double distance = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
            dir = dir.normalize();

            if (dot_product(dir, normal) > 0) {
                normal = scalarProduct(-1,normal);
            }

            vector3 start = Add(intersectionPoint,scalarProduct(1,dir));
            Ray L(start, dir);

            bool vivid = false;

            for (int j=0; j < objects.size(); j++) {

                double temp_t = objects[j]->getIntersectionT(&L);

                if(temp_t > 0 || abs(temp_t) > distance) {
                    continue;
                }

                vivid = true;
                break;
            }
            if (vivid){

                double lambert = dot_product(L.dir, normal);
                double phong = pow(dot_product(reflection, r->dir), Shine);

                lambert = lambert > 0 ? lambert : 0;
                phong = phong > 0 ? phong : 0;


                *color_0 += source_factor * lambert * co_efficients[DIFFUSE] * color[0];
                *color_0 += source_factor * phong * co_efficients[SPECULAR] * color[0];
                *color_1 += source_factor * lambert * co_efficients[DIFFUSE] * color[1];
                *color_1 += source_factor * phong * co_efficients[SPECULAR] * color[1];
                *color_2 += source_factor * lambert * co_efficients[DIFFUSE] * color[2];
                *color_2 += source_factor * phong * co_efficients[SPECULAR] * color[2];
                if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;
            }

            if (level < recursion_level) {
                vector3 start = Add(intersectionPoint, scalarProduct(1,reflection));
                Ray reflectionRay(start, reflection);

                int nearest=-1;
                double min_relection = 99999;
                double color_reflection_0, color_reflection_1, color_reflection_2;

                for (int k=0; k < objects.size(); k++) {

                    double temp_reflection = objects[k]->getIntersectionT(&reflectionRay);

                    if(temp_reflection <= 0)
                    {
                        continue;
                    }
                    else if(temp_reflection < min_relection)
                    {
                        min_relection = temp_reflection;
                        nearest = k;
                    }
                }

                if(nearest!=-1) {

                    double debug_grabage = objects[nearest]->intersect(&reflectionRay, &color_reflection_0, &color_reflection_1, &color_reflection_2, level+1);
                    if(debug_grabage==-1) return -1;
                    *color_0 += color_reflection_0 * co_efficients[REFLECTION];
                    *color_1 += color_reflection_1 * co_efficients[REFLECTION];
                    *color_2 += color_reflection_2 * co_efficients[REFLECTION];
                    if(*color_0 < 0)
                    *color_0 = 0;
                else if(*color_0 > 1)
                    *color_0 = 1;
                if(*color_0 < 0)
                    *color_1 = 0;
                else if(*color_0 > 1)
                    *color_1 = 1;
                if(*color_0 < 0)
                    *color_2 = 0;
                else if(*color_0 > 1)
                    *color_2 = 1;
                }

            }



        }



        return t;
    }





};
#endif
