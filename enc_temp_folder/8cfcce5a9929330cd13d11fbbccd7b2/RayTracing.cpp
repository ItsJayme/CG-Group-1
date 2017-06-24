#include "stdafx.h"
#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/glut.h>
#include "raytracing.h"
#include <iostream>


//temporary variables
//these are only used to illustrate 
//a simple debug drawing. A ray 
Vec3Df testRayOrigin;
Vec3Df testRayDestination;
std::vector<Vec3Df> Kd;//diffuse coefficient per vertex
std::vector<Vec3Df> Ks;//specularity coefficient per vertex
std::vector<float> Shininess;//exponent for phong and blinn-phong specularities


//use this function for any preprocessing of the mesh.
void init()
{
	//load the mesh file
	//please realize that not all OBJ files will successfully load.
	//Nonetheless, if they come from Blender, they should, if they 
	//are exported as WavefrontOBJ.
	//PLEASE ADAPT THE LINE BELOW TO THE FULL PATH OF THE dodgeColorTest.obj
	//model, e.g., "C:/temp/myData/GraphicsIsFun/dodgeColorTest.obj", 
	//otherwise the application will not load properly
   
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);

	MyMesh.loadMesh("cube.obj", true);
	Kd.resize(MyMesh.vertices.size(), Vec3Df(0.5, 0.5, 0.5));
	Ks.resize(MyMesh.vertices.size(), Vec3Df(0.5, 0.5, 0.5));
	Shininess.resize(MyMesh.vertices.size(), 3);
}

void calculateMainBox(float &Xmax, float &Xmin, float &Ymax, float &Ymin, float &Zmax, float &Zmin)
{
	Xmax = -HUGE_VALF;
	Xmin = HUGE_VALF;
	Ymax = -HUGE_VALF;
	Ymin = HUGE_VALF;
	Zmax = -HUGE_VALF;
	Zmin = HUGE_VALF;
	for (unsigned int i = 0; i < MyMesh.triangles.size(); i++) {
		Triangle currenttriangle = MyMesh.triangles[i];
		Vec3Df v0 = MyMesh.vertices[currenttriangle.v[0]].p;
		Vec3Df v1 = MyMesh.vertices[currenttriangle.v[1]].p;
		Vec3Df v2 = MyMesh.vertices[currenttriangle.v[2]].p;
			if (v0[0] > Xmax) {
				Xmax = v0[0];
			}
			if (v0[0] < Xmin) {
				Xmin = v0[0];
			}
			if (v1[0] > Xmax) {
				Xmax = v1[0];
			}
			if (v1[0] < Xmin) {
				Xmin = v1[0];
			}
			if (v2[0] > Xmax) {
				Xmax = v2[0];
			}
			if (v2[0] < Xmin) {
				Xmin = v2[0];
			}
			if (v0[1] > Ymax) {
				Ymax = v0[1];
			}
			if (v0[1] < Ymin) {
				Ymin = v0[1];
			}
			if (v1[1] > Ymax) {
				Ymax = v1[1];
			}
			if (v1[1] < Ymin) {
				Ymin = v1[1];
			}
			if (v2[1] > Ymax) {
				Ymax = v2[1];
			}
			if (v2[1] < Ymin) {
				Ymin = v2[1];
			}

			if (v0[2] > Zmax) {
				Zmax = v0[2];
			}
			if (v0[2] < Zmin) {
				Zmin = v0[2];
			}
			if (v1[2] > Zmax) {
				Zmax = v1[2];
			}
			if (v1[2] < Zmin) {
				Zmin = v1[2];
			}
			if (v2[2] > Zmax) {
				Zmax = v2[2];
			}
			if (v2[2] < Zmin) {
				Zmin = v2[2];
			}
			
		
	}
}

bool intersectBox(const Vec3Df &origin, const Vec3Df &direction, float &Xmax, float &Xmin, float&Ymax, float&Ymin, float&Zmax, float &Zmin)
{
	//t = (xmin - Ox) / dxmin
	float txmin = (Xmin - origin[0]) / direction[0];
	
	float txmax = (Xmax - origin[0]) / direction[0];

	float tymin = (Ymin - origin[1]) / direction[1];

	float tymax = (Ymax - origin[1]) / direction[1];

	float tzmin = (Zmin - origin[2]) / direction[2];

	float tzmax = (Zmax - origin[2]) / direction[2];

	float tinx = min(txmin, txmax);
	float toutx = max(txmin, txmax);
	float tiny = min(tymin, tymax);
	float touty = max(tymin, tymax);
	float tinz = min(tzmin, tzmax);
	float toutz = max(tzmin, tzmax);

	float tin = max(max(tinx, tiny),max(tiny, tinz));
	float tout = min(min(toutx, touty), min(touty, toutz));

	if (tin > tout) {
		return false;
	}
	if (tout < 0) {
		return false;
	}
	std::cout << "Tin:" << tin << " Tout: " << tout;
	return true;
}



bool intersectPlane(const Vec3Df &normal, const Vec3Df &direction, const Vec3Df &origin, float &distance, float &t, Vec3Df &planepos)
{
	// t = (dist - dot(orig, normal) / dot(direction, normal)
	float denominator = (distance - Vec3Df::dotProduct(origin, normal));
	if (Vec3Df::dotProduct(direction,normal) != 0) {
		float numerator = Vec3Df::dotProduct(direction, normal);
		t = denominator / numerator;
		planepos = origin + t * direction;
		if (t >= 0) {
			return true;
		}
	}
	return false;
}

bool rayTriangleIntersect(Vec3Df &planepos, Triangle &triangle, Vec3Df &trianglepos, Vec3Df &normal) {

	Vec3Df bary;
	Vec3Df a = MyMesh.vertices[triangle.v[0]].p;
	Vec3Df b = MyMesh.vertices[triangle.v[1]].p;
	Vec3Df c = MyMesh.vertices[triangle.v[2]].p;

	float areaABC = Vec3Df::dotProduct(normal, Vec3Df::crossProduct((b - a), (c - a)));
	float areaPBC = Vec3Df::dotProduct(normal, Vec3Df::crossProduct((b - planepos), (c - planepos)));
	float areaPCA = Vec3Df::dotProduct(normal, Vec3Df::crossProduct((c - planepos), (a - planepos)));

	bary[0] = areaPBC / areaABC;
	bary[1] = areaPCA / areaABC;
	bary[2] = 1 - bary[0] - bary[1];

	if ((bary[0] < 0) || (bary[0] > 1)) {
		return false;
	}
	if ((bary[1] < 0) || (bary[1] > 1)) {
		return false;
	}
	if (bary[0] + bary[1] > 1) {
		return false;
	}
	else {
		return true;
	}
}



Vec3Df getTriangleCenter(const Vec3Df &edge1, const Vec3Df &edge2, const Vec3Df &edge3) {
	Vec3Df centerOfTriangle = (edge1 + edge2 + edge3 / 3);
	return centerOfTriangle;
}

bool  Shade(Vec3Df shadowOrig, Vec3Df shadowDest, Vec3Df normal, float t, Vec3Df planepos, float distance, Triangle currenttriangle) {
	Vec3Df direction = shadowDest - shadowOrig;
	if (Vec3Df::dotProduct(normal,shadowDest)) {
		printf("shadedN");
		return true;
	}
	if (intersectPlane(normal, direction, shadowOrig, distance, t, planepos)) {
		Vec3Df trianglepos;
		if (rayTriangleIntersect(planepos, currenttriangle, trianglepos, normal)) {
			printf("shadedT");
			return true;
		}
	}
	return false;
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	Vec3Df direction = dest - origin;
	float Xmax; float Xmin; float Ymax; float Ymin; float Zmax; float Zmin;
	calculateMainBox(Xmax, Xmin, Ymax, Ymin, Zmax, Zmin);
	std::cout << "Xmax: " << Xmax << " Xmin:  " << Xmin << " Ymax: " << Ymax << " Ymin: " << Ymin << " Zmax: " << Zmax << " Zmin: " << Zmin;

	if (intersectBox(origin, direction, Xmax, Xmin, Ymax, Ymin, Zmax, Zmin)) {
		std::cout << "Box hit!";
			for (unsigned int i = 0; i < MyMesh.triangles.size(); i++) {
				Triangle currenttriangle = MyMesh.triangles[i];

				Vec3Df v0 = MyMesh.vertices[currenttriangle.v[0]].p;
				Vec3Df v1 = MyMesh.vertices[currenttriangle.v[1]].p;
				Vec3Df v2 = MyMesh.vertices[currenttriangle.v[2]].p;

				Vec3Df edge12 = v0 - v1;
				Vec3Df edge13 = v0 - v2;
				Vec3Df normal = Vec3Df::crossProduct(edge12, edge13);
				normal.normalize();

				float distance = Vec3Df::dotProduct(normal, v0);
				Vec3Df planepos;
				float t;

				if (intersectPlane(normal, direction, origin, distance, t, planepos)) {
					Vec3Df trianglepos;
					if (rayTriangleIntersect(planepos, currenttriangle, trianglepos, normal)) {
						std::cout << "Triangle hit!";

					}

				}
				for (int i = 0; i < MyLightPositions.size(); i++) {
						Vec3Df shadowDest = MyLightPositions[i].p;	//light source
						Vec3Df shadowOrig = getTriangleCenter(v0, v1, v2);	//triangle center
						Shade(shadowOrig, shadowDest, normal, t, planepos, distance, currenttriangle);

				}


			}
	}

	return Vec3Df(dest[0], dest[1], dest[2]);
}

void yourDebugDraw()
{
	//draw open gl debug stuff
	//this function is called every frame
	//let's draw the mesh
	MyMesh.draw();
	//let's draw the lights in the scene as points
	glPushAttrib(GL_ALL_ATTRIB_BITS); //store all GL attributes
	glDisable(GL_LIGHTING);
	glColor3f(1,1,1);
	glPointSize(10);
	glBegin(GL_POINTS);
	for (int i=0;i<MyLightPositions.size();++i)
		glVertex3fv(MyLightPositions[i].pointer());
	glEnd();
	glPopAttrib();//restore all GL attributes
	//The Attrib commands maintain the state. 
	//e.g., even though inside the two calls, we set
	//the color to white, it will be reset to the previous 
	//state after the pop.
	//as an example: we draw the test ray, which is set by the keyboard function
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3f(0,1,1);
	glVertex3f(testRayOrigin[0], testRayOrigin[1], testRayOrigin[2]);
	glColor3f(0,0,1);
	glVertex3f(testRayDestination[0], testRayDestination[1], testRayDestination[2]);
	glEnd();
	glPointSize(10);
	glBegin(GL_POINTS);
	glVertex3fv(MyLightPositions[0].pointer());
	glEnd();
	glPopAttrib();
	//draw whatever else you want...
	////glutSolidSphere(1,10,10);
	////allows you to draw a sphere at the origin.
	////using a glTranslate, it can be shifted to whereever you want
	////if you produce a sphere renderer, this 
	////triangulated sphere is nice for the preview
}


//yourKeyboardFunc is used to deal with keyboard input.
//t is the character that was pressed
//x,y is the mouse position in pixels
//rayOrigin, rayDestination is the ray that is going in the view direction UNDERNEATH your mouse position.
//
//A few keys are already reserved: 
//'L' adds a light positioned at the camera location to the MyLightPositions vector
//'l' modifies the last added light to the current 
//    camera position (by default, there is only one light, so move it with l)
//    ATTENTION These lights do NOT affect the real-time rendering. 
//    You should use them for the raytracing.
//'r' calls the function performRaytracing on EVERY pixel, using the correct associated ray. 
//    It then stores the result in an image "result.ppm".
//    Initially, this function is fast (performRaytracing simply returns 
//    the target of the ray - see the code above), but once you replaced 
//    this function and raytracing is in place, it might take a 
//    while to complete...
void yourKeyboardFunc(char t, int x, int y, const Vec3Df & rayOrigin, const Vec3Df & rayDestination)
{

	//here, as an example, I use the ray to fill in the values for my upper global ray variable
	//I use these variables in the debugDraw function to draw the corresponding ray.
	//try it: Press a key, move the camera, see the ray that was launched as a line.
	testRayOrigin=rayOrigin;	
	testRayDestination=rayDestination;
	performRayTracing(testRayOrigin, testRayDestination);
	// do here, whatever you want with the keyboard input t.
	
	//...
	
	std::cout<<t<<" pressed! The mouse was in location "<<x<<","<<y<<"!"<<std::endl;	
}

Vec3Df phongDiffuse(const Vec3Df & vertexPos, Vec3Df & normal, const Vec3Df & lightPos, unsigned int index)
{
	Vec3Df l = lightPos - vertexPos;
	return Kd.at(index) * Vec3Df::dotProduct(normal, l);
}

Vec3Df phongSpecular(const Vec3Df & vertexPos, Vec3Df & normal, const Vec3Df & lightPos, const Vec3Df & cameraPos, unsigned int index)
{
	Vec3Df spec;
	if (normal.dotProduct(normal, cameraPos) / (normal.getLength() * normal.getLength()) > 0) {
		Vec3Df view = cameraPos - vertexPos;
		Vec3Df incidence = vertexPos - lightPos;
		Vec3Df reflection = incidence - (2 * (Vec3Df::dotProduct(incidence, normal)) * normal);
		spec = Ks[index] * (std::pow(Vec3Df::dotProduct(view, reflection), Shininess.at(index)));
	}
	return spec;
}

Vec3Df phongModel(const Vec3Df & vertexPos, Vec3Df & normal, const Vec3Df & lightPos, const Vec3Df & cameraPos, unsigned int index)
{
	Vec3Df diff = phongDiffuse(vertexPos, normal, lightPos,index);
	Vec3Df spec = phongSpecular(vertexPos, normal, lightPos, cameraPos, index);
	return diff + spec;
}