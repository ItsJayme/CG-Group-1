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


//returns the exact center of a triangle by taking the average of the distance of the points
Vec3Df getTriangleCenter(const Vec3Df &edge1, const Vec3Df &edge2, const Vec3Df &edge3) {
	Vec3Df centerOfTriangle = (edge1 + edge2 + edge3 / 3);
	return centerOfTriangle;
}

bool  Shade(Vec3Df shadowOrig, Vec3Df shadowDest, Vec3Df normal, float &t, Vec3Df planepos, float distance, Triangle currenttriangle) {
	Vec3Df direction = shadowDest - shadowOrig;	//direction of the shadow, towards the light
	if (Vec3Df::dotProduct(normal,direction) < FLT_EPSILON) { //test if the triangle can even be hit by light, if not, return true
		printf("shadedN");
		return true;
	}
	if (intersectPlane(normal, direction, shadowOrig, distance, t, planepos)) {//run an intersect plane to see if the shadow ray hits a plane
		Vec3Df trianglepos;
		for (int i = 0; i < MyMesh.triangles.size(); i++) {
		Triangle currenttriangle = MyMesh.triangles[i];	//assign an easier to use name to the triangle
														//get the corners
		Vec3Df v0 = MyMesh.vertices[currenttriangle.v[0]].p;
		Vec3Df v1 = MyMesh.vertices[currenttriangle.v[1]].p;
		Vec3Df v2 = MyMesh.vertices[currenttriangle.v[2]].p;
		//calculate edges to corner 1
		Vec3Df edge12 = v0 - v1;
		Vec3Df edge13 = v0 - v2;
		//take the normal of the edges
		Vec3Df normal = Vec3Df::crossProduct(edge12, edge13);
		normal.normalize();
		//get the origin of the shadow ray
		Vec3Df shadoworig = getTriangleCenter(v0, v1, v2);
		//start iterating through the lighting positions
		if (rayTriangleIntersect(planepos, currenttriangle, trianglepos, normal)) {//test if the shadow ray hits a triangle on said plane
			printf("shadedT");
			return true;
		}
		}
	}
	return false;// return false if it doesnt intersect any triangles, this code may be faulty
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	Vec3Df direction = dest - origin;	//get the direction of the ray to be traced
	float Xmax; float Xmin; float Ymax; float Ymin; float Zmax; float Zmin;
	calculateMainBox(Xmax, Xmin, Ymax, Ymin, Zmax, Zmin);	//create the main bounding box, any rays outside this are not calculated for


	if (intersectBox(origin, direction, Xmax, Xmin, Ymax, Ymin, Zmax, Zmin)) {	//test if the ray hits a bounding box
		for (unsigned int i = 0; i < MyMesh.triangles.size(); i++) {	//if it does, run this code and 
			Triangle currenttriangle = MyMesh.triangles[i];	//assign an easier to use variable to the currently used triangle
			//get the corners of the triangle
			Vec3Df v0 = MyMesh.vertices[currenttriangle.v[0]].p;
			Vec3Df v1 = MyMesh.vertices[currenttriangle.v[1]].p;
			Vec3Df v2 = MyMesh.vertices[currenttriangle.v[2]].p;
			//calculate the edges leading to corner 1 from the other 2 corners
			Vec3Df edge12 = v0 - v1;
			Vec3Df edge13 = v0 - v2;
			//get the normal by performing a crossproduct of the 2 calculated edges from the point
			Vec3Df normal = Vec3Df::crossProduct(edge12, edge13);
			normal.normalize(); 
			
			float distance = Vec3Df::dotProduct(normal, v0); //calculate the distance between the corner and the normal
			Vec3Df planepos;
			float t;

			if (intersectPlane(normal, direction, origin, distance, t, planepos)) {	//test if the ray intersects with a plane
				Vec3Df trianglepos;
				if (rayTriangleIntersect(planepos, currenttriangle, trianglepos, normal)) { //test if it then intersects with a triangle
					return Vec3Df(1, 1, 1);	//assign a white colour if it hits
				}

			}

		}

		for (unsigned int i = 0; i < MyMesh.triangles.size(); i++) { //go through the list of triangles
			Triangle currenttriangle = MyMesh.triangles[i];	//assign an easier to use name to the triangle
			//get the corners
			Vec3Df v0 = MyMesh.vertices[currenttriangle.v[0]].p;
			Vec3Df v1 = MyMesh.vertices[currenttriangle.v[1]].p;
			Vec3Df v2 = MyMesh.vertices[currenttriangle.v[2]].p;
			//calculate edges to corner 1
			Vec3Df edge12 = v0 - v1;
			Vec3Df edge13 = v0 - v2;
			//take the normal of the edges
			Vec3Df normal = Vec3Df::crossProduct(edge12, edge13);
			normal.normalize();
			//get the origin of the shadow ray
			Vec3Df shadoworig = getTriangleCenter(v0, v1, v2);
			//start iterating through the lighting positions
			for (int i = 0; i < MyLightPositions.size(); i++) {
				float t;
				Vec3Df planepos;
				Vec3Df shadowdest = MyLightPositions[i].p;	//destination is the light source

				float shadowdistance = Vec3Df::dotProduct(normal, v0); //calculate the distance between corner 1 and the normal
				if (Shade(shadoworig, shadowdest, normal, t, planepos, shadowdistance, currenttriangle)) { //run it through shading
					return Vec3Df(0, 0, 0);	//return a value of 0.0.0 if it doesnt hit
				}
			}
		}
	}
	return Vec3Df(.1, 1, .1);	//return an arbitrary colour if no bounding boxes are hit
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