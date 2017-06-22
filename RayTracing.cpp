#include "stdafx.h"
#include <stdio.h>
#ifdef WIN32
#include <windows.h>
#endif
#include <GL/glut.h>
#include "raytracing.h"


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
    MyMesh.loadMesh("dodgeColorTest.obj", true);
	MyMesh.computeVertexNormals();

	//one first move: initialize the first light source
	//at least ONE light source has to be in the scene!!!
	//here, we set it to the current location of the camera
	MyLightPositions.push_back(MyCameraPosition);


	Kd.resize(MyMesh.vertices.size(), Vec3Df(0.5, 0.5, 0.5));
	Ks.resize(MyMesh.vertices.size(), Vec3Df(0.5, 0.5, 0.5));
	Shininess.resize(MyMesh.vertices.size(), 3);
}


//float D = N.x * v0.x + N.y * v0.y + N.z * v0.z; 

bool rayTriangleIntersect(
	const Vec3Df &orig, const Vec3Df &dir,
	const Vec3Df &edge1, const Vec3Df &edge2, const Vec3Df &edge3)
{
	// compute plane's normal
	Vec3Df line12 = edge2 - edge1;
	Vec3Df line31 = edge3 - edge1;
	// normalising is not neccesary
	Vec3Df N = Vec3Df::crossProduct(line12, line31);
	float area = N.getLength();

	// find P

	// find if ray and 
	float dotproductNandRay = Vec3Df::dotProduct(N, dir);

	if (fabs(dotproductNandRay) < FLT_EPSILON) // epsilon is 1E-5 (0.00001)
	{
		return false;
	}
	float d = Vec3Df::dotProduct(N,edge1);

	float t = (Vec3Df::dotProduct(N, orig) + d) / dotproductNandRay;

	if (t < 0) { // check whether or not the triangle is behind the ray
		return false;
	}

	Vec3Df P = orig + t * dir;

	
	// this whole calculation checks whether or not the ray lands inside the triangle from the right side
	Vec3Df diffedge21 = edge2 - edge1;
	Vec3Df vp0 = P - edge1;
	Vec3Df C = Vec3Df::crossProduct(diffedge21, vp0);
	if (Vec3Df::dotProduct(N, C) < 0) {
		return false;
	}// P is on the right side 

												   // edge 1
	Vec3Df diffedge32 = edge3 - edge2;
	Vec3Df vp1 = P - edge2;
	C = Vec3Df::crossProduct(diffedge32,vp1);
	if (Vec3Df::dotProduct(N,C) < 0) { 
		return false; 
	}// P is on the right side 

															   // edge 2
	Vec3Df diffedge13 = edge1 - edge3;
	Vec3Df vp2 = P - edge3;
	C = Vec3Df::crossProduct(diffedge13, vp2);
	if (Vec3Df::dotProduct(N, C) < 0) {
		return false;
	}// P is on the right side;

	printf("hit a triangle");
	return true;

}

bool intersectPlane(const Vec3Df &normal, const Vec3Df &planeCenter, const Vec3Df &direction, const Vec3Df &origin, float &distance)
{
	// t = (dist - dot(orig, normal) / dot(direction, normal)
	float t;
	float denominator = (distance - Vec3Df::dotProduct(origin, normal));
	if (denominator > 0) {
		Vec3Df p0l0 = planeCenter - origin;
		float numerator = Vec3Df::dotProduct(p0l0, normal);
		t = numerator / denominator;
		if (t >= 0) {
			return true;
		}
	}

	return false;

}

Vec3Df getTriangleCenter(const Vec3Df &edge1, const Vec3Df &edge2, const Vec3Df &edge3) {
	Vec3Df centerOfTriangle = (edge1 + edge2 + edge3 / 3);
	return centerOfTriangle;
}

//return the color of your pixel.
Vec3Df performRayTracing(const Vec3Df & origin, const Vec3Df & dest)
{
	printf("test");
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
		Vec3Df direction = dest - origin;

		Vec3Df centerOfTriangle = getTriangleCenter(v0, v1, v2);

		if (intersectPlane(normal,centerOfTriangle,direction,origin,distance)) {
			printf("hit plane");
		}
		else {
			printf("miss");
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