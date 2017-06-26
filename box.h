#pragma once

/*
 * Box class
 */

#include <vector>
#include "mesh.h"

class Box {
	std::list<Triangle*> t;
	float values[6];
    std::string name;

	Box* leftChild;
	Box* rightChild;
public:
	inline Box(){};

    inline void setName(std::string n) {
        name = n;
    };
    
    inline std::string getName() {
        return name;
    };    
    
	inline int getNoTriangles() {
		return t.size();
	}

	inline void addTriangle(Triangle* triangle) {
		t.push_back(triangle);
	}

	inline void removeTriangles() {
		t.clear();
	}

	inline std::list<Triangle*> getTriangles() {
		return t;
	}

	inline void setLeftChild(Box* lChild) {
		leftChild = lChild;
	}

	inline Box* getLeftChild() {
		return leftChild;
	}

	inline void setRightChild(Box* rChild) {
		rightChild = rChild;
	}

	inline Box* getRightChild() {
		return rightChild;
	}

	inline float getMinX() {
		return values[0];
	}

	inline void setMinX(float minX) {
		values[0] = minX;
	}

	inline float getMaxX() {
		return values[1];
	}

	inline void setMaxX(float maxX) {
		values[1] = maxX;
	}

	inline float getMinY() {
		return values[2];
	}

	inline void setMinY(float minY) {
		values[2] = minY;
	}

	inline float getMaxY() {
		return values[3];
	}

	inline void setMaxY(float maxY) {
		values[3] = maxY;
	}

	inline float getMinZ() {
		return values[4];
	}

	inline void setMinZ(float minZ) {
		values[4] = minZ;
	}

	inline float getMaxZ() {
		return values[5];
	}

	inline void setMaxZ(float maxZ) {
		values[5] = maxZ;
	}
};
 
