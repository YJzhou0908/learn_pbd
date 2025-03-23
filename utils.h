#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>  
#include <string>

using namespace std;

 

typedef struct Tetrahedron {
	uint32_t mId;
	uint32_t mNodes[4];
};

class Tetrahedrons {
public:
	void addPositionsByIndex(uint32_t index, glm::dvec3& pos) {
		(*mPositions)[index * 3 + 0] += pos.x;
		(*mPositions)[index * 3 + 1] += pos.y;
		(*mPositions)[index * 3 + 2] += pos.z;
	}

	glm::dvec3 getPositionsByIndex(uint32_t index) {
		double x = (*mPositions)[index * 3 + 0];
		double y = (*mPositions)[index * 3 + 1];
		double z = (*mPositions)[index * 3 + 2];
		
		return glm::dvec3(x, y, z);
	}

public:
	std::vector<Tetrahedron> mTets;
	std::vector<float>* mPositions;
	std::vector<float>* mOldPositions;
	uint32_t mNumber;
};


class Faces {
public:
	std::vector<uint32_t>* mIndices;
	std::vector<float>* mPosition;
	uint32_t mNumber;
};


Tetrahedrons ConstructTetFromFile(const string& elePath, std::vector<float> *positions, std::vector<float>* oldPostions) {
	Tetrahedrons tets;
	std::ifstream file(elePath);
	if (!file.is_open()) {
		std::cerr << "Failed to open file: " << elePath << std::endl;
		return tets;
	}

	string line;
	std::getline(file, line);
	std::istringstream iss(line);
	iss >> tets.mNumber;
	tets.mTets.resize(tets.mNumber);
	tets.mOldPositions = oldPostions;
	tets.mOldPositions->resize(positions->size());
	for (int i = 0; i < tets.mNumber; i++) {
		std::getline(file, line);
		iss.clear();
		iss.str(line);

		iss >> tets.mTets[i].mId;
		iss >> tets.mTets[i].mNodes[0];
		iss >> tets.mTets[i].mNodes[1];
		iss >> tets.mTets[i].mNodes[2];
		iss >> tets.mTets[i].mNodes[3];
	}

	tets.mPositions = positions;
	std::copy(tets.mPositions->begin(), tets.mPositions->end(), tets.mOldPositions->begin());
	return tets;
}


void ConstructFacesFromFile(const string& facePath, std::vector<uint32_t>* indices) {
	int num;
	std::ifstream file(facePath);
	if (!file.is_open()) {
		std::cerr << "Failed to open file: " << facePath << std::endl;
		return;
	}
	string line;
	std::getline(file, line);
	std::istringstream iss(line);
	iss >> num;
	indices->resize(num * 3);
	for (int i = 0; i < num; i++) {

		std::getline(file, line);
		iss.clear();
		iss.str(line);
		int index;
		iss >> index;
		iss >> (*indices)[i * 3 + 0];
		iss >> (*indices)[i * 3 + 1];
		iss >> (*indices)[i * 3 + 2];
	}


	return ;
}


void SetNodePositionsFromFile(const string& nodePath, std::vector<float>* positions) {
	std::ifstream file(nodePath);
	if (!file.is_open()) {
		std::cerr << "Failed to open file: " << nodePath << std::endl;
		return;
	}

	string line;
	std::getline(file, line);
	std::istringstream iss(line);
	int nNode;
	iss >> nNode;
	positions->resize(nNode * 3);
	for (int i = 0; i < nNode; i++) {
		std::getline(file, line);
		iss.clear();
		iss.str(line);
		int index;
		iss >> index;
		iss >> (*positions)[i * 3 + 0] ;
		iss >> (*positions)[i * 3 + 1] ;
		iss >> (*positions)[i * 3 + 2] ;
	}

}


