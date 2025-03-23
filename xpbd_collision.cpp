#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>  
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>
#include "utils.h"



using namespace std;

GLFWwindow* window;
uint32_t WIDTH = 800;
uint32_t HEIGHT = 800;
uint32_t ITER = 10;
float dt = 0.0033f;
float E = 1e2f, nu = 0.4;
float mu = E / (2 * (1 + nu));  // Young's modulus and Poisson's ratio
float la = E * nu / ((1 + nu) * (1 - 2 * nu)); // Lame parameters

static const float MODE_COMPLIANCE[] = {
  0.0f,            // Miles Macklin's blog (http://blog.mmacklin.com/2016/10/12/xpbd-slides-and-stiffness/)
  0.00000000004f, // 0.04 x 10^(-9) (M^2/N) Concrete
  0.00000000016f, // 0.16 x 10^(-9) (M^2/N) Wood
  0.000000001f,   // 1.0  x 10^(-8) (M^2/N) Leather
  0.000000002f,   // 0.2  x 10^(-7) (M^2/N) Tendon
  0.0000001f,     // 1.0  x 10^(-6) (M^2/N) Rubber
  0.00002f,       // 0.2  x 10^(-3) (M^2/N) Muscle
  0.0001f,        // 1.0  x 10^(-3) (M^2/N) Fat
};


string readShaderFile(const char* filename) {
	std::ifstream file(filename);
	if (!file) {
		std::cerr << "无法打开文件: " << filename << '\n';
		return nullptr;
	}

	std::stringstream buffer;
	buffer << file.rdbuf();
	std::string contents(buffer.str());
	file.close();
	return contents;
}

class Camera {
public:
	Camera(glm::vec3 position, glm::vec3 lookAt, glm::vec3 up) {
		mPosition = position;
		mViewMatrix = glm::lookAt(position, lookAt, up);
		mProjectionMatrix = glm::perspective<float>(glm::radians(45.0f), float(WIDTH) / float(HEIGHT), 0.1f, 100.0f);
		mProjectionViewMatrix = mProjectionMatrix * mViewMatrix;
	}

	~Camera() { };
	glm::mat4& getProjectionViewMatrix() { return mProjectionViewMatrix; }
	glm::mat4& getProjectionMatrix() { return mProjectionMatrix; }
	glm::mat4& getViewMatrix() { return mViewMatrix; }
private:
	glm::vec3 mPosition;
	glm::mat4 mProjectionViewMatrix;
	glm::mat4 mViewMatrix;
	glm::mat4 mProjectionMatrix;
};


class Constrain {

public:
	Constrain(Tetrahedrons* tets, std::vector <float>* mass_inv, uint32_t id) {

		mTets = tets;
		mDensity = 1e3f;
		mVolCompliance = 0.0f;
		mDevCompliance = 5e-5f;
		mMass_Inv = mass_inv;
		mId = id;

		uint32_t Id0 = mTets->mTets[id].mNodes[0];
		uint32_t Id1 = mTets->mTets[id].mNodes[1];
		uint32_t Id2 = mTets->mTets[id].mNodes[2];
		uint32_t Id3 = mTets->mTets[id].mNodes[3];

		glm::vec3 x0 = mTets->getPositionsByIndex(Id0);
		glm::vec3 x1 = mTets->getPositionsByIndex(Id1);
		glm::vec3 x2 = mTets->getPositionsByIndex(Id2);
		glm::vec3 x3 = mTets->getPositionsByIndex(Id3);


		glm::vec3 firstCol = x1 - x0;
		glm::vec3 secondCol = x2 - x0;
		glm::vec3 thirdCol = x3 - x0;

		glm::mat3 X(
			firstCol.x, firstCol.y, firstCol.z, secondCol.x, secondCol.y, secondCol.z, thirdCol.x, thirdCol.y, thirdCol.z
		);

		mVolume = glm::determinant(X) / 6.0f;

		mX_inv = glm::inverse(X);


		float pm = mVolume * mDensity / 4.0f;

		(*mMass_Inv)[Id0] += pm;
		(*mMass_Inv)[Id1] += pm;
		(*mMass_Inv)[Id2] += pm;
		(*mMass_Inv)[Id3] += pm;

	}

	void lambdaInit() {
		mLambda = 0.0f;
	}

	void solve(float dt) {
		// 计算constrain
		uint32_t Id0 = mTets->mTets[mId].mNodes[0];
		uint32_t Id1 = mTets->mTets[mId].mNodes[1];
		uint32_t Id2 = mTets->mTets[mId].mNodes[2];
		uint32_t Id3 = mTets->mTets[mId].mNodes[3];

		glm::vec3 x0 = mTets->getPositionsByIndex(Id0);
		glm::vec3 x1 = mTets->getPositionsByIndex(Id1);
		glm::vec3 x2 = mTets->getPositionsByIndex(Id2);
		glm::vec3 x3 = mTets->getPositionsByIndex(Id3);

		glm::vec3 firstCol = x1 - x0;
		glm::vec3 secondCol = x2 - x0;
		glm::vec3 thirdCol = x3 - x0;

		glm::mat3 Ds(
			firstCol.x, firstCol.y, firstCol.z, secondCol.x, secondCol.y, secondCol.z, thirdCol.x, thirdCol.y, thirdCol.z
		);


		glm::mat3 F = Ds * mX_inv;
		float detF = glm::determinant(F);
		//cout << detF << endl;
		float CD = std::sqrt(glm::dot(F[0], F[0]) + glm::dot(F[1], F[1]) + glm::dot(F[2], F[2]));
		float CD_inv = 1.0f / CD;
		// 计算grad
		glm::vec3 g[4];
		g[1] = glm::vec3(0.0f);
		g[1] += F[0] * CD_inv * mX_inv[0].x;
		g[1] += F[1] * CD_inv * mX_inv[1].x;
		g[1] += F[2] * CD_inv * mX_inv[2].x; // F * [I,0,0] * X_inv相当于每一列扩大一定的倍数

		g[2] = glm::vec3(0.0f);
		g[2] += F[0] * CD_inv * mX_inv[0].y;
		g[2] += F[1] * CD_inv * mX_inv[1].y;
		g[2] += F[2] * CD_inv * mX_inv[2].y;

		g[3] = glm::vec3(0.0f);
		g[3] += F[0] * CD_inv * mX_inv[0].z;
		g[3] += F[1] * CD_inv * mX_inv[1].z;
		g[3] += F[2] * CD_inv * mX_inv[2].z;

		g[0] = glm::vec3(0.0f);
		g[0] = -g[1] - g[2] - g[3];

		// sum
		float mass_inv[] = {
			(*mMass_Inv)[Id0],
			(*mMass_Inv)[Id1],
			(*mMass_Inv)[Id2],
			(*mMass_Inv)[Id3]
		};
		float w = 0.0f;
		w += glm::dot(g[0], g[0]) * mass_inv[0];
		w += glm::dot(g[1], g[1]) * mass_inv[1];
		w += glm::dot(g[2], g[2]) * mass_inv[2];
		w += glm::dot(g[3], g[3]) * mass_inv[3];

		float alpha = mDevCompliance / (dt * dt * mVolume);
		float dlambda = -CD / (w + alpha);

		mTets->addPositionsByIndex(Id0, g[0] * dlambda * mass_inv[0]);
		mTets->addPositionsByIndex(Id1, g[1] * dlambda * mass_inv[1]);
		mTets->addPositionsByIndex(Id2, g[2] * dlambda * mass_inv[2]);
		mTets->addPositionsByIndex(Id3, g[3] * dlambda * mass_inv[3]);
		
		x0 = mTets->getPositionsByIndex(Id0);
		x1 = mTets->getPositionsByIndex(Id1);
		x2 = mTets->getPositionsByIndex(Id2);
		x3 = mTets->getPositionsByIndex(Id3);

		firstCol = x1 - x0;
		secondCol = x2 - x0;
		thirdCol = x3 - x0;

		glm::mat3 Drs(
			firstCol.x, firstCol.y, firstCol.z, secondCol.x, secondCol.y, secondCol.z, thirdCol.x, thirdCol.y, thirdCol.z
		);


		F = Drs * mX_inv;

		glm::mat3 dF;
		dF[0] = glm::cross(F[1], F[2]); // 算出F的梯度
		dF[1] = glm::cross(F[2], F[0]);
		dF[2] = glm::cross(F[0], F[1]);

		g[1] = glm::vec3(0.0f);
		g[1] += dF[0] * mX_inv[0].x;
		g[1] += dF[1] * mX_inv[1].x;
		g[1] += dF[2] * mX_inv[2].x;

		g[2] = glm::vec3(0.0f);
		g[2] += dF[0] * mX_inv[0].y;
		g[2] += dF[1] * mX_inv[1].y;
		g[2] += dF[2] * mX_inv[2].y;

		g[3] = glm::vec3(0.0f);
		g[3] += dF[0] * mX_inv[0].z;
		g[3] += dF[1] * mX_inv[1].z;
		g[3] += dF[2] * mX_inv[2].z;

		float vol = glm::determinant(F); // 求F的行列式 就是体积

		float CH = vol - 1.0f;


		g[0] = -g[1] - g[2] - g[3];


		w = 0.0f;
		w += glm::dot(g[0], g[0]) * mass_inv[0];
		w += glm::dot(g[1], g[1]) * mass_inv[1];
		w += glm::dot(g[2], g[2]) * mass_inv[2];
		w += glm::dot(g[3], g[3]) * mass_inv[3];

		alpha = mVolCompliance / (dt * dt * mVolume);
		dlambda = -CH / (w + alpha);


		mTets->addPositionsByIndex(Id0, g[0] * dlambda * mass_inv[0]);
		mTets->addPositionsByIndex(Id1, g[1] * dlambda * mass_inv[1]);
		mTets->addPositionsByIndex(Id2, g[2] * dlambda * mass_inv[2]);
		mTets->addPositionsByIndex(Id3, g[3] * dlambda * mass_inv[3]);


	}
private:
	Tetrahedrons* mTets;
	glm::mat3 mX_inv;
	uint32_t mId;
	std::vector<float>* mMass_Inv;
	float mLambda;
	float mVolume;
	float mVolCompliance;
	float mDevCompliance;
	float mDensity;


};

class SoftBody {
public:
	SoftBody(
		const std::string& eleFilePath,
		const std::string& nodeFilePath,
		const std::string& faceFilePath,
		glm::vec3 center,
		glm::vec3 velocity
	) {
		SetNodePositionsFromFile(nodeFilePath, &mPositions);
		mTets = ConstructTetFromFile(eleFilePath, &mPositions, &mOldPositions);
		ConstructFacesFromFile(faceFilePath, &mIndices);
		mCenter = center;
		mAccelerate = glm::vec3(0.0f, 0.0f, 0.0f);
		mMass_inv.resize(mPositions.size() / 3);

		for (int i = 0; i < mTets.mNumber; i++) {
			// 一个四面体构建一个约束
			mConstrains.push_back(Constrain(&mTets, &mMass_inv, i));
		}
		mVelocity.resize(mPositions.size());
		std::fill(mVelocity.begin(), mVelocity.end(), 0.0f);

		for (int i = 0; i < mVelocity.size(); i += 3) {
			mVelocity[i] = velocity.x;
			mVelocity[i + 1] = velocity.y;
			mVelocity[i + 2] = velocity.z;
		}
		// do Mass Inv;
		doMassInverse();

		
	}


	void doMassInverse() {
		for (int i = 0; i < mPositions.size() / 3; i++) {
			mMass_inv[i] = 1.0f / mMass_inv[i];
		}
	}


	glm::vec3 getPositionByIndex(uint32_t index) {
		glm::vec3 pos;
		pos.x = mPositions[index * 3 + 0];
		pos.y = mPositions[index * 3 + 1];
		pos.z = mPositions[index * 3 + 2];
		return pos;
	}
	void setPostionByIndex(uint32_t index, glm::vec3& pos) {
		mPositions[index * 3 + 0] = pos.x;
		mPositions[index * 3 + 1] = pos.y;
		mPositions[index * 3 + 2] = pos.z;
	}

	void addPostionByIndex(uint32_t index, glm::vec3& pos) {
		mPositions[index * 3 + 0] += pos.x;
		mPositions[index * 3 + 1] += pos.y;
		mPositions[index * 3 + 2] += pos.z;
	}

	float getMassInvByIndex(uint32_t index) {
		return mMass_inv[index];
	}
	void addMassInvByIndex(uint32_t index, float mass_inv) {
		mMass_inv[index] += mass_inv;
	}


	void renderInit(Camera* camera) {

		mCamera = camera;

		mModelMatrix = glm::mat4(1.0f);

		mModelMatrix = glm::scale(mModelMatrix, glm::vec3(0.5f, 0.5f, 0.5f));

		//mModelMatrix = glm::rotate(mModelMatrix, glm::radians(-90.0f),glm::vec3(1.0f,0.0f,0.0f));

		mModelMatrix = glm::translate(mModelMatrix, mCenter);

		// 创建 VAO（顶点数组对象）
		glGenVertexArrays(1, &mVAO);
		glBindVertexArray(mVAO); // 必须先创建VAO，然后绑定VAO，绑定好VAO之后才能在VAO中写入VBO

		// 创建 VBO（顶点缓冲对象）
		glGenBuffers(1, &mPosVbo);
		glBindBuffer(GL_ARRAY_BUFFER, mPosVbo);
		glBufferData(GL_ARRAY_BUFFER, mPositions.size() * sizeof(float), mPositions.data(), GL_STATIC_DRAW);


		// 配置顶点属性（位置属性）
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(0);



		// 创建 EBO（索引缓冲对象）
		glGenBuffers(1, &mEBO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mEBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, mIndices.size() * sizeof(unsigned int), mIndices.data(), GL_STATIC_DRAW);


		// 解绑 VAO 和 VBO
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);


		// 处理Shader
		GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
		GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
		string vertexShaderCode = readShaderFile("./shader/basic.vert");
		string fragShaderCode = readShaderFile("./shader/basic.frag");


		const char* vertexShaderSource = vertexShaderCode.c_str();
		const char* fragShaderSource = fragShaderCode.c_str();

		glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
		glCompileShader(vertexShader);

		int success;
		char infoLog[512];
		glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
		if (!success) {
			glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
			std::cout << "ERROR::VERTEX_SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
		}

		glShaderSource(fragShader, 1, &fragShaderSource, NULL);
		glCompileShader(fragShader);

		glGetShaderiv(fragShader, GL_COMPILE_STATUS, &success);
		if (!success) {
			glGetShaderInfoLog(fragShader, 512, NULL, infoLog);
			std::cout << "ERROR::FRAGMENT_SHADER::COMPILATION_FAILED\n" << infoLog << std::endl;
		}

		mShaderProgram = glCreateProgram();
		glAttachShader(mShaderProgram, vertexShader);
		glAttachShader(mShaderProgram, fragShader);
		glLinkProgram(mShaderProgram);

		// 检查链接状态
		glGetProgramiv(mShaderProgram, GL_LINK_STATUS, &success);
		if (!success) {
			glGetProgramInfoLog(mShaderProgram, 512, NULL, infoLog);
			std::cout << "ERROR::SHADER_PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
		}

		// 删除着色器（因为它们已经链接到程序中）
		glDeleteShader(vertexShader);
		glDeleteShader(fragShader);

		mProjectionViewLoc = glGetUniformLocation(mShaderProgram, "projectionViewMatrix");
		mModelLoc = glGetUniformLocation(mShaderProgram, "modelMatrix");
		mColorLoc = glGetUniformLocation(mShaderProgram, "color");

	}

	void display() {

		glUseProgram(mShaderProgram);

		// projectionViewMatrix
		glUniformMatrix4fv(mProjectionViewLoc, 1, GL_FALSE, glm::value_ptr(mCamera->getProjectionViewMatrix()));

		// modelMatrix
		glUniformMatrix4fv(mModelLoc, 1, GL_FALSE, glm::value_ptr(mModelMatrix));


		glUniform3f(mColorLoc, 1.0f, 1.0f, 1.0f); // 白色面
		glBindVertexArray(mVAO);

		glBindBuffer(GL_ARRAY_BUFFER, mPosVbo);
		glBufferSubData(GL_ARRAY_BUFFER, 0, mPositions.size() * sizeof(float), mPositions.data());

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, mIndices.size(), GL_UNSIGNED_INT, 0);

		glBindVertexArray(0);

		glUseProgram(0);

	}

	void setOtherBall(SoftBody* other) {
		mOtherBall = other;
	}

	void update(float dt, uint32_t iteration) {
		float radius = 3.0f;
		for (int i = 0; i < mPositions.size(); i += 3) {
			glm::vec3 pos(mPositions[i], mPositions[i + 1], mPositions[i + 2]);

			glm::vec3 p_minus_center = pos - mOtherBall->mCenter;
			float len = glm::length(p_minus_center);
			if (len < radius) {
				mVelocity[i] = -mVelocity[i];
				mVelocity[i + 1]  = -mVelocity[i + 1];
				mVelocity[i + 2] = -mVelocity[i + 2];
				glm::vec3 v = glm::normalize(p_minus_center);
				mPositions[i] += v.x * (len - radius);
				mPositions[i + 1] += v.y * (len - radius);
				mPositions[i + 2] += v.z * (len - radius);
				
			}
		
		}


		for (int i = 0; i < mPositions.size(); i += 3) {
			mPositions[i] += mVelocity[i] * dt;
			mPositions[i + 1] += mVelocity[i + 1] * dt;
			mPositions[i + 2] += mVelocity[i + 2] * dt; 
		}
		mCenter += glm::vec3(mVelocity[0]*dt, mVelocity[1] * dt, mVelocity[2] * dt);


		for (int i = 0; i < mConstrains.size(); i++) {
			mConstrains[i].lambdaInit();
		}

		for (int i = 0; i < iteration; i++) {
			for (int i = 0; i < mConstrains.size(); i ++)  {
				mConstrains[i].solve(dt);
			}
			
		}

	}
public:
	glm::vec3 mCenter;
private:
	void MoveToCenter() {

		int nPoints = mPositions.size() / 3;

		float cx = 0.0f, cy = 0.0f, cz = 0.0f;
		for (int i = 0; i < nPoints; i++) {
			int idx = i * 3;
			cx += mPositions[idx];
			cy += mPositions[idx + 1];
			cz += mPositions[idx + 2];
		}
		cx /= nPoints;
		cy /= nPoints;
		cz /= nPoints;

		float offsetX = mCenter.x - cx;
		float offsetY = mCenter.y - cy;
		float offsetZ = mCenter.z - cz;

		for (int i = 0; i < nPoints; ++i) {
			int idx = i * 3;
			mPositions[idx] += offsetX;
			mPositions[idx + 1] += offsetY;
			mPositions[idx + 2] += offsetZ;
		}
	}
	
private:
	//glm::vec3 mCenter;

	std::vector<float> mOldPositions;
	std::vector<float> mPositions;
	std::vector<uint32_t> mIndices;
	std::vector<float> mColors;
	std::vector<float> mVelocity;
	std::vector <float> mMass_inv;
	glm::vec3 mAccelerate;
	std::vector<Constrain> mConstrains;


	Tetrahedrons mTets;

	Camera* mCamera;
	GLuint mVAO, mPosVbo, mColVbo, mEBO;
	GLuint mShaderProgram;
	glm::mat4 mModelMatrix;
	GLuint mProjectionViewLoc;
	GLuint mModelLoc;
	GLuint mColorLoc;

	SoftBody* mOtherBall;

};


SoftBody softbody1("GOLF.ele", "GOLF.node", "GOLF.face", glm::vec3(5.0f, 0.0f, 0.0f), glm::vec3(-8.0f,0.0f,0.0f));
SoftBody softbody2("GOLF.ele", "GOLF.node", "GOLF.face", glm::vec3(-5.0f, 0.0f, 0.0f), glm::vec3(8.0f, 0.0f, 0.0f));
Camera camera(glm::vec3(0.0f, 15.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 1.0f));




void init() {

	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	//1.2 设置OpenGL启用核心模式（非立即渲染模式）
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	window = glfwCreateWindow(WIDTH, HEIGHT, "xpbd", NULL, NULL);
	if (window == NULL) {

	}
	//设置当前窗体对象为OpenGL的绘制舞台
	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cout << "Failed to initialize GLAD." << std::endl;
	}

	softbody1.renderInit(&camera);
	softbody1.setOtherBall(&softbody2);
	softbody2.renderInit(&camera);
	softbody2.setOtherBall(&softbody1);
}


void display() {

	softbody1.display();
	softbody2.display();
}

void simulation() {

	softbody1.update(dt, ITER);
	softbody2.update(dt, ITER);
}

int main() {

	init();

	while (!glfwWindowShouldClose(window)) {
		glEnable(GL_CULL_FACE);
		glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		for (int i = 0; i < 3; i++)
			simulation();

		display();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	return 0;


}