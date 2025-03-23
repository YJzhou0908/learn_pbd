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



using namespace std;

GLFWwindow* window;
uint32_t WIDTH = 800;
uint32_t HEIGHT = 800;
uint32_t ITER = 30;
float dt = 0.033f;
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


class Ball {
public:
	class Constrain {
	public:
		Constrain(int pos0Id, int pos1Id, int pos2Id, Ball* ball, float density) {
			mIndex[0] = pos0Id;
			mIndex[1] = pos1Id;
			mIndex[2] = pos2Id;
			mBall = ball;
			mDensity = density;
			
			mGrad[0] = glm::vec3(0.0f);
			mGrad[1] = glm::vec3(0.0f);
			mGrad[2] = glm::vec3(0.0f);
			mGrad[3] = glm::vec3(0.0f);
			
			glm::vec3 x0 = mBall->getBallCenter();
			glm::vec3 x1 = mBall->getPositionByIndex(pos0Id);
			glm::vec3 x2 = mBall->getPositionByIndex(pos1Id);
			glm::vec3 x3 = mBall->getPositionByIndex(pos2Id);
	
			
			glm::vec3 firstCol = x1 - x0;
			glm::vec3 secondCol = x2 - x0;
			glm::vec3 thirdCol = x3 - x0;

			glm::mat3 X(
				firstCol.x, firstCol.y, firstCol.z, secondCol.x, secondCol.y, secondCol.z, thirdCol.x,thirdCol.y,  thirdCol.z
			);

			mVolume = glm::determinant(X) / 6.0f;

			mX_inv = glm::inverse(X);

			mX_inv_firstline = glm::vec3(mX_inv[0][0], mX_inv[1][0], mX_inv[2][0]);
			mX_inv_secondline = glm::vec3(mX_inv[0][1], mX_inv[1][1], mX_inv[2][1]);
			mX_inv_thirdline = glm::vec3(mX_inv[0][2], mX_inv[1][2], mX_inv[2][2]);

			float pm = mVolume / 4.0f * mDensity;

			mBall->addMassInvByIndex(mIndex[0], pm);
			mBall->addMassInvByIndex(mIndex[1], pm);
			mBall->addMassInvByIndex(mIndex[2], pm);
			mBall->addMassInvCenter(pm);
				
	
		}

		void lambdaInit() {
			mLambda = 0.0f;
		}

		void solve(float dt) {
			
			
			mDevCompliance = 1e-5f;
			mVolCompliance = 0.0f;
			// 计算constrain
			glm::vec3 x0 = mBall->getBallCenter();
			glm::vec3 x1 = mBall->getPositionByIndex(mIndex[0]);
			glm::vec3 x2 = mBall->getPositionByIndex(mIndex[1]);
			glm::vec3 x3 = mBall->getPositionByIndex(mIndex[2]);
			

			glm::vec3 firstCol = x1 - x0;
			glm::vec3 secondCol = x2 - x0;
			glm::vec3 thirdCol = x3 - x0;

			glm::mat3 Ds(
				firstCol.x, firstCol.y, firstCol.z, secondCol.x, secondCol.y, secondCol.z, thirdCol.x, thirdCol.y, thirdCol.z
			);


			glm::mat3 F = Ds * mX_inv;
			float detF = glm::determinant(F);
			float CD = std::sqrt(glm::dot(F[0], F[0]) + glm::dot(F[1], F[1]) + glm::dot(F[2], F[2]));
			float CD_inv = 1.0f / CD;
			// 计算grad
			mGrad[1] = glm::vec3(0.0f);
			mGrad[1] += F[0] * CD_inv * mX_inv[0].x;
			mGrad[1] += F[1] * CD_inv * mX_inv[1].x;
			mGrad[1] += F[2] * CD_inv * mX_inv[2].x; // F * [I,0,0] * X_inv相当于每一列扩大一定的倍数

			mGrad[2] = glm::vec3(0.0f);
			mGrad[2] += F[0] * CD_inv * mX_inv[0].y;
			mGrad[2] += F[1] * CD_inv * mX_inv[1].y;
			mGrad[2] += F[2] * CD_inv * mX_inv[2].y;

			mGrad[3] = glm::vec3(0.0f);
			mGrad[3] += F[0] * CD_inv * mX_inv[0].z;
			mGrad[3] += F[1] * CD_inv * mX_inv[1].z;
			mGrad[3] += F[2] * CD_inv * mX_inv[2].z;

			mGrad[0] = glm::vec3(0.0f);
			mGrad[0] = -mGrad[1] - mGrad[2] - mGrad[3];

			// sum
			float mass_inv[] = { 
				mBall->getMassInvCenter(),
				mBall->getMassInvByIndex(mIndex[0]),
				mBall->getMassInvByIndex(mIndex[1]),
				mBall->getMassInvByIndex(mIndex[2])
			};
			float w = 0.0f;
			w += glm::dot(mGrad[0], mGrad[0]) * mass_inv[0];
			w += glm::dot(mGrad[1], mGrad[1]) * mass_inv[1];
			w += glm::dot(mGrad[2], mGrad[2]) * mass_inv[2];
			w += glm::dot(mGrad[3], mGrad[3]) * mass_inv[3];

			float alpha = mDevCompliance / (dt * dt * mVolume);
			float dlambda = -CD / (w + alpha);

			mBall->addPostionByIndex(mIndex[0], mGrad[1] * dlambda * mass_inv[1]);
			mBall->addPostionByIndex(mIndex[1], mGrad[2] * dlambda * mass_inv[2]);
			mBall->addPostionByIndex(mIndex[2], mGrad[3] * dlambda * mass_inv[3]);
			mBall->addCenterPosition(mGrad[0] * dlambda * mass_inv[0]);

			x0 = mBall->getBallCenter();
			x1 = mBall->getPositionByIndex(mIndex[0]);
			x2 = mBall->getPositionByIndex(mIndex[1]);
			x3 = mBall->getPositionByIndex(mIndex[2]);

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

			mGrad[1] = glm::vec3(0.0f);
			mGrad[1] += dF[0] * mX_inv[0].x;
			mGrad[1] += dF[1] * mX_inv[1].x;
			mGrad[1] += dF[2] * mX_inv[2].x;

			mGrad[2] = glm::vec3(0.0f);
			mGrad[2] += dF[0] * mX_inv[0].y;
			mGrad[2] += dF[1] * mX_inv[1].y;
			mGrad[2] += dF[2] * mX_inv[2].y;

			mGrad[3] = glm::vec3(0.0f);
			mGrad[3] += dF[0] * mX_inv[0].z;
			mGrad[3] += dF[1] * mX_inv[1].z;
			mGrad[3] += dF[2] * mX_inv[2].z;

			float vol = glm::determinant(F); // 求F的行列式 就是体积

			float CH = vol - 1.0f ;


			mGrad[0] = -mGrad[1] - mGrad[2] - mGrad[3];

			
			w = 0.0f;
			w += glm::dot(mGrad[0], mGrad[0]) * mass_inv[0];
			w += glm::dot(mGrad[1], mGrad[1]) * mass_inv[1];
			w += glm::dot(mGrad[2], mGrad[2]) * mass_inv[2];
			w += glm::dot(mGrad[3], mGrad[3]) * mass_inv[3];

			alpha = mVolCompliance / (dt * dt * mVolume);
			dlambda = -CH / (w + alpha);
			
			mBall->addPostionByIndex(mIndex[0], mGrad[1] * dlambda* mass_inv[1]);
			mBall->addPostionByIndex(mIndex[1], mGrad[2] * dlambda* mass_inv[2]);
			mBall->addPostionByIndex(mIndex[2], mGrad[3] * dlambda* mass_inv[3]);
			mBall->addCenterPosition(mGrad[0] * dlambda* mass_inv[0]);
			
		}
	public:
		uint32_t mIndex[3]; // 三角形的三个顶点索引
		float mLambda; // 一个约束一次性算出两个mLamda，两个矫正向量
		float mVolume;
		glm::mat3 mX_inv;
		glm::vec3 mX_inv_firstline;
		glm::vec3 mX_inv_secondline;
		glm::vec3 mX_inv_thirdline;
		float mDevCompliance;
		float mVolCompliance;
		glm::vec3 mGrad[4];

		Ball* mBall;
		float mDensity;

	};

public:
	Ball(float radius, glm::vec3& center, int resolution) {
		mRadius = radius;
		mCenter = center;
		mResolution = resolution;
		mCenterMass_inv = 0.0;
		createIcosphere();


		for (int i = 0; i < mIndices.size(); i += 3) {
			mConstrains.push_back(Constrain(mIndices[i], mIndices[i + 1], mIndices[i + 2], this, 1000.0f));
		}
	}

	~Ball() {};
	
	void doMassInverse() {
		for (int i = 0; i < mPositions.size() / 3; i ++) {
			mMass_inv[i] = 1.0f / mMass_inv[i];
		}
		mCenterMass_inv = 1.0f / mCenterMass_inv;
	}

	

	void setCenterMassInv(float mass_inv) {
		mCenterMass_inv = mass_inv;
	}


	glm::vec3 getPositionByIndex(uint32_t index) {
		glm::vec3 pos;
		pos.x = mPositions[index * 3 + 0];
		pos.y = mPositions[index * 3 + 1];
		pos.z = mPositions[index * 3 + 2];
		return pos;
	}

	glm::vec3 getBallCenter() {
		return mCenter;
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

	void addMassInvCenter(float mass_inv) {
		mCenterMass_inv += mass_inv;
	}


	float getMassInvCenter() {
		return mCenterMass_inv;
	}

	void addCenterPosition(glm::vec3 pos) {
		mCenter += pos;
	}

	void renderInit(Camera* camera) {

		mCamera = camera;

		mModelMatrix = glm::mat4(1.0f);

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

	void update(float sdt, uint32_t iteration) {
		for (int i = 0; i < mPositions.size(); i += 3) {
			mVelocity[i] += mAccelerate.x * sdt;
			mVelocity[i + 1] += mAccelerate.y * sdt;
			mVelocity[i + 2] += mAccelerate.z * sdt;
			mOldPositions[i] = mPositions[i];
			mOldPositions[i + 1] = mPositions[i + 1];
			mOldPositions[i + 2] = mPositions[i + 2];
			float dx = mVelocity[i + 1] * sdt + mAccelerate.y * sdt * sdt;
			mPositions[i] += mVelocity[i] * sdt + mAccelerate.x * sdt * sdt; //glm::vec3(0.0f, -0.01f, 0.0f);//
			mPositions[i + 1] += mVelocity[i + 1] * sdt + mAccelerate.y * sdt * sdt;
			mPositions[i + 2] += mVelocity[i + 2] * sdt + mAccelerate.z * sdt * sdt;
		}
		mVelCenter = mAccelerate * dt;
		mOldCenter =  mCenter;
		mCenter += mVelCenter * dt + mAccelerate * dt * dt;


		for (int i = 0; i < mConstrains.size(); i++) {
			mConstrains[i].lambdaInit();
			doMassInverse();
		}

		for (int i = 0; i < iteration; i++) {
			for (int i = 0; i < mPositions.size(); i += 3) {
				if (mPositions[i + 1] < -5.0f) {
					mPositions[i + 1] = -5.0f;
				}
			}
			
			for (auto constraint = mConstrains.begin(); constraint != mConstrains.end(); constraint++) {
				(*constraint).solve(dt);
			}
			
		}
		//cout << mPositions[3] << " " << mPositions[4] << " " << mPositions[5] << endl;

	}
private:
	void createIcosphere() {
		mAccelerate = glm::vec3(0.0f, -9.8f, 0.0f);
		float phi = (1.0f + std::sqrt(5.0f)) / 2.0f;

		mPositions.insert(mPositions.end(), {
			-1.0f, phi, 0.0f, 1.0f, phi, 0.0f,-1.0f, -phi, 0.0f, 1.0f, -phi, 0.0f,
			 0.0f,-1.0f, phi, 0.0f, 1.0f, phi, 0.0f, -1.0f, -phi, 0.0f, 1.0f, -phi,
			 phi, 0.0f, -1.0f, phi, 0.0f, 1.0f, -phi, 0.0f, -1.0f,-phi, 0.0f, 1.0f
		});


		float length = std::sqrt(1.0f + phi * phi) / mRadius;

		for (int i = 0; i < mPositions.size() / 3; i++) {
			mPositions[i * 3 + 0] /= length;
			mPositions[i * 3 + 1] /= length;
			mPositions[i * 3 + 2] /= length;
		}

		std::vector<uint32_t> indices = {
			0, 11, 5,  0, 5, 1,  0, 1, 7,  0, 7, 10,  0, 10, 11,
			1, 5, 9,  5, 11, 4,  11, 10, 2,  10, 7, 6,  7, 1, 8,
			3, 9, 4,  3, 4, 2,  3, 2, 6,  3, 6, 8,  3, 8, 9,
			4, 9, 5,  2, 4, 11, 6, 2, 10, 8, 6, 7,  9, 8, 1
		};


		std::vector<uint32_t> finalIndices = indices;

		for (int i = 0; i < mResolution; i++) {
			finalIndices = divide(finalIndices);
		}

		mIndices = finalIndices;

		for (int i = 0; i < mPositions.size(); i += 3) {
			mPositions[i] += mCenter.x;
			mPositions[i + 1] += mCenter.y;
			mPositions[i + 2] += mCenter.z;
		}

		mOldPositions.resize(mPositions.size());
		std::copy(mPositions.begin(), mPositions.end(), mOldPositions.begin());
		mVelocity.resize(mPositions.size());
		std::fill(mVelocity.begin(), mVelocity.end(), 0.0f);
		mMass_inv.resize(mPositions.size() / 3);
		std::fill(mVelocity.begin(), mVelocity.end(), 0.0f);
		mVelCenter = glm::vec3(0.0f);
		mOldCenter = mCenter;
	}

	std::vector<uint32_t> divide(std::vector<uint32_t>& oldIndices) {
		std::vector<uint32_t> newIndices;

		for (int i = 0; i < oldIndices.size(); i += 3) {
			uint32_t v1 = oldIndices[i];
			uint32_t v2 = oldIndices[i + 1];
			uint32_t v3 = oldIndices[i + 2];

			uint32_t a = getMidPoint(v1, v2);
			uint32_t b = getMidPoint(v2, v3);
			uint32_t c = getMidPoint(v3, v1);


			newIndices.insert(newIndices.end(), {
				v1, a, c,
				a, v2, b,
				a, b, c,
				c, b, v3
				});
		}
		return newIndices;
	}

	uint32_t getMidPoint(int v1, int v2) {
		float x = (mPositions[v1 * 3 + 0] + mPositions[v2 * 3 + 0]) * 0.5f;
		float y = (mPositions[v1 * 3 + 1] + mPositions[v2 * 3 + 1]) * 0.5f;
		float z = (mPositions[v1 * 3 + 2] + mPositions[v2 * 3 + 2]) * 0.5f;

		float length = std::sqrt(x * x + y * y + z * z) / mRadius;

		x /= length;
		y /= length;
		z /= length;

		uint32_t index = mPositions.size() / 3;
		mPositions.insert(mPositions.end(), { x, y, z });

		return index;
	}

	void saveObj(const std::string& filename) {
		std::ofstream file(filename);

		if (!file.is_open()) {
			std::cerr << "Failed to open file: " << filename << std::endl;
			return;
		}

		for (size_t i = 0; i < mPositions.size(); i += 3) {
			file << "v " << mPositions[i] << " " << mPositions[i + 1] << " " << mPositions[i + 2] << "\n";
		}

		for (size_t i = 0; i < mIndices.size(); i += 3) {
			file << "f "
				<< mIndices[i] + 1 << " "
				<< mIndices[i + 1] + 1 << " "
				<< mIndices[i + 2] + 1 << "\n";
		}

		file.close();
		std::cout << "Mesh saved to " << filename << std::endl;
	}
private:
	float mRadius;
	glm::vec3 mCenter;
	int mResolution;

	std::vector<float> mPositions;
	std::vector<float> mOldPositions;
	std::vector<uint32_t> mIndices;
	std::vector<float> mColors;
	std::vector<float> mVelocity;
	std::vector <float> mMass_inv;
	float mCenterMass_inv;
	glm::vec3 mVelCenter;
	glm::vec3 mOldCenter;
	glm::vec3 mAccelerate;
	std::vector<Constrain> mConstrains;

	GLuint mVAO, mPosVbo, mColVbo, mEBO;
	GLuint mShaderProgram;
	glm::mat4 mModelMatrix;

	Camera* mCamera;
	GLuint mProjectionViewLoc;
	GLuint mModelLoc;
	GLuint mColorLoc;

};


Ball ball(10.0f, glm::vec3(0.0f, 15.0f, 0.0f), 2);
Camera camera(glm::vec3(0.0f, 0.0f, 100.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));




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

	ball.renderInit(&camera);
}


void display() {

	ball.display();

}

void simulation() {
	ball.update(dt / ITER, ITER);
}

int main() {

	init();

	while (!glfwWindowShouldClose(window)) {
		glEnable(GL_CULL_FACE);
		glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		for (int i = 0; i < ITER; i++) {
			simulation();
		}
		

		display();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	return 0;


}