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

using namespace std;
GLFWwindow* window;
uint32_t WIDTH = 800;
uint32_t HEIGHT = 800;
float DT = 0.033f;
uint32_t ITERATION = 10;


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
		mProjectionMatrix = glm::perspective<float>(glm::radians(45.0f), WIDTH / HEIGHT, 0.1f, 100.0f);
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


class Particle {
public:
	Particle(float inv_mass, glm::vec3& position, glm::vec3& acceleration) :
		mInvMass(inv_mass),mPosition(position), mAcceleration(acceleration) {
		mVelocity = glm::vec3(0.0f, 0.0f, 0.0f);
		mOldPosition = position;
	}
	Particle() {}
	~Particle() {}

	glm::vec3 getPosition() {
		return mPosition;
	}

	void addPosition(glm::vec3& translation) {
		if (mInvMass > 1e-5f)
			mPosition += translation;
	}

	void update(float dt) {
		if (mInvMass <= 1e-5f) return;
		mVelocity =( mPosition - mOldPosition) * dt;
		mOldPosition = mPosition;
		mPosition +=  mVelocity * dt + mAcceleration * dt * dt; //glm::vec3(0.0f, -0.01f, 0.0f);//
	}

	float getInvMass() { return mInvMass; }

	glm::vec3& getVelocity() { return mVelocity; }
private:
	float mInvMass;
	glm::vec3 mPosition;
	glm::vec3 mOldPosition;
	glm::vec3 mAcceleration;
	glm::vec3 mVelocity;

};

class Ball {

public:
	Ball(glm::vec3& position, glm::vec3& velocity, float radius) : // 球在往内部走
		mPosition(position),
		mVelocity(velocity),
		mRadius(radius) {}

	Ball() { }
	~Ball() { }

	glm::vec3& getPosition() { return mPosition; }
	float      getRadius() { return mRadius; }

	void update(float dt) {
		mPosition += mVelocity * dt;

	}

private:
	glm::vec3     mVelocity;
	glm::vec3 mPosition;
	float     mRadius;
};


class Constrain {
public:
	Constrain(Particle* p0, Particle* p1) :
		mRestLength(0.0f), // 静态长度
		mParticle1(p0), // 粒子1
		mParticle2(p1), // 粒子1
		mStiffness(0.1f), // 刚度
		mCompliance(0.0f), // 柔度
		mLambda(0.0f)  // lambda
	{
		glm::vec3 p0_to_p1 = mParticle2->getPosition() - mParticle1->getPosition();
		mRestLength = glm::length(p0_to_p1); // 更新静态长度
	}

	void lambdaInit() {
		mLambda = 0.0f;
	}

	void solve(float dt) {
		GLfloat   inv_mass1 = mParticle1->getInvMass();
		GLfloat   inv_mass2 = mParticle2->getInvMass();
		GLfloat   sum_mass = inv_mass1 + inv_mass2;
		if (sum_mass == 0.0f) { return; }
		glm::vec3 p1_minus_p2 = mParticle1->getPosition() - mParticle2->getPosition();
		GLfloat   distance = glm::length(p1_minus_p2);
		GLfloat   constraint = distance - mRestLength; // Cj(x)
		glm::vec3 correction_vector;
	
		mCompliance = MODE_COMPLIANCE[1]; // 查看布料的柔度系数
		mCompliance /= dt * dt;    // a~
		GLfloat dlambda = (-constraint - mCompliance * mLambda) / (sum_mass + mCompliance); // eq.18
		correction_vector = dlambda * p1_minus_p2 / (distance + FLT_EPSILON);                    // eq.17
		mLambda += dlambda;

		mParticle1->addPosition(+inv_mass1 * correction_vector);
		mParticle2->addPosition(-inv_mass2 * correction_vector);

	}


private:
	float    mRestLength;
	Particle* mParticle1;
	Particle* mParticle2;
	float    mStiffness;   
	float    mCompliance;  
	float    mLambda;      

};

class Cloth {
public:
	Cloth(float width, float height, uint32_t nWidth, uint32_t nHeight) {
		mWidth = nWidth;
		mHeight = nHeight;
		mParticles.resize(nWidth * nHeight);

		// push particle
		for (int w = 0; w < mWidth; w++) {
			for (int h = 0; h < mHeight; h++) {
				glm::vec3 pos(width * ((float)w / (float)mWidth) - width * 0.5f,
					-height * ((float)h / (float)mHeight) +  height * 0.5f ,0.0f );
				glm::vec3 gravity(0.0f,-9.8f, 0.0f);
				GLfloat inv_mass = 0.1f;
				if ((h == 0) && (w == 0) ||
					(h == 0) && (w == mWidth - 1)) {
					inv_mass = 0.0f; //fix only edge point
				}
				mParticles[h * mWidth + w] = Particle(inv_mass, pos, gravity);
			}
		}

		// 构建indices
		mIndices.resize((nWidth - 1) * (nHeight - 1) * (nWidth - 1) * (nHeight - 1) * 2);

		for (int h = 0; h < mHeight; h++) {
			for (int w = 0; w < mWidth; w++) {
				if (w == mWidth - 1 && h != mHeight - 1) {
					mIndices[(h * mWidth + w) * 4 + 0] = h * mWidth + w;
					mIndices[(h * mWidth + w) * 4 + 1] = (h + 1) * mWidth + w;
				}
				else if (h == mHeight - 1 && w != mWidth - 1) {
					mIndices[(h * mWidth + w) * 4 + 0] = h * mWidth + w;
					mIndices[(h * mWidth + w) * 4 + 1] = h * mWidth + w + 1;
				}
				else if (h == mHeight - 1 && w == mWidth - 1) {
					
				}
				else {
					mIndices[(h * mWidth + w) * 4 + 0] = h * mWidth + w;
					mIndices[(h * mWidth + w) * 4 + 1] = h * mWidth + w + 1;
					mIndices[(h * mWidth + w) * 4 + 2] = h * mWidth + w;
					mIndices[(h * mWidth + w) * 4 + 3] = (h + 1) * mWidth + w;
				}
 			}
		}

		// 构建约束
		for (int w = 0; w < mWidth; w++) {
			for (int h = 0; h < mHeight; h++) {           // structual constraint
				if (w < mWidth - 1) { 
					mConstrains.push_back(Constrain(&mParticles[h * mWidth + w], &mParticles[h * mWidth + w + 1])); 
				}
				if (h < mWidth - 1) { 
					mConstrains.push_back(Constrain(&mParticles[h * mWidth + w], &mParticles[(h + 1) * mWidth + w])); 
				}
				if (w < mWidth - 1 && h < mHeight - 1) { // shear constraint
					mConstrains.push_back(Constrain(&mParticles[h * mWidth + w], &mParticles[(h + 1) * mWidth + w + 1]));
					mConstrains.push_back(Constrain(&mParticles[h * mWidth + w + 1], &mParticles[(h + 1) * mWidth + w ]));
				}
			}
		}

		for (int w = 0; w < mWidth; w++) {
			for (int h = 0; h < mHeight; h++) {           // bend constraint
				if (w < mWidth - 2) {
					mConstrains.push_back(Constrain(&mParticles[h * mWidth + w], &mParticles[h * mWidth + w + 2]));
				}
				if (h < mHeight - 2) { 
					mConstrains.push_back(Constrain(&mParticles[h * mWidth + w], &mParticles[(h + 2)* mWidth + w ]));
				}
				if (w < mWidth - 2 && h < mHeight - 2) {
					mConstrains.push_back(Constrain(&mParticles[h * mWidth + w], &mParticles[(h + 2) * mWidth + w + 2]));
					mConstrains.push_back(Constrain(&mParticles[h * mWidth + w + 2], &mParticles[(h + 2) * mWidth + w]));
				}
			}
		}
		

	}

	void display() {

		glLineWidth(2.0f);



		glUseProgram(mShaderProgram);


		// projectionViewMatrix
		glUniformMatrix4fv(mProjectionViewLoc, 1, GL_FALSE, glm::value_ptr(mCamera->getProjectionViewMatrix()));

		// modelMatrix
		glUniformMatrix4fv(mModelLoc, 1, GL_FALSE, glm::value_ptr(mModelMatrix));

		glBindVertexArray(mVAO);

		for (int i = 0; i < mParticles.size(); i++) {
			glm::vec3 pos = mParticles[i].getPosition();
			mPositions[i * 3 + 0] = pos.x;
			mPositions[i * 3 + 1] = pos.y;
			mPositions[i * 3 + 2] = pos.z;
		}

		float maxV = 1.5e-3f;

		for (int i = 0; i < mParticles.size(); i++) {

			float vec = glm::length(mParticles[i].getVelocity()) ;
			mColors[i * 3 + 0] = vec / maxV;
			mColors[i * 3 + 1] = 0.0f;
			mColors[i * 3 + 2] = 0.0f;
		}

		glBindBuffer(GL_ARRAY_BUFFER, mPosVbo);
		glBufferSubData(GL_ARRAY_BUFFER, 0, mPositions.size() * sizeof(float), mPositions.data());

		glBindBuffer(GL_ARRAY_BUFFER, mColVbo);
		glBufferSubData(GL_ARRAY_BUFFER, 0, mColors.size() * sizeof(float), mColors.data());

		glDrawElements(GL_LINES, mIndices.size(), GL_UNSIGNED_INT, 0);

		glBindVertexArray(0);

		glUseProgram(0);
	
	}

	void renderInit(Camera *camera) {

		mCamera = camera;

		mModelMatrix = glm::mat4(1.0f);
		mColors.resize(mWidth * mHeight * 3);
		mPositions.resize(mWidth * mHeight * 3);
		for (int i = 0; i < mWidth; i++) {
			for (int j = 0; j < mHeight; j++) {
				uint32_t base = i * mHeight + j;
				mPositions[base * 3 + 0] = mParticles[base].getPosition().x;
				mPositions[base * 3 + 1] = mParticles[base].getPosition().y;
				mPositions[base * 3 + 2] = mParticles[base].getPosition().z;
			}
		}

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
		

		// 创建 VBO（顶点缓冲对象）
		glGenBuffers(1, &mColVbo);
		glBindBuffer(GL_ARRAY_BUFFER, mColVbo);
		glBufferData(GL_ARRAY_BUFFER, mColors.size() * sizeof(float), mColors.data(), GL_STATIC_DRAW);


		// 配置顶点属性（位置属性）
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
		glEnableVertexAttribArray(1);


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



	
	}
	

	void update(float dt, int iteration) {
		
		for (auto ball = mBalls.begin(); ball != mBalls.end(); ball++) {
			if ((*ball).getPosition().z < -2.0f) {
				continue;
			}
			(*ball).update(dt);
			
		}

		for (auto particle = mParticles.begin(); particle != mParticles.end(); particle++) {
			(*particle).update(dt); // predict position
		}
		
		for (auto constraint = mConstrains.begin(); constraint != mConstrains.end(); constraint++) {
			(*constraint).lambdaInit();
		}
		
		for (int i = 0; i < iteration; i++) {
			for (auto particle = mParticles.begin(); particle != mParticles.end(); particle++) {
				for (auto ball = mBalls.begin(); ball != mBalls.end(); ball++) {
					if ((*ball).getPosition().z < -2.0f) {
						continue;
					}
					glm::vec3 vec = (*particle).getPosition() - ball->getPosition();
					float     length = glm::length(vec); // 检查质点到球的距离
					float     radius = ball->getRadius() * 1.2f; // fake radius
					if (length < radius) {
						(*particle).addPosition(glm::normalize(vec) * (radius - length)); // 简单的把位置修正到球外面
					}
				}	
			}

			for (auto constraint = mConstrains.begin(); constraint != mConstrains.end(); constraint++) {
				(*constraint).solve(dt);
			}
		
		}
		

	}
	
	void launchBall(glm::vec3 ballPosition) {
		Ball ball(ballPosition, glm::vec3(0.0f, 0.0f, -1.0f), 1.0f);
		mBalls.push_back(ball);
	}


private:
	uint32_t mWidth;
	uint32_t mHeight;
	std::vector<Particle> mParticles;
	std::vector<Constrain> mConstrains;
	std::vector<Ball> mBalls;


	GLuint mVAO, mPosVbo, mColVbo, mEBO, mShaderProgram;
	std::vector<uint32_t> mIndices;
	std::vector<float> mPositions;
	std::vector<float> mColors;
	glm::mat4 mModelMatrix;

	Camera* mCamera;
	GLuint mProjectionViewLoc;
	GLuint mModelLoc;

};



Cloth cloth(3.0f, 3.0f, 30, 30);
Camera camera(glm::vec3(0.0f,0.0f,5.0f), glm::vec3(0.0f, 0.0f, 0.0f),  glm::vec3(0.0f, 1.0f, 0.0f));


void simulation() {
	cloth.update(DT, ITERATION);
}

void display() {

	cloth.display();

}

void init() {

	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	//1.2 设置OpenGL启用核心模式（非立即渲染模式）
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	window = glfwCreateWindow(800, 800, "xpbd", NULL, NULL);
	if (window == NULL) {

	}
	//设置当前窗体对象为OpenGL的绘制舞台
	glfwMakeContextCurrent(window);

	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cout << "Failed to initialize GLAD" << std::endl;

	}

	cloth.renderInit(&camera);
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		int width, height;
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos); // 获取鼠标点击的屏幕坐标
		glfwGetFramebufferSize(window, &width, &height);

		// 将鼠标坐标转换到 NDC（范围为 -1 到 1）
		float x = (2.0f * xpos) / width - 1.0f;
		float y = 1.0f - (2.0f * ypos) / height; // Y 轴需要反转
		float z = 1.0f; // 假设鼠标点击在近平面

		// 构造 NDC 坐标
		glm::vec4 ndcPos = glm::vec4(x, y, z, 1.0f);

		// 逆投影矩阵（Projection^-1 * View^-1）
		glm::mat4 invProj = glm::inverse(camera.getProjectionMatrix());
		glm::mat4 invView = glm::inverse(camera.getViewMatrix());

		// 将 NDC 转换为视图坐标
		glm::vec4 viewPos = invProj * ndcPos;
		viewPos = glm::vec4(viewPos.x, viewPos.y, -1.0f, 0.0f); // 设置方向向量

		// 将视图坐标转换为世界坐标
		glm::vec4 worldPos = invView * viewPos;

		// 归一化方向
		glm::vec3 worldDir = glm::normalize(glm::vec3(worldPos));


		cloth.launchBall(glm::vec3(worldDir.x, worldDir.y, 0.5f));
	}
}



int main() {

	init();

	glfwSetMouseButtonCallback(window, mouseButtonCallback);

	while (!glfwWindowShouldClose(window)) {
		glDisable(GL_CULL_FACE);
		glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);


		simulation();
		

		display();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	return 0;


}