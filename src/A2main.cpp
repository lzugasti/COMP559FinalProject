/**
 * Provided code for particle system simulator.
 * This code provides the mouse interface for clicking and dragging particles, and the
 * code to draw the system.  When the simulator is running system.advanceTime is called
 * to numerically integrate the system forward.
 * @author kry
 */
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "GLSL.h"
#include "MatrixStack.h"
#include "Program.h"

#include "Text.hpp"

#include "VoxelSystem.hpp"
#include "TestSystems.hpp"
#include "SanityCheck.hpp"
#include "RecordVideo.hpp"
#include "Voxel.hpp"

using namespace std;

VoxelSystem voxelSystem;
SanityCheck sanityCheck;
TestSystems* testSystems;
bool run = false;
bool sanity = true;
float stepsize = 1.0f / 60.0f; // 60 Hz screen refresh
int substeps = 1;

// parameters for interacting with particles
float maxDist = 150;
float minDist = 50;
float grabThresh = 10;

// for grabbing particles
Particle* p1 = NULL;
Particle* p2 = NULL;
float d1 = 0;
float d2 = 0;

int xdown = 0;
int ydown = 0;
int xcurrent = 0;
int ycurrent = 0;
bool mouseDown = false;
bool mouseInWindow = false;
bool grabbed = false;
bool wasPinned = false;
bool closeToParticlePairLine = false;
bool canEdit = true; // when simulation hasn't been run/stepped yet

string scene = "";
int testSceneID = 0;
int showKeyBindings = 0;
// for openGL
GLFWwindow* window; // Main application window
string RES_DIR = ""; // Where data files live
shared_ptr<Program> progIM; // immediate mode

ImageRecorder imageRecorder;
bool recordFrames = false; // toggled by keyboard
bool recordFrame = false; // set to record when stepped in display


static void error_callback(int error, const char* description) {
	cerr << description << endl;
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if ((action != GLFW_PRESS) && (action != GLFW_REPEAT)) return;
	if (key == GLFW_KEY_ESCAPE) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	} 
	else if (key == GLFW_KEY_G)
	{
		voxelSystem.useGravity = !voxelSystem.useGravity;
	}
	else if (key == GLFW_KEY_C)
	{
		voxelSystem.useCollisions = !voxelSystem.useCollisions;
	}
	else if (key == GLFW_KEY_SPACE) {
		run = !run;
	}
	else if (key == GLFW_KEY_R) {
		cout << "Resetting" << endl;
		voxelSystem.reset();
	}
	else if (key == GLFW_KEY_P) {
		if (mods & GLFW_MOD_SHIFT) {
			voxelSystem.peak += 500;
		}
		else
		{
			voxelSystem.peak -= 500;
			if (voxelSystem.peak < 500)
			{
				voxelSystem.peak = 500;
			}
		}
	}
	else if (key == GLFW_KEY_T) {
		if (mods & GLFW_MOD_SHIFT) {
			voxelSystem.T += 0.1;
		}
		else
		{
			voxelSystem.T -= 0.1;
			if (voxelSystem.T < 0.1)
			{
				voxelSystem.T = 0.1;
			}
		}
	}
	else if (key == GLFW_KEY_B) {
		if (mods & GLFW_MOD_SHIFT) {
			voxelSystem.b += 1;
		}
		else
		{
			voxelSystem.b -= 1;
			if (voxelSystem.b < 1)
			{
				voxelSystem.b = 1;
			}
		}
	}
	else if (key == GLFW_KEY_S) {
		if (mods & GLFW_MOD_SHIFT) {
			voxelSystem.speed += 50;
		}
		else
		{
			voxelSystem.speed -= 50;
			if (voxelSystem.speed < 50)
			{
				voxelSystem.speed = 50;
			}
		}
	}
	else if (key == GLFW_KEY_D) {
		voxelSystem.decay = (voxelSystem.decay + 1) % 4;
	}
	else if (key == GLFW_KEY_L) {
		if (mods & GLFW_MOD_SHIFT) {
			voxelSystem.threshold += 500;
		}
		else
		{
			voxelSystem.threshold -= 500;
			if (voxelSystem.threshold < 0)
			{
				voxelSystem.threshold = 0;
			}
		}
	}
	else if (key == GLFW_KEY_W) {
		if (mods & GLFW_MOD_SHIFT) {
			voxelSystem.rectangleWidth += 1;
		}
		else
		{
			voxelSystem.rectangleWidth -= 1;
			if (voxelSystem.rectangleWidth < 1)
			{
				voxelSystem.rectangleWidth = 1;
			}
		}
	}
	else if (key == GLFW_KEY_H) {
		if (mods & GLFW_MOD_SHIFT) {
			voxelSystem.rectangleHeight += 1;
		}
		else
		{
			voxelSystem.rectangleHeight -= 1;
			if (voxelSystem.rectangleHeight < 1)
			{
				voxelSystem.rectangleHeight = 1;
			}
		}
	}
	else if (key == GLFW_KEY_O) {
		voxelSystem.spawnType = (voxelSystem.spawnType + 1) % 2;
	}
	else if (key == GLFW_KEY_I) {
		voxelSystem.displayLinks = !voxelSystem.displayLinks;
	}
}

void mouse_pos_callback(GLFWwindow* window, double x, double y) {
	xcurrent = floor(x);
	ycurrent = floor(y);
	if (mouseDown) { // dragged
	}
}

void mouse_callback(GLFWwindow* window, int button, int action, int mods) {

	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		double x;
		double y;
		glfwGetCursorPos(window, &x, &y);
		xcurrent = floor(x);
		ycurrent = floor(y);

		bool pressed = false;
		bool released = false;
		if (button == GLFW_MOUSE_BUTTON_LEFT) {
			if (GLFW_PRESS == action) {
				voxelSystem.createShape(xcurrent, ycurrent, /*rand() % 20 - 10*/0.0, 0.0, /*(rand() % 100) / 10.0*/0.0, /**/0.0, (rand() % 100) / 100.0, (rand() % 100) / 100.0, (rand() % 100) / 100.0, 5, rand() % 4 + 2);
				pressed = !mouseDown;
				mouseDown = true;
			} else if (GLFW_RELEASE == action) {
				released = mouseDown;
				mouseDown = false;
			}
		}
	}
	if (button == GLFW_MOUSE_BUTTON_RIGHT) {
		if (GLFW_PRESS == action) {
			voxelSystem.createExplosion(xcurrent, ycurrent);
		}
	}
}

static void init() {
	GLSL::checkVersion();

	// Check how many texture units are supported in the vertex shader
	int tmp;
	glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &tmp);
	cout << "GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS = " << tmp << endl;
	// Check how many uniforms are supported in the vertex shader
	glGetIntegerv(GL_MAX_VERTEX_UNIFORM_COMPONENTS, &tmp);
	cout << "GL_MAX_VERTEX_UNIFORM_COMPONENTS = " << tmp << endl;
	glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &tmp);
	cout << "GL_MAX_VERTEX_ATTRIBS = " << tmp << endl;

	// Set background color.
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	// Enable z-buffer test.
	glEnable(GL_DEPTH_TEST);

	// Initialize the GLSL programs.
	progIM = make_shared<Program>();
	progIM->setVerbose(true);
	progIM->setShaderNames(RES_DIR + "simple_vert.glsl", RES_DIR + "simple_frag.glsl");
	progIM->init();
	progIM->addUniform("P");
	progIM->addUniform("MV");
	progIM->setVerbose(false);

	initTextRender(RES_DIR);


	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glDisable(GL_DEPTH_TEST);
	voxelSystem.init();
	
	
	GLSL::checkError(GET_FILE_LINE);
}

/** draws a line from the given point to the given particle */
void drawLineToParticle(double x, double y, Particle* p, double d) {
	if (p == NULL) return;
	if (d > maxDist) return;
	double col = d < minDist ? 1 : (maxDist - d) / (maxDist - minDist);
	glColor4d(1 - col, 0, col, 0.75f);
	glBegin(GL_LINES);
	glVertex2d(x, y);
	glVertex2d(p->p.x, p->p.y);
	glEnd();
}

void display() {
	// set up projection for drawing in pixel units...

	if (run) {
		for (int i = 0; i < substeps; i++) {
			voxelSystem.advanceTime(stepsize / substeps);
		}
		if ( recordFrames ) recordFrame = true;
	}
	voxelSystem.display();

	
}

static void render() {
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	float aspect = width / (float)height;
	glViewport(0, 0, width, height);
	voxelSystem.width = width;
	voxelSystem.height = height;

	
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Create matrix stacks.
	auto P = make_shared<MatrixStack>();
	auto MV = make_shared<MatrixStack>();

	// Setup projection for drawing in pixel coordinates
	float scale = 1;
	glm::mat4 projection = glm::ortho(0.0f, static_cast<float>(width * scale), static_cast<float>(height * scale), 0.0f);
	glm::mat4 modelview = glm::identity<glm::mat4>();

	progIM->bind();
	glUniformMatrix4fv(progIM->getUniform("P"), 1, GL_FALSE, &projection[0][0]);
	glUniformMatrix4fv(progIM->getUniform("MV"), 1, GL_FALSE, &modelview[0][0]);
	display();
	progIM->unbind();

	stringstream ss;
	ss << "\n";
		ss << "---------- Controls ----------" << "\n";
		ss << "\n";
		ss << "Right Click: Explosion!!!\n";
		ss << "Left Click: Spawn Objects (\"o\" to change type)\n";
		ss << "\n";
		ss << "---------- Simulation Settings ----------" << "\n";
		ss << "\n";
		ss << "[g] use gravity = " << voxelSystem.useGravity << "\n";
		ss << "[c] use collisions = " << voxelSystem.useCollisions << "\n";
		ss << "[space] start/pause = " << run << "\n";
		ss << "[r] restart simulation\n";
		ss << "[i] display links\n";
		ss << "\n";
		ss << "---------- Explosion Settings ----------" << "\n";
		ss << "\n";
		ss << "[p,P] peak - (strength of shock front) = " << voxelSystem.peak << "\n";
		ss << "[t,T] T - (length of positive pressure phase) = " << voxelSystem.T << "\n";
		ss << "[b,B] b - (larger means shorter negative phase) = " << voxelSystem.b << "\n";
		ss << "[s,S] speed - (How fast shock wave moves) = " << voxelSystem.speed << "\n";
		ss << "[d] distance decay (0 -> 3: no decay -> cubic decay) = " << voxelSystem.decay << "\n";
		ss << "\n";
		ss << "---------- Object Spawning Settings ----------" << "\n";
		ss << "\n";
		ss << "[l,L] break threshold = " << voxelSystem.threshold << "\n";
		ss << "[o] spawn type (0 = random shapes | 1 = rectangles) = " << voxelSystem.spawnType << "\n";
		ss << "[w,W] rectangle width = " << voxelSystem.rectangleWidth << "\n";
		ss << "[h,H] rectangle height = " << voxelSystem.rectangleHeight << "\n";
	string text = ss.str();
	RenderString(projection, modelview, 600, 50, 0.4, text);

	if (recordFrame) {
		recordFrame = false;
		imageRecorder.writeCurrentFrameToFile(window);
	}

	GLSL::checkError(GET_FILE_LINE);
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "Please specify the resource directory." << endl;
		return 0;
	} else if (argc < 3) {
		cout << "Unspecified scene to run, using default value " << scene << endl;
	} else {
		scene = string(argv[2]);
	}
	RES_DIR = argv[1] + string("/");

	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if (!glfwInit()) {
		return -1;
	}
	// https://en.wikipedia.org/wiki/OpenGL
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	// glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
	// glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	// glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(1280, 720, "COMP 559 W23 - PROJECT - LEOPOLDO ZUGASTI", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}
	glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	// Set vsync.
	glfwSwapInterval(1);
	// Set interaction callbacks.
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouse_callback);
	glfwSetCursorPosCallback(window, mouse_pos_callback);

	// Initialize scene.
	init();
	// Loop until the user closes the window.
	while (!glfwWindowShouldClose(window)) {
		// Render scene.
		render();
		// Swap front and back buffers.
		glfwSwapBuffers(window);
		// Poll for and process events.
		glfwPollEvents();
	}
	// Quit program.
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}

