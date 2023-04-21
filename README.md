# COMP559FinalProject

## Installation
### 1. Requirements
I ran this project on Windows 10 - Visual Studio 2022

Nearly identical to our other assignments, particularly assignment 2.
Requires:
* glew 2.1.0
* glfw 3.3.8
* freetype 2.12.1
* glm 0.9.9

**You must set an environment variable named GLM_INCLUDE_DIR to point to the installation of glm
You must set an environment variable named THIRDPARTY_DIR to point to the other installations (freetype, glfw, glew)**

### 2. CMAKE
Use CMAKE: create a new folder called build within the COMP559FinalProject directory and tell CMAKE to build there. Then configure and generate

### 3. Debugging
If any errors occur try deleting the CMAKE cache, restarting CMAKE and deleting the build folder. Other errors can show up from 
