#version 330 core
layout(location = 0) in vec3 aPos;


uniform mat4 projectionViewMatrix;
uniform mat4 modelMatrix;

void main(){
	vec4 tranPos = vec4(aPos, 1.0);

	gl_Position  = projectionViewMatrix * modelMatrix * tranPos;
	
}