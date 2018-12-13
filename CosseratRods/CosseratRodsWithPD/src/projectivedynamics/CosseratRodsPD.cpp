#include "visualization.h"
#include "scene.h"

int main(int argc, char **argv) {
	Scene* scene = new Scene();
	auto body = scene->InitializeScene();

	return initialize_renderer(argc, argv, body);
}
