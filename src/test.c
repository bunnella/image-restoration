#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lib/png/png.h"

int main(int argc, char *argv[]) {
	printf("HELLO FROM THE OTHER SIZE.\n");
	if (argc == 1) {
		printf("Give me a file to open, son.\n");
		return 1;
	}

	png_image image;
	memset(&image, 0, sizeof(image));
	image.version = PNG_IMAGE_VERSION;
	image.format = PNG_FORMAT_RGBA;

	if (!png_image_begin_read_from_file(&image, argv[1])) {
		// error handling
		return 1;
	}

	png_bytep buffer = malloc(PNG_IMAGE_SIZE(image));

	if (!png_image_finish_read(&image, NULL, buffer, 0, NULL)) {
		// error handling...
		return 1;
	}

	free(buffer);
	return 0;
}
