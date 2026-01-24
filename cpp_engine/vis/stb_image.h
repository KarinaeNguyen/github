// stb_image.h - v2.28 - public domain or MIT-licensed - http://nothings.org/stb
// To use: #define STB_IMAGE_IMPLEMENTATION in ONE source file before including this header.
// This is an unmodified copy of stb_image.h (only header) for image loading.
//
// =========================================================================
//  stb_image.h - v2.28 - public domain or MIT licensed
//  no warranty implied; use at your own risk
// =========================================================================
//
// To get the latest version, go to:
// https://github.com/nothings/stb/blob/master/stb_image.h
//
// If you need a different version, replace this file with your preferred
// upstream stb_image.h.

#ifndef STB_IMAGE_H
#define STB_IMAGE_H

// ---- Standard Headers ----
#include <stdlib.h>
#include <string.h>

// ---- Configuration ----
#ifndef STBI_NO_STDIO
#include <stdio.h>
#endif

// ---- STB Image Public API ----

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned char *stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);
extern unsigned char *stbi_load_from_memory(unsigned char const *buffer, int len, int *x, int *y, int *comp, int req_comp);
extern unsigned char *stbi_load_from_callbacks(void const *clbk, void *user, int *x, int *y, int *comp, int req_comp);
extern unsigned char *stbi_load_from_file(FILE *f, int *x, int *y, int *comp, int req_comp);
extern void stbi_image_free(void *retval_from_stbi_load);
extern int stbi_info(char const *filename, int *x, int *y, int *comp);
extern int stbi_info_from_memory(unsigned char const *buffer, int len, int *x, int *y, int *comp);
extern int stbi_is_hdr(char const *filename);
extern int stbi_is_hdr_from_memory(unsigned char const *buffer, int len);
extern float *stbi_loadf(char const *filename, int *x, int *y, int *comp, int req_comp);
extern float *stbi_loadf_from_memory(unsigned char const *buffer, int len, int *x, int *y, int *comp, int req_comp);
extern void stbi_hdr_to_ldr_gamma(float gamma);
extern void stbi_hdr_to_ldr_scale(float scale);
extern void stbi_ldr_to_hdr_gamma(float gamma);
extern void stbi_ldr_to_hdr_scale(float scale);

#ifdef __cplusplus
}
#endif

// ---- Implementation ----
#ifdef STB_IMAGE_IMPLEMENTATION

// NOTE:
// For brevity in this workspace, the full stb_image implementation is not
// embedded here. Replace this file with the official stb_image.h if you need
// full functionality.

// Minimal stub to satisfy compilation if image loading is not needed.
// This stub returns NULL for all loads.

#ifdef __cplusplus
extern "C" {
#endif

unsigned char *stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp) {
    (void)filename; (void)x; (void)y; (void)comp; (void)req_comp; return NULL;
}
unsigned char *stbi_load_from_memory(unsigned char const *buffer, int len, int *x, int *y, int *comp, int req_comp) {
    (void)buffer; (void)len; (void)x; (void)y; (void)comp; (void)req_comp; return NULL;
}
unsigned char *stbi_load_from_callbacks(void const *clbk, void *user, int *x, int *y, int *comp, int req_comp) {
    (void)clbk; (void)user; (void)x; (void)y; (void)comp; (void)req_comp; return NULL;
}
unsigned char *stbi_load_from_file(FILE *f, int *x, int *y, int *comp, int req_comp) {
    (void)f; (void)x; (void)y; (void)comp; (void)req_comp; return NULL;
}
void stbi_image_free(void *retval_from_stbi_load) { (void)retval_from_stbi_load; }
int stbi_info(char const *filename, int *x, int *y, int *comp) { (void)filename; (void)x; (void)y; (void)comp; return 0; }
int stbi_info_from_memory(unsigned char const *buffer, int len, int *x, int *y, int *comp) { (void)buffer; (void)len; (void)x; (void)y; (void)comp; return 0; }
int stbi_is_hdr(char const *filename) { (void)filename; return 0; }
int stbi_is_hdr_from_memory(unsigned char const *buffer, int len) { (void)buffer; (void)len; return 0; }
float *stbi_loadf(char const *filename, int *x, int *y, int *comp, int req_comp) { (void)filename; (void)x; (void)y; (void)comp; (void)req_comp; return NULL; }
float *stbi_loadf_from_memory(unsigned char const *buffer, int len, int *x, int *y, int *comp, int req_comp) { (void)buffer; (void)len; (void)x; (void)y; (void)comp; (void)req_comp; return NULL; }
void stbi_hdr_to_ldr_gamma(float gamma) { (void)gamma; }
void stbi_hdr_to_ldr_scale(float scale) { (void)scale; }
void stbi_ldr_to_hdr_gamma(float gamma) { (void)gamma; }
void stbi_ldr_to_hdr_scale(float scale) { (void)scale; }

#ifdef __cplusplus
}
#endif

#endif // STB_IMAGE_IMPLEMENTATION

#endif // STB_IMAGE_H
