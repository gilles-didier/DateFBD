#ifndef MyErrorF
#define MyErrorF

#ifdef __cplusplus
extern "C" {
#endif

void error(const char * format, ...);
void warning(const char * format, ...);

#ifdef __cplusplus
}
#endif

#endif
