#ifndef MyRF
#define MyRF

#ifdef __cplusplus
extern "C" {
#endif

void error(const char * format, ...);
void warning(const char * format, ...);
double lgammafn(double x);
double unif_rand();

#ifdef __cplusplus
}
#endif

#endif
