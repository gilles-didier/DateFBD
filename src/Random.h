#ifndef RandomF
#define RandomF


#ifdef __cplusplus
extern "C" {
#endif
void initRandom();
/*return an uniform  random value in [0,1]*/
double getUniformStd();
/*return an uniform  random value in [0, m]*/
double getUniformCont(double m);
/*return an uniform discrete random value in {0, 1, m}*/
double getUniformDisc(int n);
/*return an exponential distributed random value in {0, 1, m}*/
double getExponential(double l);

#ifdef __cplusplus
}
#endif

#endif
