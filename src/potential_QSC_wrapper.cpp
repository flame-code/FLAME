#include <stdio.h>
#include "potential_QSC_class.h"

QSC test;

extern "C" void qsc_energy_forces_init_(void);
extern "C" void qsc_energy_forces_final_(void);
extern "C" void qsc_energy_forces_(int *nat,double *rat,int *matomicnum,double *fat,
	double *epot,double *cellvec);
//****************************************************************************************
//int main () {
//	int iat;
//	int nat,matomicnum[2];
//	double rat[2][3],fat[2][3],epot,cellvec[3][3];
//	nat=2;
//	rat[0][0]=1.0; rat[0][1]=1.0; rat[0][2]=1.0;
//	rat[1][0]=3.0; rat[1][1]=1.0; rat[1][2]=1.0;
//	matomicnum[0]=79;
//	matomicnum[1]=79;
//	cellvec[0][0]=20.0;cellvec[0][1]= 0.0;cellvec[0][2]= 0.0;
//	cellvec[1][0]= 0.0;cellvec[1][1]=20.0;cellvec[1][2]= 0.0;
//	cellvec[2][0]= 0.0;cellvec[2][1]= 0.0;cellvec[2][2]=20.0;
//	qsc_energy_forces_(&nat,rat[0],matomicnum,fat[0],&epot,cellvec[0]);
//	printf("epot= %.15f \n",epot);
//	for(iat=0;iat<nat;iat++) {
//		printf("%20.15f %20.15f %20.15f \n",fat[iat][0],fat[iat][1],fat[iat][2]);
//	}
//	return 1;
//}
//****************************************************************************************
void qsc_energy_forces_init_(void) {
	test.QSC_constructor();
}
//****************************************************************************************
void qsc_energy_forces_final_(void) {
	test.QSC_destructor();
}
//****************************************************************************************
void qsc_energy_forces_(int *nat,double *rat,int *matomicnum,double *fat,
	double *epot,double *cellvec) {
	//printf("nat= %d\n",*nat);
	//printf("rat= %20.15f %20.15f %20.15f\n",rat[0],rat[1],rat[2]);
	//printf("rat= %20.15f %20.15f %20.15f\n",rat[3],rat[4],rat[5]);
	//printf("cutoff= %20.5f\n",test.cutoff);
	test.force(*nat,rat,matomicnum,fat,epot,cellvec);
	//printf("epot= %20.15f\n",*epot);
}
//****************************************************************************************
