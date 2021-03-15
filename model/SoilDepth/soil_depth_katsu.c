/*--
soil_depth_katsu.c
Soil depth pattern prediction using Dietrich et al. 1995 model
by Yasushi TANAKA
This program is developed on the 
[Devis2dim.c]
[diffusion.c]
[soil_production_1995]
Geomorphic Development Simulation Program  (Devis Model, 2Dimension, 2001/05/15) program.
Note that the NV treatment is difficult, especially when it is in the certificate data
2007/5/21, -> 2018/3/10修正，2021/3/13修正
--*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
float calc_cvt(int,int,float);
float spf(float);

//----------------- input parameter ----------------------
#define IFN "katsu01m_dem.flt"                  // Input DEM data (float matrix)
#define OFN "katsu01m_SoilDepth.flt"            // Output Soil Depth (float matrix, h2)
#define OFN2 "katsu01m_ElevationChange.flt"     // Output Elevation Change (float matrix, dem2-dem0)
#define WE 1536                                 // Number of pixel (Column)
#define NS 2358                                 // Number of pixel (row)
#define GD 1.0                                  // Grid Distance (meter)
#define BL 0                                    // Base Lebel Height (m)
#define K 0.005                     // Diffusivity (m2/year), see Nogami(1980), Reneau(1988)
#define M_TS 10                     // Meaning of 1 Step (year: 10 year)
#define TS 100                      // Time Step: number of iterataion (100 times)
#define CI 10                       // value check interval
#define cWE 1010                    // value check at WE [site2 in Kawano(2020)]
#define cNS 921                     // value check at NS [site2 in Kawano(2020)]
#define NV -9999                    // No data value
#define ISD 0.3                     // Initial Soil Depth (m)
#define Rho_RS 1.7                  // density ratio of rock to soil
//----------------- input parameter end -------------------

float dem0[NS][WE];      // [Row][Col]: initial DEM(z)
float dem1[NS][WE];      // [Row][Col]: DEM(z) before changed
float dem2[NS][WE];      // [Row][Col]: DEM(z) after changed
float h0[NS][WE];        // [Row][Col]: initial soil thickness (h=z-e)
float h1[NS][WE];        // [Row][Col]: soil thickness before changed
float h2[NS][WE];        // [Row][Col]: soil thickness after changed
float e0[NS][WE];        // [Row][Col]: initial bedrock height (e=z-h)
float e1[NS][WE];        // [Row][Col]: bedrock height before changed
float e2[NS][WE];        // [Row][Col]: bedrock height after changed

int i,j;

int main(){

	double t0, t2;
	t0 = clock();
	int t,st;
	FILE *infile,*outfile,*outfile2;
	float de_dt=0;
	float gd,cvt,diff;

	infile = fopen(IFN,"rb");
	outfile = fopen(OFN,"wb");
	outfile2 = fopen(OFN2,"wb");
	
	gd=GD*1.0;
	printf("x & y Interval = %f m\n",gd);
	
	// read DEM
	fread(&dem0[0][0],WE*NS*sizeof(dem0[0][0]),1,infile);        // initial topography
	memcpy(dem1,dem0,sizeof(dem0));                              // copy from dem0[][] to dem1[][]
	printf("t=0, dem1[%d][%d]=%f\n",cNS,cWE,dem1[cNS][cWE]);
	// set initial soil depth and bedrock value;
	for(i=0;i<NS;i++){
		for(j=0;j<WE;j++){
			if(dem0[i][j]==NV){
				h0[i][j]=NV;
				e0[i][j]=NV;
			}
			else{
				h0[i][j]=ISD;
				e0[i][j]=dem0[i][j]-ISD;
			}
		}
	}
	memcpy(h1,h0,sizeof(h0));                              // copy from h0[][] to h1[][]
	memcpy(e1,e0,sizeof(e0));                              // copy from e0[][] to e[][]
	printf("t=0, h1[%d][%d]=%f\n",cNS,cWE,h1[cNS][cWE]);
	printf("t=0, e1[%d][%d]=%f\n",cNS,cWE,e1[cNS][cWE]);
	//printf("End reading file\n");
	
	printf("Start Calculation\n");
	
	for(t=1;t<=TS;t++){
		
		if(t!=1){
			memcpy(dem1,dem2,sizeof(dem2));               // copy from dem1[][] to dem2[][]
			memcpy(h1,h2,sizeof(h2));                     // copy from h1[][] to h2[][]
			memcpy(e1,e2,sizeof(e2));                     // copy from e1[][] to e2[][]
		}
		
		for(i=0;i<NS;i++){
			for(j=0;j<WE;j++){
				// frame is NV
				if(i==0 || i==NS-1 || j==0 || j==WE-1 ){
					dem2[i][j]=NV; h2[i][j]=NV; e2[i][j]=NV;
				}
				else if(dem1[i][j]==NV){
					dem2[i][j]=NV; h2[i][j]=NV; e2[i][j]=NV;
				}
				else if(dem1[i-1][j-1]==NV || dem1[i-1][j  ]==NV || dem1[i-1][j+1]==NV || 
			          dem1[i  ][j-1]==NV                       || dem1[i  ][j+1]==NV || 
			          dem1[i+1][j-1]==NV || dem1[i+1][j  ]==NV || dem1[i+1][j+1]==NV){
					dem2[i][j]=dem1[i][j]; h2[i][j]=h1[i][j]; e2[i][j]=e1[i][j];
				}
				else{
					// Soil production term
					if(h1[i][j]<=0){
						de_dt=0;
						e2[i][j]=e1[i][j];
						h2[i][j]=0;
						dem2[i][j]=e2[i][j];
					}
					else{
						//de_dt=spf(h1[i][j]);            // de_dt(m/y)
						de_dt=spf(h1[i][j])*M_TS;         // de_dt
						e2[i][j]=e1[i][j]-de_dt;          // elevation of bedrock(m)
						h2[i][j]=h1[i][j]+Rho_RS*de_dt;   // soil thickness
						dem2[i][j]=e2[i][j]+h2[i][j];     // elevation(z)
						
						// Transportation term
						cvt = calc_cvt(i,j,gd);           // calculate Topographic Curvature
						h2[i][j]=h2[i][j]+(cvt*K*M_TS);   // topographic change (soil movement)
						                                  // もしh2がここで負の値になっても，次のステップの前のif分でe2の値に戻る
						dem2[i][j]=e2[i][j]+h2[i][j];     // caculate elevation (again)　soilのみ流す
					}
					/*if(dem2[i][j]<BL){
						dem2[i][j]=BL;
					}*/
					
				}
				// check result
				if(t%CI==0 && i==cNS && j==cWE){
					//printf("curvature[%d][%d]=%f\n",i,j,cvt);
					//printf("K=%f(m2/year): Diffusivity(constant)\n",K);
					//printf("K*∇2z[%d][%d]=%f\n",i,j,K*cvt);
					//printf("dz/dt=%f(m/%dyears)\n",TS*K*cvt,TS);
					//printf("t=%d, dem2[%d][%d]=%f\n",t,i,j,dem2[i][j]);
					//printf("Soil production rate at %.1f m=%f (m/y)\n",h1[i][j],de_dt);
					//printf("Soil production rate at %.3f m=%f (m/100y)\n",h1[i][j],de_dt);
					//printf("t=%d, e2[%d][%d]=%f\n",t,i,j,e2[i][j]);
					printf("%d [year], SoilDepth@[%d][%d]=%.4f [m]\n",t*M_TS,i,j,h2[i][j]);
					//printf("t=%d, dem2[%d][%d]=%f\n",t,i,j,dem2[i][j]);
				}
			}
		}
	}
	printf("End Calculation\n");
	
	// Write changed landform DEM
	//fwrite(&dem2[0][0],WE*NS*sizeof(dem2[0][0]),1,outfile);
	fwrite(&h2[0][0],WE*NS*sizeof(h2[0][0]),1,outfile);
	printf("Wrote Last Soil Depth Map: %s\n",OFN);
	
	// Write difference DEM
	for(i=0;i<NS;i++){
		for(j=0;j<WE;j++){
			if(dem0[i][j]==NV) diff=NV;
			else{
				diff=dem2[i][j]-dem0[i][j];
			}
			fwrite(&diff,sizeof(diff),1,outfile2);
		}
	}
	//printf("Wrote difference DEM\n");
	
	fclose(infile);
	fclose(outfile);
	fclose(outfile2);

	t2 = clock();
	printf("Calculate time; %.2f[s]\n", (t2 - t0) / CLOCKS_PER_SEC);
}


// function of topographic curvature calculation
float calc_cvt(int i, int j, float gd){
	float cvt;
	float e2468,e1379,e5;
	
	e2468=dem1[i-1][j  ]+dem1[i  ][j-1]+dem1[i  ][j+1]+dem1[i+1][j  ];
	e1379=dem1[i-1][j-1]+dem1[i-1][j+1]+dem1[i+1][j-1]+dem1[i+1][j+1];
	e5=dem1[i][j];
	cvt=(2*e2468+e1379-12*e5)/(4*(GD*GD));

	return cvt;
}

// function of soil production rate
// input (m), output(m)
float spf(float h){
	float h_cm=h*100;
	float de_dt_cm;
	float de_dt_m;
	float P0=0.019;      // Production rate (cm/yr) at zero soil thickness (Constant) 
	float m=0.05;        // rate constant (1/cm) (Constant)
	
	de_dt_cm=P0*exp(-m*h_cm);
	de_dt_m=de_dt_cm/100;
	
	return de_dt_m;
}

