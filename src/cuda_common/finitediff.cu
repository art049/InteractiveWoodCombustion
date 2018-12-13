/* finitediff.cu
 * finite difference methods on a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160707
 */
 #include "finitediff.cuh"
 #include "errors.h"
 
 __constant__ float3 dev_cnus[4]; // i = x,y; j = 1,2,3,4
 
 void set1DerivativeParameters(const float hd_i[3] )
 {
	 float unscaled_cnu[4] { 1.f, 0.f, 0.f, 0.f };
	 
	 float3 *cnus = new float3[4];
 
	 for (int nu = 0; nu < 4; ++nu ) {
		 cnus[nu].x = unscaled_cnu[nu]*(1.f/(hd_i[0]*hd_i[0]) );
		 cnus[nu].y = unscaled_cnu[nu]*(1.f/(hd_i[1]*hd_i[1]) );
		 cnus[nu].z = unscaled_cnu[nu]*(1.f/(hd_i[2]*hd_i[2]) );
	 }
	 
	 HANDLE_ERROR(
		 cudaMemcpyToSymbol( dev_cnus, cnus, sizeof(float3)*4, 0, cudaMemcpyHostToDevice) 
	 ); // offset from start is 0
	 
	 delete[] cnus;		
 }
 
 void set2DerivativeParameters(const float hd_i[3]  )
 {
	 float unscaled_cnu[4] { 4.f/3.f, -1.f/12.f, 0.f, 0.f, };
	 
	 float3 *cnus = new float3[4];
 
	 for (int nu = 0; nu < 4; ++nu ) {
		 cnus[nu].x = unscaled_cnu[nu]*(1.f/(hd_i[0]*hd_i[0]) );
		 cnus[nu].y = unscaled_cnu[nu]*(1.f/(hd_i[1]*hd_i[1]) );
		 cnus[nu].z = unscaled_cnu[nu]*(1.f/(hd_i[2]*hd_i[2]) );
	 }
	 
	 HANDLE_ERROR(
		 cudaMemcpyToSymbol( dev_cnus, cnus, sizeof(float3)*4, 0, cudaMemcpyHostToDevice) 
	 ); // offset from start is 0
	 
	 delete[] cnus;		
 }
 
 
 
 void set3DerivativeParameters(const float hd_i[3] )
 {
	 float unscaled_cnu[4] { 3.f/2.f , -3.f/20.f, 1.f/90.f, 0.f };
	 
	 float3 *cnus = new float3[4];
 
	 for (int nu = 0; nu < 4; ++nu ) {
		 cnus[nu].x = unscaled_cnu[nu]*(1.f/(hd_i[0]*hd_i[0]) );
		 cnus[nu].y = unscaled_cnu[nu]*(1.f/(hd_i[1]*hd_i[1]) );
		 cnus[nu].z = unscaled_cnu[nu]*(1.f/(hd_i[2]*hd_i[2]) );
	 }
	 
	 HANDLE_ERROR(
		 cudaMemcpyToSymbol( dev_cnus, cnus, sizeof(float3)*4, 0, cudaMemcpyHostToDevice) 
	 ); // offset from start is 0
	 
	 delete[] cnus;		
	 
 }
 
 void set4DerivativeParameters(const float hd_i[3]  )
 {
	 float unscaled_cnu[4] { 8.f/5.f, -1.f/5.f, 8.f/315.f, -1.f/560.f  };
	 
	 float3 *cnus = new float3[4];
 
	 for (int nu = 0; nu < 4; ++nu ) {
		 cnus[nu].x = unscaled_cnu[nu]*(1.f/(hd_i[0]*hd_i[0]) );
		 cnus[nu].y = unscaled_cnu[nu]*(1.f/(hd_i[1]*hd_i[1]) );
		 cnus[nu].z = unscaled_cnu[nu]*(1.f/(hd_i[2]*hd_i[2]) );
	 }
	 
	 HANDLE_ERROR(
		 cudaMemcpyToSymbol( dev_cnus, cnus, sizeof(float3)*4, 0, cudaMemcpyHostToDevice) 
	 ); // offset from start is 0
	 
	 delete[] cnus;		
 }
 
 __device__ float dev_dirdblder1(float centerval, float stencil[1][2], float c_nus[4]) {
	 float tempvalue {0.f};
 
	 tempvalue += c_nus[0]*( stencil[0][1] + stencil[0][0] -2*centerval );
 
	 return tempvalue;
 }
 
 __device__ float dev_dirdblder2(float centerval, float stencil[2][2], float c_nus[4]) {
	 int NU {2};
	 float tempvalue {0.f};
 
	 for (int nu = 0; nu < NU; ++nu ) {
		 tempvalue += c_nus[nu]*( stencil[nu][1] + stencil[nu][0] - 2*centerval );
	 }
	 return tempvalue;
 }
 
 __device__ float dev_dirdblder3(float centerval, float stencil[3][2], float c_nus[4]) {
	 int NU {3};
	 float tempvalue {0.f};
		 
	 for (int nu = 0; nu < NU; ++nu ) {
		 tempvalue += c_nus[nu]*( stencil[nu][1] + stencil[nu][0] -2*centerval);
	 }
	 return tempvalue;
 }
 
 __device__ float dev_dirdblder4(float centerval, float stencil[4][2], float c_nus[4]) {
	 int NU {4};
	 float tempvalue {0.f};
 
	 for (int nu = 0; nu < NU; ++nu ) {
		 tempvalue += c_nus[nu]*( stencil[nu][1] + stencil[nu][0] -2*centerval);
	 }
	 return tempvalue;
 }
 
 // LAPLACIAN (LAP) in 1,2,3,4 stencils, central difference
 
 __device__ float dev_lap1( float centerval, float3 stencil[1][2]  ) {
	 float stencilx[1][2] { { stencil[0][0].x, stencil[0][1].x } };
	 float stencily[1][2] { { stencil[0][0].y, stencil[0][1].y } };
	 float stencilz[1][2] { { stencil[0][0].z, stencil[0][1].z } };
 
	 float c_nusx[4] { dev_cnus[0].x, dev_cnus[1].x, dev_cnus[2].x, dev_cnus[3].x };
	 float c_nusy[4] { dev_cnus[0].y, dev_cnus[1].y, dev_cnus[2].y, dev_cnus[3].y };
	 float c_nusz[4] { dev_cnus[0].z, dev_cnus[1].z, dev_cnus[2].z, dev_cnus[3].z };
	 
	 float lap_value { dev_dirdblder1(centerval, stencilx, c_nusx ) } ;
	 lap_value += dev_dirdblder1( centerval, stencily, c_nusy ) ;
	 lap_value += dev_dirdblder1( centerval, stencilz, c_nusz ) ;
 
	 return lap_value;
 }
 
 __device__ float dev_lap2(float centerval, float3 stencil[2][2]  ) {
	 float stencilx[2][2] { { stencil[0][0].x, stencil[0][1].x }, { stencil[1][0].x, stencil[1][1].x } };
	 float stencily[2][2] { { stencil[0][0].y, stencil[0][1].y }, { stencil[1][0].y, stencil[1][1].y } };
	 float stencilz[2][2] { { stencil[0][0].z, stencil[0][1].z }, { stencil[1][0].z, stencil[1][1].z } };
 
	 float c_nusx[4] { dev_cnus[0].x, dev_cnus[1].x, dev_cnus[2].x, dev_cnus[3].x };
	 float c_nusy[4] { dev_cnus[0].y, dev_cnus[1].y, dev_cnus[2].y, dev_cnus[3].y };
	 float c_nusz[4] { dev_cnus[0].z, dev_cnus[1].z, dev_cnus[2].z, dev_cnus[3].z };
	 
	 float lap_value { dev_dirdblder2( centerval, stencilx, c_nusx ) } ;
	 lap_value += dev_dirdblder2( centerval, stencily, c_nusy ) ;
	 lap_value += dev_dirdblder2( centerval, stencilz, c_nusz ) ;
 
	 return lap_value;
 }
 
 __device__ float dev_lap3( float centerval, float3 stencil[3][2]  ) {
	 float stencilx[3][2] { { stencil[0][0].x, stencil[0][1].x }, { stencil[1][0].x, stencil[1][1].x }, { stencil[2][0].x, stencil[2][1].x } };
	 float stencily[3][2] { { stencil[0][0].y, stencil[0][1].y }, { stencil[1][0].y, stencil[1][1].y }, { stencil[2][0].y, stencil[2][1].y } };
	 float stencilz[3][2] { { stencil[0][0].z, stencil[0][1].z }, { stencil[1][0].z, stencil[1][1].z }, { stencil[2][0].z, stencil[2][1].z } };
 
	 float c_nusx[4] { dev_cnus[0].x, dev_cnus[1].x, dev_cnus[2].x, dev_cnus[3].x };
	 float c_nusy[4] { dev_cnus[0].y, dev_cnus[1].y, dev_cnus[2].y, dev_cnus[3].y };
	 float c_nusz[4] { dev_cnus[0].z, dev_cnus[1].z, dev_cnus[2].z, dev_cnus[3].z };
	 
	 float lap_value { dev_dirdblder3(centerval, stencilx, c_nusx ) } ;
	 lap_value += dev_dirdblder3( centerval,stencily, c_nusy ) ;
	 lap_value += dev_dirdblder3( centerval,stencilz, c_nusz ) ;
 
	 return lap_value;
 }
 
 __device__ float dev_lap4( float centerval, float3 stencil[4][2]  ) {
	 float stencilx[4][2] { { stencil[0][0].x, stencil[0][1].x }, { stencil[1][0].x, stencil[1][1].x }, { stencil[2][0].x, stencil[2][1].x }, { stencil[3][0].x, stencil[3][1].x } };
	 float stencily[4][2] { { stencil[0][0].y, stencil[0][1].y }, { stencil[1][0].y, stencil[1][1].y }, { stencil[2][0].y, stencil[2][1].y }, { stencil[3][0].y, stencil[3][1].y } };
	 float stencilz[4][2] { { stencil[0][0].z, stencil[0][1].z }, { stencil[1][0].z, stencil[1][1].z }, { stencil[2][0].z, stencil[2][1].z }, { stencil[3][0].z, stencil[3][1].z } };
 
	 float c_nusx[4] { dev_cnus[0].x, dev_cnus[1].x, dev_cnus[2].x, dev_cnus[3].x };
	 float c_nusy[4] { dev_cnus[0].y, dev_cnus[1].y, dev_cnus[2].y, dev_cnus[3].y };
	 float c_nusz[4] { dev_cnus[0].z, dev_cnus[1].z, dev_cnus[2].z, dev_cnus[3].z };
	 
	 float lap_value { dev_dirdblder4( centerval, stencilx, c_nusx ) } ;
	 lap_value += dev_dirdblder4( centerval, stencily, c_nusy ) ;
	 lap_value += dev_dirdblder4( centerval, stencilz, c_nusz ) ;
 
	 return lap_value;
 }
 