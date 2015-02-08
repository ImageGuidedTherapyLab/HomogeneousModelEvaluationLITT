/*
 * Example Matlab cuda kernel interface.
 */


__device__
void pointSource(double rVar, double r1, double r2, double wPerf, double cblood, double kCond, double mueff, double u0, double ua, double Power, double *temperature )
{

   double pi = 3.141592653589793;
//*temperature  =ua-(exp(-mueff*rVar)*powf(mueff,2)*Power)/(4*kCond*powf(mueff,2)*pi*rVar-4*cblood*pi*rVar*wPerf)+(exp(-mueff*(r1+rVar)+(r1+r2-rVar)*sqrt((cblood*wPerf)/kCond))*(exp(r1*(mueff+sqrt((cblood*wPerf)/kCond)))*powf(mueff,2)*Power*powf(r2,2)*(1+mueff*rVar);

*temperature=ua-(exp(-mueff*rVar)*powf(mueff,2)*Power)/(4*kCond*powf(mueff,2)*pi*rVar-4*cblood*pi*rVar*wPerf)+(exp(-mueff*(r1+rVar)+(r1+r2-rVar)*sqrt((cblood*wPerf)/kCond))*(exp(r1*(mueff+sqrt((cblood*wPerf)/kCond)))*powf(mueff,2)*Power*powf(r2,2)*(1+mueff*rVar)-exp(mueff*rVar+r2*sqrt((cblood*wPerf)/kCond))*powf(mueff,2)*Power*powf(rVar,2)*(-1+r2*sqrt((cblood*wPerf)/kCond))+4*exp(mueff*(r1+rVar)+r2*sqrt((cblood*wPerf)/kCond))*pi*r1*powf(rVar,2)*(u0-ua)*(-kCond*powf(mueff,2)+cblood*wPerf)*(-1+r2*sqrt((cblood*wPerf)/kCond))))/(4*pi*powf(rVar,3)*(-kCond*powf(mueff,2)+cblood*wPerf)*(exp(2*r2*sqrt((cblood*wPerf)/kCond))*(-1+r2*sqrt((cblood*wPerf)/kCond))+exp(2*r1*sqrt((cblood*wPerf)/kCond))*(1+r2*sqrt((cblood*wPerf)/kCond))))+(exp(-mueff*(r1+rVar)+(2*r1+rVar)*sqrt((cblood*wPerf)/kCond))*(-exp(mueff*r1+r2*sqrt((cblood*wPerf)/kCond))*powf(mueff,2)*Power*powf(r2,2)*(1+mueff*rVar)-exp(mueff*rVar+r1*sqrt((cblood*wPerf)/kCond))*powf(mueff,2)*Power*powf(rVar,2)*(1+r2*sqrt((cblood*wPerf)/kCond))-4*exp(mueff*(r1+rVar)+r1*sqrt((cblood*wPerf)/kCond))*pi*r1*powf(rVar,2)*(u0-ua)*(kCond*powf(mueff,2)-cblood*wPerf)*(1+r2*sqrt((cblood*wPerf)/kCond))))/(4*pi*powf(r1,2)*powf(rVar,3)*(-kCond*powf(mueff,2)+cblood*wPerf)*(exp(2*r2*sqrt((cblood*wPerf)/kCond))*(-1+r2*sqrt((cblood*wPerf)/kCond))+exp(2*r1*sqrt((cblood*wPerf)/kCond))*(1+r2*sqrt((cblood*wPerf)/kCond))));

//	*temperature = ua-(exp(-mueff*rVar)*mueff*mueff*Power)/(4*kCond*mueff*mueff*pi*rVar-4*cblood*pi*rVar*wPerf)+(exp(-mueff*(r1+rVar)+(r1+r2-rVar)*sqrt((cblood*wPerf)/kCond))*(exp(r1*(mueff+sqrt((cblood*wPerf)/kCond)))*mueff*mueff*Power*r2*r2*(1+mueff*rVar)-exp(mueff*rVar+r2*sqrt((cblood*wPerf)/kCond))*mueff*mueff*Power*rVar*rVar*(-1+r2*sqrt((cblood*wPerf)/kCond))+4*exp(mueff*(r1+rVar)+r2*sqrt((cblood*wPerf)/kCond))*pi*r1*rVar*rVar*(u0-ua)*(-kCond*mueff*mueff+cblood*wPerf)*(-1+r2*sqrt((cblood*wPerf)/kCond))))/(4*pi*rVar*rVar*rVar*(-kCond*mueff*mueff+cblood*wPerf)*(exp(2*r2*sqrt((cblood*wPerf)/kCond))*(-1+r2*sqrt((cblood*wPerf)/kCond))+exp(2*r1*sqrt((cblood*wPerf)/kCond))*(1+r2*sqrt((cblood*wPerf)/kCond))))+(exp(-mueff*(r1+rVar)+(2*r1+rVar)*sqrt((cblood*wPerf)/kCond))*(-exp(mueff*r1+r2*sqrt((cblood*wPerf)/kCond))*mueff*mueff*Power*r2*r2*(1+mueff*rVar)-exp(mueff*rVar+r1*sqrt((cblood*wPerf)/kCond))*mueff*mueff*Power*rVar*rVar*(1+r2*sqrt((cblood*wPerf)/kCond))-4*exp(mueff*(r1+rVar)+r1*sqrt((cblood*wPerf)/kCond))*pi*r1*rVar*rVar*(u0-ua)*(kCond*mueff*mueff-cblood*wPerf)*(1+r2*sqrt((cblood*wPerf)/kCond))))/(4*pi*r1*r1*rVar*rVar*rVar*(-kCond*mueff*mueff+cblood*wPerf)*(exp(2*r2*sqrt((cblood*wPerf)/kCond))*(-1+r2*sqrt((cblood*wPerf)/kCond))+exp(2*r1*sqrt((cblood*wPerf)/kCond))*(1+r2*sqrt((cblood*wPerf)/kCond))));

//   *temperature = ua+(P*PI_Var*(mueff*mueff)*exp(-mueff*r)*(1.0/4.0))/(r*(w-k*(mueff*mueff)))-(exp(-R1*mueff-R2*mueff)*exp(r*sqrt(w/k))*(P*PI_Var*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R2*mueff)-P*PI_Var*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R1*mueff)-P*PI_Var*R2*(mueff*mueff*mueff)*exp(R2*sqrt(w/k))*exp(R1*mueff)-R1*u0*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*ua*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*k*(mueff*mueff)*u0*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*k*(mueff*mueff)*ua*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+P*PI_Var*R2*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R2*mueff)*sqrt(w/k)-R1*R2*u0*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*ua*w*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*k*(mueff*mueff)*u0*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0-R1*R2*k*(mueff*mueff)*ua*exp(R1*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0)*(1.0/4.0))/(r*(w-k*(mueff*mueff))*(exp(R1*sqrt(w/k)*2.0)-exp(R2*sqrt(w/k)*2.0)+R2*exp(R1*sqrt(w/k)*2.0)*sqrt(w/k)+R2*exp(R2*sqrt(w/k)*2.0)*sqrt(w/k)))-(exp(R1*sqrt(w/k))*exp(R2*sqrt(w/k))*exp(-r*sqrt(w/k))*exp(-R1*mueff)*exp(-R2*mueff)*(P*PI_Var*(mueff*mueff)*exp(R1*sqrt(w/k))*exp(R1*mueff)-P*PI_Var*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R2*mueff)+P*PI_Var*R2*(mueff*mueff*mueff)*exp(R1*sqrt(w/k))*exp(R1*mueff)+R1*u0*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*ua*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0-R1*k*(mueff*mueff)*u0*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+R1*k*(mueff*mueff)*ua*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*4.0+P*PI_Var*R2*(mueff*mueff)*exp(R2*sqrt(w/k))*exp(R2*mueff)*sqrt(w/k)-R1*R2*u0*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*ua*w*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0+R1*R2*k*(mueff*mueff)*u0*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0-R1*R2*k*(mueff*mueff)*ua*exp(R2*sqrt(w/k))*exp(R1*mueff)*exp(R2*mueff)*sqrt(w/k)*4.0)*(1.0/4.0))/(r*(w-k*(mueff*mueff))*(exp(R1*sqrt(w/k)*2.0)-exp(R2*sqrt(w/k)*2.0)+R2*exp(R1*sqrt(w/k)*2.0)*sqrt(w/k)+R2*exp(R2*sqrt(w/k)*2.0)*sqrt(w/k)));
}
__device__
void DebugWrite(int idx,int idmat,double rad,double omega, double conduction, double mueff,double temp)
{
   printf("%d %d %12.5e %12.5e %12.5e %12.5e %12.5e\n",idx,idmat,rad,omega,conduction,mueff,temp);
   //int j,k;

   //for (j=0;j<n;j++) {
   //   for (k=0;k<n+1;k++) {
   //      printf("%d %d %12.5e %12.5e ",k,j,a[k][j].real(),a[k][j].imag());
   //   }
   //   printf(" | %d  %12.5e %12.5e \n",j,x[j].real(),x[j].imag());
   //}
   //printf("\n");
}

/*
 * Device code
 */
__global__ 
void steadyStatePennesLaser(
         int const NTissue,
         const    int* MaterialID,
         const double* Perfusion,
         const double* ThermalConduction,
         const double* EffectiveAttenuation,
         double const innerRadius,
         double const outerRadius,
         int const NSource,
         double const Power,
         const double* SourceXloc,
         const double* SourceYloc,
         const double* SourceZloc,
         double const InitialTemperature,
         double const ArterialTemperature,
         double const SpecificHeatBlood,
	 double const SpacingX,
	 double const SpacingY,
	 double const SpacingZ,
         int const NpixelX,
         int const NpixelY,
         int const NpixelZ,
         double* d_TemperatureArray)
{

//     double SpacingX=0.00078;
    /*
      grid stride loop design pattern, 1-d grid
      http://devblogs.nvidia.com/parallelforall/cuda-pro-tip-write-flexible-kernels-grid-stride-loops/
         - By using a loop, you can support any problem size even if it exceeds the largest grid size your CUDA device supports. Moreover, you can limit the number of blocks you use to tune performance.
    */
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x; 
         idx < NpixelX * NpixelY * NpixelZ;
         idx += blockDim.x * gridDim.x) 
      {
        // compute indices
        int index = idx; // use dummy variable
        int kkk = index/(NpixelX*NpixelY); 
        index -= kkk*NpixelX*NpixelY; 
        
        int jjj = index/NpixelX; 
        index -= jjj*NpixelX; 
        
        int iii = index/1;

        /* get material parameters */
        int const idmaterial =  MaterialID[idx];
        double omega      = Perfusion[idmaterial];
        double conduction = ThermalConduction[idmaterial];
        double mueff      = EffectiveAttenuation[idmaterial];
//	printf("%d",mueff);
        // linear superpostion of temperature sources
        double temperature = 0.0;
        for (int lll=0;lll<NSource;lll++) 
          {
//           double radiusSQ = (iii * SpacingX + 0.13281 - SourceXloc[lll])*(iii * SpacingX + 0.13281 - SourceXloc[lll])
//                           + (jjj * SpacingY + 0.10547 - SourceYloc[lll])*(jjj * SpacingY + 0.10547 - SourceYloc[lll])
//                           + (kkk * SpacingZ + 0.06000 - SourceZloc[lll])*(kkk * SpacingZ + 0.06000- SourceZloc[lll]);

	   double radiusSQ=powf(iii*SpacingX-SourceXloc[lll],2)
			  +powf(jjj*SpacingY-SourceYloc[lll],2)
			  +powf(kkk*SpacingZ-SourceZloc[lll],2);//SourceXloc[0]*SourceXloc[0];
           double radius   = sqrt(radiusSQ);

           // call GF code 
	   double sourcetemperature;
           pointSource(radius, innerRadius, outerRadius, omega , SpecificHeatBlood, conduction , mueff, InitialTemperature, ArterialTemperature, Power , &sourcetemperature);

	   if (radius <= innerRadius && NSource ==1)
		{
                    sourcetemperature = InitialTemperature;
		}
           if (radius <= innerRadius && NSource == 10)
		{
                    sourcetemperature = InitialTemperature+55;
		}
           if (radius <= innerRadius && NSource > 1)
		{
                   sourcetemperature = InitialTemperature;
		}
           // DebugWrite(idx,idmaterial,radius,omega,conduction,mueff,sourcetemperature);
           // superposition
	   if (idmaterial==0)
	 	{
		   temperature=0;
		}
	   else
		{
                   temperature = temperature + sourcetemperature/((double)NSource); 		
		}	 
          }
        // store temperature in array
        d_TemperatureArray[idx] = temperature;
      }
}


