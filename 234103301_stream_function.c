#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int main()
{
    FILE *omega; // declaration of pointers
    FILE *psi1;
    FILE *u1;
    FILE *v1;
    FILE *vv;

// declaration of variables
	int i,j,m,n,k=0;
	double delx,dely,e_omega,e_psi,e_omega1=1.0,e_psi1=1.0,T,Re;
	printf("enter the values of m,n,Re");
	scanf("%d%d%lf",&m,&n,&Re);  //taking input
	delx=1.0/(m-1);
	dely=1.0/(n-1);
	double B;
	double psi[m][n],psiold[m][n],w[m][n],wold[m][n],u[m][n],v[m][n];
    B=delx/dely;
    T=B*B;
//opening the file
	omega=fopen("omega1.plt","w");
	psi1=fopen("psi_1.plt","w");
	u1=fopen("u_1.plt","w");
	v1=fopen("v_1.plt","w");
	vv=fopen("vv2.plt","w");

//initializing the interior points
	for(j=1;j<n-1;j++)
    {
        for(i=1;i<m-1;i++)
        {
          w[i][j]=0;
          psi[i][j]=0;
          u[i][j]=0;
          v[i][j]=0;
        }
    }
//boundary conditions
	for(i=0;i<m;i++)
    {
        psi[i][0]=0;
        psi[i][n-1]=0;
        u[i][0]=0;
        u[i][n-1]=1;
        v[i][0]=0;
        v[i][n-1]=0;
    }

    for(j=0;j<n;j++)
    {
        psi[0][j]=0;
        psi[m-1][j]=0;
        u[0][j]=0;
        u[m-1][j]=0;
        v[0][j]=0;
        v[m-1][j]=0;
    }
//boundary conditions for vorticity
    for(j=0;j<n;j++)
    {
        for(i=0;i<m;i++)
        {
            w[i][0]=0;

            w[i][n-1]=-((2/dely));

            w[0][j]=0;

            w[m-1][j]=0;
        }
    }
 //iteration loop
    do
    {
        e_omega=0;
        e_psi=0;
//storing the values to calculate error for vorticity and stream function
        for(j=0;j<n;j++)
   {
        for(i=0;i<m;i++)
        {
            psiold[i][j]=psi[i][j];
            wold[i][j]=w[i][j];
        }
   }
//calculation for stream function
         for(j=1;j<n-1;j++)
   {
        for(i=1;i<m-1;i++)
        {
            psi[i][j]=(psi[i+1][j]+psi[i-1][j]+T*(psi[i][j+1]+psi[i][j-1])+w[i][j]*(pow(delx,2)))/(2*(1+T));
        }
   }
//calculation for vorticity
     for(j=1;j<n-1;j++)
   {
        for(i=1;i<m-1;i++)
        {
          w[i][j]=((1-(psi[i][j+1]-psi[i][j-1])*((B*Re)/4))*w[i+1][j]+(1+(psi[i][j+1]-psi[i][j-1])*((B*Re)/4))*w[i-1][j]+(1+(psi[i+1][j]-psi[i-1][j])*((Re)/(B*4)))*T*w[i][j+1]+(1-(psi[i+1][j]-psi[i-1][j])*((Re)/(B*4)))*T*w[i][j-1])/(2*(1+T));
        }
   }

 //updating the value of omega at the boundary
  for(j=0;j<n;j++)
    {
        for(i=0;i<m;i++)
        {
            w[i][0]=2*(psi[i][0]-psi[i][1])/(pow(dely,2));
            w[i][n-1]=((2/pow(dely,2))*(psi[i][n-1]-psi[i][n-2]))-((2/dely));
            w[0][j]=-2*(psi[1][j]-psi[0][j])/(pow(delx,2));
            w[m-1][j]=-2*(psi[m-2][j]-psi[m-1][j])/(pow(delx,2));
        }
    }

   for(j=1;j<n-1;j++)
	{
		for(i=1;i<m-1;i++)
		{
        e_omega+=pow((w[i][j]-wold[i][j]),2);//computing sum for omega
        e_psi+=pow((psi[i][j]-psiold[i][j]),2); //computing sum for psi
       }
   }

    e_omega1=sqrt(e_omega/((m-2)*(n-2))); // calculating error for vorticity
    e_psi1=sqrt(e_psi/((m-2)*(n-2)));// calculating error for stream function

    k++;  //incrementing the no of iterations

}while(e_omega1>pow(10,-6)||e_psi1>pow(10,-6));  //comparing if the error is within the given limit for the maximum one

printf("%d\t%lf\t%lf",k,e_omega,e_psi);

//updating the values of velocity for both u and v
for(j=1;j<n-1;j++)
	{
		for(i=1;i<m-1;i++)
		{
        u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2*dely);
        v[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2*delx);
       }
   }
//storing the results in a file to get the contours
for(i=0;i<m;i++)
{
    fprintf(v1,"%lf\t%lf\n",i*delx,v[i][50]);//velocity for horizontal centerline

    for (j=0;j<n;j++)
    {//
       fprintf(psi1,"%lf\t%lf\t%lf\n",i*delx,j*dely,psi[i][j]);
       fprintf(omega,"%lf\t%lf\t%lf\n",i*delx,j*dely,w[i][j]);
       fprintf(vv,"%lf\t%lf\t%lf\t%lf\n",i*delx,j*delx,u[i][j],v[i][j]);
    }
}

 for (j=0;j<n;j++)
    {
      fprintf(u1,"%lf\t%lf\n",u[50][j],j*dely); //velocity for vertical centerline
    }

    fclose(omega);
    fclose(psi1);
    fclose(u1);
    fclose(v1);
    fclose(vv);

return 0;
}
