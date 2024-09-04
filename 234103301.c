//Abhishek Gautam
//234103301


#include<stdio.h>
#include<math.h>
// #include<conio.h>

int main()
{
    FILE *fp;
    FILE *fp1;
    FILE *fp2;
    
    fp=fopen("Output.dat","w");
    fp1=fopen("Centerline_u_velocity.dat","w");
    fp2=fopen("Centerline_v_velocity.dat","w");
    
    int i,j,m,n,p,Re;
    int iteration ;
    float delx,dely,x,y,B;
    //B = beta
    //m = no of grid points along x direction
    //n = no of grid points along y direction
    printf("\n Enter the value of m : ");
    scanf("%d",&m);
    printf("\n Enter the value of n : ");
    scanf("%d",&n);
    printf("\nEnter the Reynolds Number : ");
    scanf("%d",&Re);

    p=(m-2)*(n-2);   //total no. of interior points
    delx=1.0/(m-1);
    dely=1.0/(n-1);

    B=(delx/dely);

    float psi[m][n],psi_old[m][n],omg[m][n],omg_old[m][n],u[m][n],v[m][n];

    float error_psi=0.0,error_omg=0.0;

    iteration=0;

    //boundary condition initialization
    for(j=0; j<n; j++) {
        for(i=0; i<m; i++) {

            v[i][j]=0.0;
            psi[i][j]=0.0;
            if(j==(n-1)) {
                u[i][j]=1.0;
            } else
                u[i][j]=0.0;
        }
    }

    for(j=0; j<n; j++) {
        for(i=0; i<m; i++) {
            if(j==0)
                omg[i][j]=(2.0/dely*dely)*(psi[i][j]-psi[i][j+1]);
            else if(i==0)
                omg[i][j]=(2.0/delx*delx)*(psi[i][j]-psi[i+1][j]);
            else if(i==(m-1))
                omg[i][j]=(2.0/delx*delx)*(psi[i][j]-psi[i-1][j]);
            else if(j==(n-1))
                omg[i][j]=(2.0/dely*dely)*(psi[i][j]-psi[i][j-1])-((2.0/dely)*u[i][j]);
            else
                omg[i][j]=0.0;

        }
    }

    //GAUSS-SIEDEL METHOD

    do {
        for(j=0; j<n; j++) {
            for(i=0; i<m; i++) {
                psi_old[i][j]=psi[i][j];
                omg_old[i][j]=omg[i][j];
            }
        }
        //solving for stream function
        for(j=1; j<(n-1); j++) {
            for(i=1; i<(m-1); i++) {
                psi[i][j]=(1.0/(2*(1.0+B*B)))*(psi[i+1][j]+psi[i-1][j]+B*B*(psi[i][j+1]+psi[i][j-1])+delx*delx*omg[i][j]);
            }
        }

        //solving for vorticity Equation
        for(j=1; j<(n-1); j++) {
            for(i=1; i<(m-1); i++) {
                omg[i][j] = (0.5/(1.0+pow(B,2)))*((1.0-(psi[i][j+1]-psi[i][j-1])*((B*Re)/4.0))*omg[i+1][j]+ (1.0+(psi[i][j+1]-psi[i][j-1])*((B*Re)/4.0))*omg[i-1][j]
                                                  + (1.0+(psi[i+1][j]-psi[i-1][j])*(Re/(4.0*B)))*(pow(B,2)*omg[i][j+1])+ (1.0-(psi[i+1][j]-psi[i-1][j])*(Re/(4.0*B)))*(pow(B,2)*omg[i][j-1]));;
            }
        }
        // UPDATE VORTICITY AT BOUNDARIES

        for(j=0; j<n; j++) {
            omg[0][j]=(-2*(psi[1][j]-psi[0][j]))/(delx*delx);
            omg[m-1][j]=(-2*(psi[m-2][j]-psi[m-1][j]))/(delx*delx);
        }
        for(i=0; i<m; i++) {
            omg[i][0]=(-2*(psi[i][1]-psi[i][0]))/(dely*dely);
            omg[i][n-1]=(-2*(psi[i][n-2]-psi[i][n-1]+dely))/(dely*dely);
        }
        //ERROR CALCULATION for psi and omega
        error_psi=0.0;
        error_omg=0.0;
        for(j=1; j<(n-1); j++) {
            for(i=1; i<(m-1); i++) {

                error_psi=error_psi+pow((psi[i][j]-psi_old[i][j]),2.0);
                error_omg=error_omg+pow((omg[i][j]-omg_old[i][j]),2.0);
            }
        }
        error_psi=sqrt(error_psi/p);
        error_omg=sqrt(error_omg/p);

        printf("iteration=%d\t",iteration);
        printf("error_psi=%.9lf\terror_omg=%.9lf\n",error_psi,error_omg);
        iteration++;
    } while(error_psi>pow(10,-6) || error_omg>pow(10,-6));

    //UPDATING VELOCITIES

    for (j=1; j<(n-1); j++) {
        for(i=1; i<(m-1); i++) {
            u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2*dely);
            v[i][j]=(psi[i+1][j]-psi[i-1][j])/(-2.0*delx);
        }
    }
    fprintf(fp,"ZONE I=%d, J=%d\n",m,n);
    for(i = 0; i < m; i++) {
        x=i*delx;
        for(j = 0; j < n; j++) {
            y=j*dely;
            fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\n",x,y,u[i][j],v[i][j],psi[i][j],omg[i][j]);
        }
    }
    for(i=0; i<m; i++)
        fprintf(fp1,"%f \t %f \n",u[n/2][i],i*dely);
    for(j=0; j<n; j++)
        fprintf(fp2,"%lf \t %lf \n",j*delx,v[j][m/2]);
    fclose(fp);
    fclose(fp1);
    fclose(fp2);
    return 0;
}