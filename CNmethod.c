#include<stdio.h> 
#include<math.h>
#include<stdlib.h>
#define f(x) (2.-1.5*(x)+sin(M_PI*(x))) //initial heat distribution
#define T 1000 //time steps
void lusolver(double *e2, double *a2, double *c2, double *b2, double *x2, int m);

/*
Solución a la ecuación de calor:
u_{t}=1.44u_{xx}
con condición inicial 
u(x,0)=f(x)
y condiciones de Dirichlet u(0,t)=2., u(L,t)=0.5 para todo t>0.
A partir del Método Crank-Nicolson
*/

int main()
{
    int i, j;
    FILE *archivo=fopen("heateqnCN.txt","w");
    double L=1., dx=0.005, dt=0.001, D=1.44;
    double T_l=2., T_r=0.5; //Boundary conditions
    int sp=ceil(L/dx);
    double x[sp+1],t[T+1];
    double u[sp+1];

    double x1[sp-1];
    double a1[sp-1], e1[sp-1], c1[sp-1];
    double b1[sp-1];

    double alpha=D*dt/(dx*dx);

    //Construimos la condición inicial
    x[0]=0.;
    t[0]=0.;
    u[0]=T_l;
    fprintf(archivo, "%lf %lf %lf\n",x[0], t[0], u[0]);
    for ( i = 1; i < sp ; i++)
    {
     x[i]=i*dx;
     u[i]=f(x[i]);
     fprintf(archivo, "%lf %lf %lf\n",x[i],t[0], u[i]);
    }
    x[sp]=L;
    u[sp]=T_r;
    fprintf(archivo, "%lf %lf %lf\n",x[sp],t[0],u[sp]);
    fprintf(archivo, "\n" );

    //Resolvemos para los instantes posteriores

    //Construimos los vectores que caracterizan a la matriz A:
    //Vector e superior
    for(i=0;i<sp-2;i++)
    {
         e1[i]=-alpha;
    }
    e1[sp-2]=0.;
    //Vector a diagonal 
    for(i=0;i<sp-1;i++)
    {
         a1[i]=2*(1+alpha);
    }
    //Vector c inferior
    c1[0]=0.;
    for(i=1;i<sp-1;i++)
    {
        c1[i]=-alpha;
    }

    for ( j = 1; j <= T; j++)
    {
        x[0]=0.;
        t[j]=j*dt;
        u[0]=T_l;
        fprintf(archivo, "%lf %lf %lf\n", x[0], t[j], u[0]);
        x[sp]=L;
        u[sp]=T_r;
        //Formamos el vector b de la ecuación Ax=b
        b1[0]=alpha*T_l+alpha*u[2]+2*(1-alpha)*u[1]+alpha*u[0];
        for(i=1;i<sp-2;i++)
        {
            b1[i]=alpha*u[i+2]+2*(1-alpha)*u[i+1]+alpha*u[i];
        }
        b1[sp-2]=alpha*T_r+alpha*u[sp]+2*(1-alpha)*u[sp-1]+alpha*u[sp-2];
        lusolver(e1, a1, c1, b1, x1, sp-1);
        for(i=0; i<sp-1; i++)
        {
            u[i+1]=x1[i];
        }
        for(i=1; i<sp; i++)
        {
            fprintf(archivo, "%lf %lf %lf\n", x[i], t[j], u[i]);
        }
        fprintf(archivo, "%lf %lf %lf\n", x[sp], t[j], u[sp]);
        fprintf(archivo, "\n" );
    }
    fclose(archivo);
}

/*Función que resuelve Ax=b para A una matriz triagonal.
La función recibe los vectores e, a y c que caracterizan a la 
matriz triagonal, b un vector y x la variable a resolver.
Devuelve los valores de la solución en el vector x. 
*/
void lusolver(double *e2, double *a2, double *c2, double *b2, double *x2, int m)
{
    int i, j;
    double a[m], e[m], c[m]; //de los vectores e,c no ocupamos todas sus entradas, pero sirve para indexar adecuadamente
    double w[m], u[m], b[m];
    double x[m], y[m];
    for(i=0; i<m; i++)
    {
        a[i]=*a2;
        b[i]=*b2;
        a2++;
        b2++;
    }
    for(i=0; i<m-1; i++)
    {
        e[i]=*e2;
        e2++;
    }
    for(i=1; i<m; i++)
    {
        c2++;
        c[i]=*c2;
    }
    //Factorización LU de la matriz m
    //Calculamos las entradas w asociadas a L y u asociadas a U:
    w[0] = a[0];
    for (i = 1; i < m; i++)
    {
        u[i - 1] = e[i - 1] / w[i - 1];
        w[i] = a[i] - u[i - 1] * c[i]; // Es importante ver cómo está indexado c.
    }


    //Resolvemos para y en Ly=b
    y[0] = b[0] / w[0];
    for (i = 1; i < m; i++)
    {
        y[i] = (b[i] - c[i] * y[i - 1]) / (w[i]);
    }
    //Resolvemos para x en Ux=y
    x[m - 1] = y[m - 1];
    for (i = m - 1; i > 0; i--)
    {
        x[i - 1] = y[i - 1] - u[i - 1] * x[i];
    }
    for(i=0; i<m ; i++)
    {
        *x2=x[i];
        x2++;
    }
}