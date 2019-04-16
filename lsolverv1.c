#include<stdio.h>
#include<math.h>
#include "HYPRE_struct_ls.h"
#include<string.h>

#define G 4.0	/*coeff of grid points*/
#define S -1.0	/*coeff of stencil points (neighbours)*/
#define B 0	/*value of RHS in Poisson's equation (if 0, then Laplacian)*/
#define Dl 0	
#define Dr 0
#define Dt 0
#define Db 0
#define Dtt 100
#define Drt 0
#define Nl 0
#define Nr 0
#define Nt 0
#define Nb 0
#define Ntt 0
#define Nrt 0

int main()
{
	HYPRE_StructGrid grid;
	HYPRE_StructStencil stencil;
	HYPRE_StructMatrix A;
	HYPRE_StructVector b;
	HYPRE_StructVector x;
	HYPRE_StructSolver solver;
	HYPRE_StructSolver precond;
	
	int n_pre, n_post;
	int rap, relax, skip, sym;
	int time_index;

	int num_iterations;
	double final_res_norm;

	char bottom[4]="Neu";
	char top[4]="Neu";
	char left[4]="Neu";
	char right[4]="Dir";
	char toptop[4]="Dir";
	char righttop[4]="Neu";
	
	int i,j;
	int n1=128;	/*length of box1*/
	int m1=64;	/*breadth of box1*/
	int n2=64;	/*length of box2*/
	int m2=64;	/*breadth of box2*/
	
	/*GRID*/
	{
		/*Create empty grid*/
		HYPRE_StructGridCreate(MPI_COMM_WORLD,2,&grid);
		
		{
			/*Set lower corner and upper corners of box1 and box2*/
			int ilower[2][2]={{0,0},{0,m1}};
			int iupper[2][2]={{n1-1,m1-1},{n2-1,m1+m2-1}};
			HYPRE_StructGridSetExtents(grid,ilower[0],iupper[0]);
			HYPRE_StructGridSetExtents(grid,ilower[1],iupper[1]);			
		}
		
		/*Assemble grid*/
		HYPRE_StructGridAssemble(grid);
	}
	
	
	/*STENCIL*/
	{
		/*Create empty stencil*/
		HYPRE_StructStencilCreate(2,5,&stencil);
		
		{
			/*Set stencil elements. Here we choose a 5-pt stencil*/
			int entry;
			int offset[5][2]={{0,0},{-1,0},{1,0},{0,-1},{0,1}};
			for(entry=0;entry<5;entry++)
			{
				HYPRE_StructStencilSetElement(stencil,entry,offset[entry]);
			}
		}
	
	}
	
		
	/*MATRIX*/
	{	
		/*Create empty matrix*/
		HYPRE_StructMatrixCreate(MPI_COMM_WORLD,grid,stencil,&A);
	

			/*Initialize matrix*/
			HYPRE_StructMatrixInitialize(A);
			int ilower[2][2]={{0,0},{0,m1}};
			int iupper[2][2]={{n1-1,m1-1},{n2,m1+m2-1}};
			
			int nvalues=5*(n1*m1+n2*m2);	/*Number of total points (stencil pts*grid points*)*/
			int nentry=5;	/*Number of stencil points*/
			
			double *values;
			values=calloc(nvalues,sizeof(double));
			
			int stencil_indices[5]={0,1,2,3,4}; /*Stencil indices in an array*/
			
			/*Add values to matrix*/
			
			for(i=0;i<nvalues;i+=nentry)
			{
				values[i]=G;
				for(j=1;j<nentry;j++)
				{
					values[i+j]=S;
				}
				
			}
			
			HYPRE_StructMatrixSetBoxValues(A,ilower[0],iupper[0],nentry,stencil_indices,values);
			HYPRE_StructMatrixSetBoxValues(A,ilower[1],iupper[1],nentry,stencil_indices,values);
			
			/*boundary conditions*/
			double *bvalues1;
			double *bvalues2;
			double *bvalues3;
			double *bvalues4;
			double *bvalues5;
			double *bvalues6;
			
			
			//~ bvalues1=calloc(n1*2,sizeof(double));
			//~ bvalues2=calloc((n1-n2)*2,sizeof(double));
			//~ bvalues3=calloc((m1+m2)*2,sizeof(double));
			//~ bvalues4=calloc(m1*2,sizeof(double));
			//~ bvalues5=calloc(n2*2,sizeof(double));
			//~ bvalues6=calloc((m2-m1)*2,sizeof(double));
			
			bvalues1=calloc((m1+m2)*2,sizeof(double));
			bvalues2=calloc((m1+m2)*2,sizeof(double));
			bvalues3=calloc((m1+m2)*2,sizeof(double));
			bvalues4=calloc((m1+m2)*2,sizeof(double));
			bvalues5=calloc((m1+m2)*2,sizeof(double));
			bvalues6=calloc((m1+m2)*2,sizeof(double));
			
			/*set stencils to emulate dirichlet or neumann boundaries. 3-(-1)-(-1)-(-1)-0 for neumann and 5-(-1)-(-1)-(-1)-0 for dirichlet*/
			for(i=0;i<(m1+m2)*2;i+=2)
			{
				/*bottom*/
				if(strcmp(bottom,"Neu")==0)
				{
					bvalues1[i]=3;
					bvalues1[i+1]=0;
				}
				if(strcmp(bottom,"Dir")==0)
				{
					bvalues1[i]=5;
					bvalues1[i+1]=0;
				}
				
				/*top1*/
				if(strcmp(top,"Neu")==0)
				{
					bvalues2[i]=3;
					bvalues2[i+1]=0;
				}
				if(strcmp(top,"Dir")==0)
				{
					bvalues2[i]=5;
					bvalues2[i+1]=0;
				}
				
				/*left*/
				if(strcmp(left,"Neu")==0)
				{
					bvalues3[i]=3;
					bvalues3[i+1]=0;
				}
				if(strcmp(left,"Dir")==0)
				{
					bvalues3[i]=5;
					bvalues3[i+1]=0;
				}
				
				/*right*/
				if(strcmp(right,"Neu")==0)
				{
					bvalues4[i]=3;
					bvalues4[i+1]=0;
				}
				if(strcmp(right,"Dir")==0)
				{
					bvalues4[i]=5;
					bvalues4[i+1]=0;
				}
				
				/*top-top*/
				if(strcmp(toptop,"Neu")==0)
				{
					bvalues5[i]=3;
					bvalues5[i+1]=0;
				}
				if(strcmp(toptop,"Dir")==0)
				{
					bvalues5[i]=5;
					bvalues5[i+1]=0;
				}
				
				/*right-top*/
				if(strcmp(righttop,"Neu")==0)
				{
					bvalues6[i]=3;
					bvalues6[i+1]=0;
				}
				if(strcmp(righttop,"Dir")==0)
				{
					bvalues6[i]=5;
					bvalues6[i+1]=0;
				}
		
			}

		
				
				/*bottom*/
				{
					int ilower[2]={0,0};
					int iupper[2]={n1-1,0};
					int stencil_indices[2]={0,3};
					HYPRE_StructMatrixSetBoxValues(A,ilower,iupper,2,stencil_indices,bvalues1);
					
				}
				/*top*/
				{
					int ilower[2]={n2,m1-1};
					int iupper[2]={n1-1,m1-1};
					int stencil_indices[2]={0,4};
					HYPRE_StructMatrixSetBoxValues(A,ilower,iupper,2,stencil_indices,bvalues2);
					
				}
				/*left*/
				{
					int ilower[2][2]={{0,0},{0,m1}};
					int iupper[2][2]={{0,m1-1},{0,m1+m2-1}};
					int stencil_indices[2]={0,1};
					HYPRE_StructMatrixSetBoxValues(A,ilower[0],iupper[0],2,stencil_indices,bvalues3);
					HYPRE_StructMatrixSetBoxValues(A,ilower[1],iupper[1],2,stencil_indices,bvalues3);
				}
				/*right*/
				{
					int ilower[2]={n1-1,0};
					int iupper[2]={n1-1,m1-1};
					int stencil_indices[2]={0,2};
					HYPRE_StructMatrixSetBoxValues(A,ilower,iupper,2,stencil_indices,bvalues4);					
				}
				/*top-top*/
				{
					int ilower[2]={0,m1+m2-1};
					int iupper[2]={n2-1,m1+m2-1};
					int stencil_indices[2]={0,4};
					HYPRE_StructMatrixSetBoxValues(A,ilower,iupper,2,stencil_indices,bvalues5);	
				}
				/*right-top*/
				{
					int ilower[2]={n2-1,m1};
					int iupper[2]={n2-1,m1+m2-1};
					int stencil_indices[2]={0,2};
					HYPRE_StructMatrixSetBoxValues(A,ilower,iupper,2,stencil_indices,bvalues6);
				}	

				/*set corner point stencils to 2-(-1)-(-1)-0-0 for neumann and 6-(-1)-(-1)-0-0 for dirichlet*/

				{
					
					double cvalues1[5];
					double cvalues2[5];
					double cvalues3[5];
					double cvalues4[5];
					double cvlaues5[5];
					
					/*bottom left*/
					{
						if(strcmp(bottom,"Neu")==0&&strcmp(left,"Neu")==0)
							cvalues1[0]=2;
						else if(strcmp(bottom,"Dir")==0&&strcmp(left,"Dir")==0)
							cvalues1[0]=6;
						else
							cvalues1[0]=4;
						cvalues1[1]=-1;
						cvalues1[2]=-1;
						cvalues1[3]=0;
						cvalues1[4]=0;
						int corner[2]={0,0};
						int stencil_indices[5]={0,2,4,1,3};
						
						HYPRE_StructMatrixSetBoxValues(A,corner,corner,nentry,stencil_indices,cvalues1);
					}
					/*bottom right*/
					{
						if(strcmp(bottom,"Neu")==0&&strcmp(right,"Neu")==0)
							cvalues2[0]=2;
						else if(strcmp(bottom,"Dir")==0&&strcmp(right,"Dir")==0)
							cvalues2[0]=6;
						else
							cvalues2[0]=4;
						cvalues2[1]=-1;
						cvalues2[2]=-1;
						cvalues2[3]=0;
						cvalues2[4]=0;
						int corner[2]={n1-1,0};
						int stencil_indices[5]={0,1,4,2,3};
						HYPRE_StructMatrixSetBoxValues(A,corner,corner,nentry,stencil_indices,cvalues2);
					}
					/*top-top left*/
					{
						if(strcmp(toptop,"Neu")==0&&strcmp(left,"Neu")==0)
							cvalues3[0]=2;
						else if(strcmp(toptop,"Dir")==0&&strcmp(left,"Dir")==0)
							cvalues3[0]=6;
						else
							cvalues3[0]=4;
						cvalues3[1]=-1;
						cvalues3[2]=-1;
						cvalues3[3]=0;
						cvalues3[4]=0;
						int corner[2]={0,m1+m2-1};
						int stencil_indices[5]={0,2,3,1,4};
						HYPRE_StructMatrixSetBoxValues(A,corner,corner,nentry,stencil_indices,cvalues3);
					}
					/*top right*/
					{
						if(strcmp(top,"Neu")==0&&strcmp(right,"Neu")==0)
							cvalues4[0]=2;
						else if(strcmp(top,"Dir")==0&&strcmp(right,"Dir")==0)
							cvalues4[0]=6;
						else
							cvalues4[0]=4;
						cvalues4[1]=-1;
						cvalues4[2]=-1;
						cvalues4[3]=0;
						cvalues4[4]=0;
						int corner[2]={n1-1,m1-1};
						int stencil_indices[5]={0,1,3,2,4};
						HYPRE_StructMatrixSetBoxValues(A,corner,corner,nentry,stencil_indices,cvalues4);
					}
					/*top-top right-top*/
					{
						if(strcmp(toptop,"Neu")==0&&strcmp(righttop,"Neu")==0)
							cvalues5[0]=2;
						else if(strcmp(toptop,"Dir")==0&&strcmp(righttop,"Dir")==0)
							cvalues5[0]=6;
						else
							cvalues5[0]=4;
						cvalues5[1]=-1;
						cvalues5[2]=-1;
						cvalues5[3]=0;
						cvalues5[4]=0;
						int corner[2]={n2-1,m1+m2-1};
						int stencil_indices[5]={0,1,3,2,4};
						HYPRE_StructMatrixSetBoxValues(A,corner,corner,nentry,stencil_indices,cvalues5);
					}
					
				}	
		/*Assemble matrix*/
		HYPRE_StructMatrixAssemble(A);
	}	

	/*VECTORS*/
	{
		/*Create empty vectors*/
		HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,&b);
		HYPRE_StructVectorCreate(MPI_COMM_WORLD,grid,&x);
		
		{
			/*Initialize vectors*/
			HYPRE_StructVectorInitialize(b);
			HYPRE_StructVectorInitialize(x);
			
			/*Add values to vectors*/
			int ilower[2][2]={{0,0},{0,m1}};
			int iupper[2][2]={{n1-1,m1-1},{n2,m1+m2-1}};
			
			/*for b vector*/
			double *values;
			values=calloc(n1*m1+n2*m2,sizeof(double));
			double h=(double)1/n1;
			for(i=0;i<(n1*m1+n2*m2);i++)
			{
				values[i]=-B*((double)h*h);
				
			}
			HYPRE_StructVectorSetBoxValues(b,ilower[0],iupper[0],values);
			HYPRE_StructVectorSetBoxValues(b,ilower[1],iupper[1],values);
			/*boundaries of b-vector*/
			double *bvalues1;
			double *bvalues2;
			double *bvalues3;
			double *bvalues4;
			double *bvalues5;
			double *bvalues6;
			
			//~ bvalues1=calloc(n1,sizeof(double));
			//~ bvalues2=calloc((n1-n2),sizeof(double));
			//~ bvalues3=calloc((m1+m2),sizeof(double));
			//~ bvalues4=calloc(m1,sizeof(double));
			//~ bvalues5=calloc(n2,sizeof(double));
			//~ bvalues6=calloc((m2-m1),sizeof(double));
			
			bvalues1=calloc((m1+m2)*2,sizeof(double));
			bvalues2=calloc((m1+m2)*2,sizeof(double));
			bvalues3=calloc((m1+m2)*2,sizeof(double));
			bvalues4=calloc((m1+m2)*2,sizeof(double));
			bvalues5=calloc((m1+m2)*2,sizeof(double));
			bvalues6=calloc((m1+m2)*2,sizeof(double));
			
			for(i=0;i<(m1+m2);i++)
			{
				
				/*bottom*/
				{
					if(strcmp(bottom,"Neu")==0)
						bvalues1[i]=(-B*((double)h*h))-(h*Nb);
					if(strcmp(bottom,"Dir")==0)
						bvalues1[i]=(-B*((double)h*h))+2*Db;	/*set boundary value on bottom boundary*/
					int ilower[2]={0,0};
					int iupper[2]={n1-1,0};
					HYPRE_StructVectorSetBoxValues(b,ilower,iupper,bvalues1);
					
					
				}
				/*top*/
				{
					if(strcmp(top,"Neu")==0)
						bvalues2[i]=(-B*((double)h*h))-(h*Nt);
					if(strcmp(top,"Dir")==0)
						bvalues2[i]=(-B*((double)h*h))+2*Dt;	/*set boundary value on upper boundary*/
					int ilower[2]={n2,m1-1};
					int iupper[2]={n1-1,m1-1};
					HYPRE_StructVectorSetBoxValues(b,ilower,iupper,bvalues2);
				}
				/*left*/
				{
					if(strcmp(left,"Neu")==0)
						bvalues3[i]=(-B*((double)h*h))-(h*Nl);
					if(strcmp(left,"Dir")==0)
						bvalues3[i]=(-B*((double)h*h))+2*Dl;	/*set boundary value on left boundary*/
					int ilower[2][2]={{0,0},{0,m1}};
					int iupper[2][2]={{0,m1-1},{0,m1+m2-1}};
					HYPRE_StructVectorSetBoxValues(b,ilower[0],iupper[0],bvalues3);
					HYPRE_StructVectorSetBoxValues(b,ilower[0],iupper[0],bvalues3);
				}
				/*right*/
				{
					if(strcmp(right,"Neu")==0)
						bvalues4[i]=(-B*((double)h*h))-(h*Nl);
					if(strcmp(right,"Dir")==0)
						bvalues4[i]=(-B*((double)h*h))+2*Dl;	/*set boundary value on right boundary*/
					int ilower[2]={n1-1,0};
					int iupper[2]={n1-1,m1-1};
					HYPRE_StructVectorSetBoxValues(b,ilower,iupper,bvalues4);
					
				}
				/*top-top*/
				{
					if(strcmp(toptop,"Neu")==0)
						bvalues5[i]=(-B*((double)h*h))-(h*Ntt);
					if(strcmp(toptop,"Dir")==0)
						bvalues5[i]=(-B*((double)h*h))+2*Dtt;	/*set boundary value on upper boundary*/
					int ilower[2]={0,m1+m2-1};
					int iupper[2]={n2-1,m1+m2-1};
					HYPRE_StructVectorSetBoxValues(b,ilower,iupper,bvalues5);
				}
				/*right-top*/
				{
					if(strcmp(righttop,"Neu")==0)
						bvalues6[i]=(-B*((double)h*h))-(h*Nrt);
					if(strcmp(righttop,"Dir")==0)
						bvalues6[i]=(-B*((double)h*h))+2*Drt;	/*set boundary value on right boundary*/
					int ilower[2]={n2-1,m1};
					int iupper[2]={n2-1,m1+m2-1};
					HYPRE_StructVectorSetBoxValues(b,ilower,iupper,bvalues6);
					
				}
				
					
			}
			
			
			/*corner points of b-vector*/
			double cvalues1[1];
			double cvalues2[1];
			double cvalues3[1];
			double cvalues4[1];
			double cvalues5[1];
			/*bottom left*/
			{
				if(strcmp(bottom,"Neu")==0&&strcmp(left,"Neu")==0)
					cvalues1[0]=(-B*((double)h*h))-(h*Nl)-(h*Nb);
				if(strcmp(bottom,"Dir")==0&&strcmp(left,"Dir")==0)
					cvalues1[0]=(-B*((double)h*h))+2*(Dl+Db);
				if(strcmp(bottom,"Neu")==0&&strcmp(left,"Dir")==0)
					cvalues1[0]=(-B*((double)h*h))+2*Dl-(h*Nb);
				if(strcmp(bottom,"Dir")==0&&strcmp(left,"Neu")==0)
					cvalues1[0]=(-B*((double)h*h))+2*Db-(h*Nl);
					
				int corner[2]={0,0};
				HYPRE_StructVectorSetBoxValues(b,corner,corner,cvalues1);
			}
			/*bottom right*/
			{
				if(strcmp(bottom,"Neu")==0&&strcmp(right,"Neu")==0)
					cvalues2[0]=(-B*((double)h*h))-(h*Nr)-(h*Nb);
				if(strcmp(bottom,"Dir")==0&&strcmp(right,"Dir")==0)
					cvalues2[0]=(-B*((double)h*h))+2*(Dr+Db);
				if(strcmp(bottom,"Neu")==0&&strcmp(right,"Dir")==0)
					cvalues2[0]=(-B*((double)h*h))+2*Dr-(h*Nb);
				if(strcmp(bottom,"Dir")==0&&strcmp(right,"Neu")==0)
					cvalues2[0]=(-B*((double)h*h))+2*Db-(h*Nr);
					
				int corner[2]={n1-1,0};
				HYPRE_StructVectorSetBoxValues(b,corner,corner,cvalues2);
			}
			/*top-top left*/
			{
				if(strcmp(top,"Neu")==0&&strcmp(left,"Neu")==0)
					cvalues3[0]=(-B*((double)h*h))-(h*Nl)-(h*Ntt);
				if(strcmp(top,"Dir")==0&&strcmp(left,"Dir")==0)
					cvalues3[0]=(-B*((double)h*h))+2*(Dl+Dtt);
				if(strcmp(top,"Neu")==0&&strcmp(left,"Dir")==0)
					cvalues3[0]=(-B*((double)h*h))+2*Dl-(h*Ntt);
				if(strcmp(top,"Dir")==0&&strcmp(left,"Neu")==0)
					cvalues3[0]=(-B*((double)h*h))+2*Dtt-(h*Nl);
				
				int corner[2]={0,m1+m2-1};
				HYPRE_StructVectorSetBoxValues(b,corner,corner,cvalues3);
			}
			/*top right*/
			{
				if(strcmp(top,"Neu")==0&&strcmp(right,"Neu")==0)
					cvalues4[0]=(-B*((double)h*h))-(h*Nr)-(h*Nt);
				if(strcmp(top,"Dir")==0&&strcmp(right,"Dir")==0)
					cvalues4[0]=(-B*((double)h*h))+2*(Dr+Dt);
				if(strcmp(top,"Neu")==0&&strcmp(right,"Dir")==0)
					cvalues4[0]=(-B*((double)h*h))+2*Dr-(h*Nt);
				if(strcmp(top,"Dir")==0&&strcmp(right,"Neu")==0)
					cvalues4[0]=(-B*((double)h*h))+2*Dt-(h*Nr);
				
				int corner[2]={n1-1,m1-1};
				HYPRE_StructVectorSetBoxValues(b,corner,corner,cvalues4);
			}
			/*top-top right-top*/
			{
				if(strcmp(top,"Neu")==0&&strcmp(right,"Neu")==0)
					cvalues4[0]=(-B*((double)h*h))-(h*Nrt)-(h*Ntt);
				if(strcmp(top,"Dir")==0&&strcmp(right,"Dir")==0)
					cvalues4[0]=(-B*((double)h*h))+2*(Drt+Dtt);
				if(strcmp(top,"Neu")==0&&strcmp(right,"Dir")==0)
					cvalues4[0]=(-B*((double)h*h))+2*Drt-(h*Ntt);
				if(strcmp(top,"Dir")==0&&strcmp(right,"Neu")==0)
					cvalues4[0]=(-B*((double)h*h))+2*Dtt-(h*Nrt);
				
				int corner[2]={n2-1,m1+m2-1};
				HYPRE_StructVectorSetBoxValues(b,corner,corner,cvalues5);
			}
			/*for x vector*/
			for(i=0;i<(n1*m1+n2*m2);i++)
			{
				values[i]=0.0;
				HYPRE_StructVectorSetBoxValues(x,ilower[0],iupper[0],values);
				HYPRE_StructVectorSetBoxValues(x,ilower[1],iupper[1],values);
			}
			
		}
		HYPRE_StructVectorAssemble(b);
		HYPRE_StructVectorAssemble(x);
	}
	
	/*SOLVER*/
	{
		/*Solver 1-PCG*/
		/*Create solver*/
		HYPRE_StructPCGCreate(MPI_COMM_WORLD,&solver);
		
		/*Set parameters for solver*/
		HYPRE_StructPCGSetTol(solver,1e-08);
		HYPRE_StructPCGSetPrintLevel(solver,2);
		HYPRE_StructPCGSetMaxIter(solver, 50);
		
	      /* Use symmetric SMG as preconditioner */
	      HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
	      HYPRE_StructSMGSetMaxIter(precond, 1);
	      HYPRE_StructSMGSetTol(precond, 0.0);
	      HYPRE_StructSMGSetZeroGuess(precond);
	      HYPRE_StructSMGSetNumPreRelax(precond, 1);
	      HYPRE_StructSMGSetNumPostRelax(precond, 1);

		
		/*Setup the solver and solve for x*/
		HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve,
		HYPRE_StructSMGSetup, precond);
		HYPRE_StructPCGSetup(solver,A,b,x);
		HYPRE_StructPCGSolve(solver,A,b,x);
		

	}

	/*PRINT RESULT TO FILE*/
	{
		
		int Dim = 2;
		int nentry=5;
		int i, j,k;
		double h = 1.0/n1;
		int stencil_indices[5]={0,1,2,3,4};
		int ilower[2][2]={{0,0},{0,m1}};
		int iupper[2][2]={{n1-1,m1-1},{n2-1,m1+m2-1}};
		double *values;
		values = (double*) calloc(n1*m1+n2*m2,sizeof(double));
		
		{
			FILE *file=fopen("Lmesh1","w");	
			HYPRE_StructVectorGetBoxValues(x, ilower[0], iupper[0], values);
				for (j = 0; j < m1; j++)
				{
					for (i = 0; i < n1; i++)
					{
						fprintf(file, "%.5e\t", values[i + j*n1]);	
					}
					fprintf(file,"\n");
				}
			fflush(file);
			fclose(file);
		}
		
		{
			FILE *file=fopen("Lmesh2","w");
			HYPRE_StructVectorGetBoxValues(x, ilower[1], iupper[1], values);
				for (j = 0; j < m2; j++)
				{
					for (i = 0; i < n2; i++)
					{
						fprintf(file, "%.5e\t", values[i + j*n2]);	
					}
					fprintf(file,"\n");
				}
			fflush(file);
			fclose(file);
			
			
		double *mvalues;
		mvalues = (double*) calloc((n1*m1+n2*m2)*5, sizeof(double));
			
		HYPRE_StructMatrixGetBoxValues(A,ilower[0],iupper[0],nentry,stencil_indices, mvalues);
			for (j = 0; j < (n1*m1*5); j++)
			{
				
				fprintf(file, "%.5e\n", mvalues[j]);	
				
			}
			fprintf(file,"\n\n\n");
		HYPRE_StructMatrixGetBoxValues(A,ilower[1],iupper[1],nentry,stencil_indices, mvalues);
			for (j = 0; j < (n2*m2*5); j++)
			{
				
				fprintf(file, "%.5e\n", mvalues[j]);	
				
			}	
			
		}	
			
		
	}

	/* Free memory */
	HYPRE_StructGridDestroy(grid);
	HYPRE_StructStencilDestroy(stencil);
	HYPRE_StructMatrixDestroy(A);
	HYPRE_StructVectorDestroy(b);
	HYPRE_StructVectorDestroy(x);
	HYPRE_StructPCGDestroy(solver);
	HYPRE_StructSMGDestroy(precond);

	/* Finalize MPI */
	MPI_Finalize();

	//~ for(i=0;i<n*n*5;i++)
	//~ {
		//~ printf("%g\n",debug[i]);
	//~ }

	return (0);
}
		
			