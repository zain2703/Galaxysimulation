#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include<pthread.h>
#include <omp.h>

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

typedef struct galaxy{
	double posX;
	double posY;
	double mass;
	double velX;
	double velY;
	double brightness;
}gal,*link;

typedef struct tree
{
    int flag;        //indicate node is empty or not 
    double massX;      //Center of mass X
    double massY;      //  Center of mass Y
    double m;          //Total mass for a box containing many particals; for a partical it's its mass
    double centX;         //X center of box
    double centY;         //Y center of box
    double width;      //width of the box
    struct tree* ne;   // Northeast 
    struct tree* nw;   // Northwest
    struct tree* sw;   // Southeast
    struct tree* se;   // Southwest
} treenode;

struct Targ{
	 int T_N,T_n_threads,T_i;
	double  T_G;
	double T_thetamax;
	double T_delta_t;
	//double *T_xcomp;
	//double *T_ycomp;
	treenode* T_tree;
	gal *T_obj,*T_obj2;	
};


void tree_insert(treenode **node, double x, double y, double m);
void tree_split(treenode **node);
void cal_BH(double x, double y, double theta_max, treenode *node, double* Xcomp, double* Ycomp);
void cal_P(gal *obj,gal * obj2,int N,int nsteps,double theta_max, double delta_t,int n_threads);

int main(int argc, const char *argv[])
{
  if (argc != 8)
	  {
		printf("Enter 6 input args\n");
		return -1;
	  }
	  const int N = atoi(argv[1]);
	  const char* fileName1 = argv[2];
	  const int nsteps = atoi(argv[3]);
	  const double delta_t = atof(argv[4]);
	  const double theta_max = atof(argv[5]);
	  int graphics = atoi(argv[6]);
	  int n_threads = atoi(argv[7]);
  /*Read the file and save in 2 objects */
      gal *obj= (gal*) malloc (N*sizeof(gal));
	  gal *obj2= (gal*) malloc (N*sizeof(gal));
	  FILE *f=fopen(fileName1,"rb");
  if(fread(obj,N*sizeof(gal),1,f)==0)
	  {
		  printf("Error opneing the file\n");
		  exit -1;	  
	  }
	  fseek(f, 0L, SEEK_SET);
  if(fread(obj2,N*sizeof(gal),1,f)==0)
	  {
		  printf("Error opneing the second obj2 file\n");
		  exit -1;	  
	  }
	  printf("sizeof structure =%ld \n",sizeof(gal));  
     double t1 = get_wall_seconds();
     cal_P(obj,obj2,N,nsteps,theta_max,delta_t,n_threads);
	 double t2 = get_wall_seconds();
	 printf("Total time taken=%f\n",t2-t1);
  /*write the  file*/
 FILE *f2=fopen("result.gal", "w");
	  if(fwrite(obj,N*sizeof(gal),1,f2)==0)
	  {
		  printf("Error opneing the file for write\n");
		  exit -1;	  
	  }
	  if (graphics==1)
		  {
			  printf("Graphics file not added'\n");
		  }
	  fclose(f2);
  return 0;
}

void *T_1(void *arg)
{      
        //const double G = 100.0 / N;
		int N,n_threads,chunksize,i;
		double accX,  accY;
	    struct Targ *mydata;
		double theta_max, delta_t,G; //posX, posY;
		double Xcomp,Ycomp;
		struct galaxy *obj2, *obj;
		struct Tree *tree;
		////Type casting the structure
		mydata=(struct Targ*) arg;
		////Storing data from strucutre in variables
		i=mydata->T_i;
		N=mydata->T_N;
		n_threads=mydata->T_n_threads;
		theta_max=mydata->T_thetamax;
		G=mydata->T_G;
		delta_t=mydata->T_delta_t;
		obj=mydata->T_obj;
		obj2=mydata->T_obj2;
		tree=mydata->T_tree;
		/* Xcomp=mydata->T_xcomp;
		Ycomp=mydata->T_ycomp; */
		//dividing threads into equal parts
		chunksize=N/n_threads;
		int start= chunksize*i;
		int end = start + chunksize;
		//for ood case
		if( i == n_threads-1 )
        {end += (N % n_threads);}
		
		//printf("Thread=%d, start=%d, end=%d \n",i, start,(end-1));
		//printf("In thread function Sorting of elements from structure into variable\n");
		// calculation for force
		 for(int i=start;i<end;i++)
		{
			//printf("In thread function before BH\n");
			Xcomp = Ycomp = 0.0;
			cal_BH(obj[i].posX, obj[i].posY, theta_max, tree, &Xcomp, &Ycomp);
			//locks on mutex
	//		pthread_mutex_lock (&mutexsum);
			accX = -G * (Xcomp);
			accY = -G * (Ycomp);
			obj2[i] .velX += delta_t * accX;
			obj2[i] .velY += delta_t * accY;
			obj2[i].posX += delta_t * obj2[i] .velX;
			obj2[i].posY += delta_t * obj2[i] .velY;
	//		pthread_mutex_unlock (&mutexsum);
		} 
}

void cal_P(gal *obj,gal * obj2,int N,int nsteps,double theta_max, double delta_t, int n_threads)
{
	//// Initiazing  structure for threads
	  struct Targ thread_data_array[n_threads];
	  pthread_t threads[n_threads];
	  pthread_attr_t attr;
	  pthread_attr_init(&attr);
	  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	   // mutex initlizaion 
      //pthread_mutex_init(&mutexsum, NULL);
	  //
	   void *status;
		const double G = 100.0 / N; 
	  //  double *Xcomp = (double *)malloc(sizeof(double));
	  //  double *Ycomp = (double *)malloc(sizeof(double));
	    double accX,  accY;
		//printf("print numer of threads=%d /n",n_threads);
   for(int t = 0; t < nsteps; t++)
  {
    treenode* tree = NULL;
    for(int i = 0; i < N; i++) {
			  tree_insert(&tree, obj[i].posX, obj[i].posY, obj[i].mass);
			}
    for (int i = 0; i < n_threads; i++)
			{
			  thread_data_array[i].T_i=i;	
			  thread_data_array[i].T_N=N;
			  thread_data_array[i].T_n_threads=n_threads;
			  thread_data_array[i].T_thetamax=theta_max;
			  thread_data_array[i].T_delta_t=delta_t;
			  thread_data_array[i].T_G=G;
		//	  thread_data_array[i].T_xcomp=Xcomp;
		//     thread_data_array[i].T_ycomp=Ycomp;
			  thread_data_array[i].T_tree=tree;
		      thread_data_array[i].T_obj=obj;
		      thread_data_array[i].T_obj2=obj2;
			  
		      pthread_create(&threads[i], &attr, T_1, (void *)&thread_data_array[i]);
			}
			//printf("After thread ends in cal_P\n");
		for(int t=0; t<n_threads; t++) 
		{
			pthread_join(threads[t], &status);
	    }
		    //printf("join ends in cal_P\n");
    //after each step, update buf array*/
     for (int i = 0; i < N; i++)
			{
			  obj[i].posX = obj2[i].posX;
			  obj[i].posY = obj2[i].posY;
			  obj[i].velX = obj2[i].velX;
			  obj[i].velY = obj2[i].velY;
			}
    free(tree);
  }
  //free(Xcomp);
  //free(Ycomp); 
}

void tree_insert(treenode **node, double x, double y, double m)
{
    
    if (!*node)
    {
        //when tree is emtpy
        *node = (treenode *)malloc(sizeof(treenode));
        (*node)->flag = 1;
        (*node)->massX = x;
        (*node)->massY = y;            
        (*node)->m = m;
        (*node)->width = 1;           // width to 1 to make equal parts 
        (*node)->centX = 0.5 * ((*node)->width);
        (*node)->centY = 0.5 * ((*node)->width);
        (*node)->ne = NULL;
        (*node)->nw = NULL;
        (*node)->sw = NULL;
        (*node)->se = NULL;
        return;
    }

    if (*node)
    {
        if ((*node)->flag == 1)
        {
            
            tree_split(node);
        }

        if((*node)->flag == 0) {
           
            if((x >= (*node)->centX) && (y >= (*node)->centY)) {
                    if(!(*node)->ne) {
                        treenode *tree = (treenode *)malloc(sizeof(treenode));
                        tree->flag = 1;
                        tree->massX = x;
                        tree->massY = y;
                        tree->m = m;
                        tree->ne = NULL;
                        tree->nw = NULL;
                        tree->sw = NULL;
                        tree->se = NULL;
                        tree->width = (*node)->width * 0.5;
                        tree->centX =(*node)->centX + tree->width * 0.5;
                        tree->centY =(*node)->centY + tree->width * 0.5;
                        (*node)->ne = tree;
                    }
                    else {
                        tree_insert(&((*node)->ne), x, y, m);
                    }
                }

                else if ((x >= (*node)->centX) && (y < (*node)->centY)) {
                    if(!(*node)->se) {
                        treenode *tree = (treenode *)malloc(sizeof(treenode));
                        tree->flag = 1;
                        tree->massX = x;
                        tree->massY = y;
                        tree->m = m;
                        tree->ne = NULL;
                        tree->nw = NULL;
                        tree->sw = NULL;
                        tree->se = NULL;
                        tree->width = (*node)->width * 0.5;
                        tree->centX =(*node)->centX + tree->width * 0.5;
                        tree->centY =(*node)->centY - tree->width * 0.5;
                        (*node)->se = tree;
                    }
                    else {
                        tree_insert(&((*node)->se), x, y, m);
                    }
                }
            
            else if  ((x < (*node)->centX) && (y >= (*node)->centY)){
                    if(!(*node)->nw) {
                        treenode *tree = (treenode *)malloc(sizeof(treenode));
                        tree->flag = 1;
                        tree->massX = x;
                        tree->massY = y;
                        tree->m = m;
                        tree->ne = NULL;
                        tree->nw = NULL;
                        tree->sw = NULL;
                        tree->se = NULL;
                        tree->width = (*node)->width * 0.5;
                        tree->centX =(*node)->centX - tree->width * 0.5;
                        tree->centY =(*node)->centY + tree->width * 0.5;
                        (*node)->nw = tree;
                    }
                    else {
                        tree_insert(&((*node)->nw), x, y, m);
                    }
                }
                else {
                    if(!(*node)->sw) {
                        treenode *tree = (treenode *)malloc(sizeof(treenode));
                        tree->flag = 1;
                        tree->massX = x;
                        tree->massY = y;
                        tree->m = m;
                        tree->ne = NULL;
                        tree->nw = NULL;
                        tree->sw = NULL;
                        tree->se = NULL;
                        tree->width = (*node)->width * 0.5;
                        tree->centX =(*node)->centX - tree->width * 0.5;
                        tree->centY =(*node)->centY - tree->width * 0.5;
                        (*node)->sw = tree;
                    }
                    else {
                        tree_insert(&((*node)->sw), x, y, m);
                    }
                } 
             //Storing everything in head node
            (*node)->massX = ((*node)->massX * (*node)->m + x * m)/((*node)->m+m);
            (*node)->massY = ((*node)->massY * (*node)->m + y * m)/((*node)->m+m);
            (*node)->m += m;
		}	
    }
 }

void tree_split(treenode **node) {
    
    treenode *tree = (treenode *)malloc(sizeof(treenode));
    tree->flag = 1;
    tree->massX = (*node)->massX;
    tree->massY = (*node)->massY;
    tree->m = (*node)->m;
    tree->width = 0.5*((*node)->width);
    tree->ne = NULL;
    tree->nw = NULL;
    tree->sw = NULL;
    tree->se = NULL;

  
		if((tree->massX >= (*node)->centX)&&(tree->massY >= (*node)->centY)) {
				(*node)->ne = tree;
				tree->centX = (*node)->centX + 0.5*(tree->width);   
				tree->centY = (*node)->centY + 0.5*(tree->width);
				(*node)->flag = 0;
			}
		 else if ((tree->massX >= (*node)->centX)&&(tree->massY < (*node)->centY)){
				(*node)->se = tree;
				tree->centX = (*node)->centX + 0.5*(tree->width);
				tree->centY = (*node)->centY - 0.5*(tree->width);
				(*node)->flag = 0;
			}
		else if ((tree->massX < (*node)->centX)&&(tree->massY >= (*node)->centY)){
				(*node)->nw = tree;
				tree->centX = (*node)->centX - 0.5*(tree->width);
				tree->centY = (*node)->centY + 0.5*(tree->width);
				(*node)->flag = 0;
			}
        else {
            (*node)->sw = tree;
            tree->centX = (*node)->centX - 0.5*(tree->width);
            tree->centY = (*node)->centY - 0.5*(tree->width);
            (*node)->flag = 0;
        }
}

void cal_BH(double x, double y, double theta_max, treenode *node, double* Xcomp, double* Ycomp){
      const double e = 1e-3;
	  double x_j,y_j,m_j, theta,r,p; 
	  x_j = node->massX;
	  y_j = node->massY;
	  m_j = node->m;	  
	  r = sqrt((x-x_j)*(x-x_j)+(y-y_j)*(y-y_j));
	  if(node->flag == 1) {
		//flag 1 means a partical exist
		if(r != 0) {
		  p = m_j / ((r + e) * (r + e) * (r + e));
		  *Xcomp += p*(x-x_j);
		  *Ycomp += p*(y-y_j);
					}
		}
		else if(node->flag == 0) {
					theta = node->width / r;
					if(theta <= theta_max) {
						  p = m_j / ((r + e) * (r + e) * (r + e));
						  *Xcomp += p*(x-x_j);
						  *Ycomp += p*(y-y_j);
					}
					else {
					  //DFS to calculate others
					  if(node->ne)
						cal_BH(x, y, theta_max, node->ne, Xcomp, Ycomp);
					  if(node->nw)
						cal_BH(x, y, theta_max, node->nw, Xcomp, Ycomp);
					  if(node->sw)
						cal_BH(x, y, theta_max, node->sw, Xcomp, Ycomp);
					  if(node->se)
						cal_BH(x, y, theta_max, node->se, Xcomp, Ycomp);
					}
		}	
	}